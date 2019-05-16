from SimuMoleScripts.simulation_main_script import Simulation
from formtools.wizard.views import CookieWizardView
from .models import UploadForm
from django.http import HttpResponseRedirect, JsonResponse, HttpResponse
from django.urls import reverse
from django.contrib import messages
from django.shortcuts import render
from django.core.files.storage import FileSystemStorage
from django.conf import settings
import os
import threading
import fnmatch

################################
#   Home, News, Contact, About
################################

def home(request):
    some_dict = {}
    return render(request, 'home.html', some_dict)


def news(request):
    some_dict = {}
    return render(request, 'news.html', some_dict)


def contact(request):
    some_dict = {}
    return render(request, 'contact.html', some_dict)


def about(request):
    some_dict = {}
    return render(request, 'about.html', some_dict)


################################
#   Simulation Result
################################

def update_simulation_status(request):
    f = open(os.path.join(settings.MEDIA_ROOT, 'files', 'simulation_status.txt'), 'r')
    simulation_status = f.read()
    f.close()

    f = open(os.path.join(settings.MEDIA_ROOT, 'files', 'simulation_status_during_run.txt'), 'r')
    simulation_status_during_run_lines = f.readlines()
    simulation_status_during_run = '' if len(simulation_status_during_run_lines) == 0 \
        else simulation_status_during_run_lines[-1]
    f.close()

    context = {'simulation_status': simulation_status, 'simulation_status_during_run': simulation_status_during_run}
    return JsonResponse(context)


################################
#   Create Simulation
################################

class SimulationWizard(CookieWizardView):
    template_name = 'create_simulation.html'

    # file_storage:
    file_storage = FileSystemStorage(location=os.path.join(settings.MEDIA_ROOT, 'files'))

    def delete_temp_files(self):
        path_to_temp_dir = self.file_storage.base_location
        listdir = os.listdir(path_to_temp_dir)
        for file in listdir:
            os.remove(os.path.join(path_to_temp_dir, file))

    def clean_form_dict(self, dict_):
        """
        get dictionary that represent the form data, and return clean dictionary data (without contradictions)
        """
        clean_dict = {}
        first_pdb_type, first_pdb_id, first_pdb_file, second_pdb_type, second_pdb_id, second_pdb_file = '', '', '', '', '', ''
        x1, y1, z1, x2, y2, z2 = '0', '0', '0', '0', '0', '0'
        degXY_1, degYZ_1, degXY_2, degYZ_2 = '0', '0', '0', '0'

        num_of_proteins = dict_.get('num_of_proteins')

        first_pdb_type = dict_.get('first_pdb_type')
        if first_pdb_type == 'by_id':
            first_pdb_id = dict_.get('first_pdb_id')
            first_pdb_file = ''
        elif first_pdb_type == 'by_file':
            first_pdb_id = ''
            first_pdb_file = dict_.get('first_pdb_file')

        if num_of_proteins == '2':
            second_pdb_type = dict_.get('second_pdb_type')
            if second_pdb_type == 'by_id':
                second_pdb_id = dict_.get('second_pdb_id')
                second_pdb_file = ''
            elif first_pdb_type == 'by_file':
                second_pdb_id = ''
                second_pdb_file = dict_.get('second_pdb_file')
            x2, y2, z2 = dict_.get('x2', 0), dict_.get('y2', 0), dict_.get('z2', 0)
            degXY_2, degYZ_2 = dict_.get('degXY_2', 0), dict_.get('degYZ_2', 0)

        x1, y1, z1 = dict_.get('x1', 0), dict_.get('y1', 0), dict_.get('z1', 0)
        degXY_1, degYZ_1 = dict_.get('degXY_1', 0), dict_.get('degYZ_1', 0)
        temperature = dict_.get('temperature', '')
        production_steps = dict_.get('production_steps', '')

        clean_dict['num_of_proteins'] = num_of_proteins
        clean_dict['first_pdb_type'] = first_pdb_type
        clean_dict['first_pdb_id'] = first_pdb_id
        clean_dict['first_pdb_file'] = first_pdb_file
        clean_dict['second_pdb_type'] = second_pdb_type
        clean_dict['second_pdb_id'] = second_pdb_id
        clean_dict['second_pdb_file'] = second_pdb_file
        clean_dict['x1'] = x1
        clean_dict['y1'] = y1
        clean_dict['z1'] = z1
        clean_dict['x2'] = x2
        clean_dict['y2'] = y2
        clean_dict['z2'] = z2
        clean_dict['degXY_1'] = degXY_1
        clean_dict['degYZ_1'] = degYZ_1
        clean_dict['degXY_2'] = degXY_2
        clean_dict['degYZ_2'] = degYZ_2
        clean_dict['temperature'] = temperature
        clean_dict['production_steps'] = production_steps

        return clean_dict

    def create_simulation_thread(self, form_dict):
        s = Simulation(form_dict['num_of_proteins'],
                       form_dict['first_pdb_type'], form_dict['first_pdb_id'],
                       form_dict['second_pdb_type'], form_dict['second_pdb_id'],
                       form_dict['x1'], form_dict['y1'], form_dict['z1'],
                       form_dict['x2'], form_dict['y2'], form_dict['z2'],
                       form_dict['degXY_1'], form_dict['degYZ_1'],
                       form_dict['degXY_2'], form_dict['degYZ_2'],
                       form_dict['temperature'], form_dict['production_steps'])
        s.create_simulation()

    def done(self, form_list, **kwargs):
        """
        override "done": this function is called when the form is submitted
        """
        # Processing the input parameters:
        form_data = [form.cleaned_data for form in form_list]
        form_dict = {k: v for d in form_data for k, v in d.items()}  # convert list of dictionaries to one dictionary
        form_dict = self.clean_form_dict(form_dict)

        # Create a new thread responsible for creating the simulation:
        t = threading.Thread(target=self.create_simulation_thread, args=(form_dict,))
        t.setDaemon(True)
        t.start()

        # Initialize the status file:
        with open(os.path.join(settings.MEDIA_ROOT, 'files', 'simulation_status.txt'), "w+") as f:
            f.write("Processing your parameters...")
        with open(os.path.join(settings.MEDIA_ROOT, 'files', 'simulation_status_during_run.txt'), "w+") as f:
            f.write("")

        # Render 'create_simulation_result.html' without waiting until the simulation is complete:
        return render(self.request, 'create_simulation_result.html', {'form_data': form_dict})

    def get_form_initial(self, step):
        """
        override "get_form_initial" for getting data from previous steps
        """
        initial_data = self.initial_dict.get(step, {})  # initial data of current step

        # SimulationForm0_LoadPdb
        step_0_prev_data = self.storage.get_step_data('0')
        step_0_prev_data = {} if step_0_prev_data is None \
            else {'num_of_proteins': step_0_prev_data.get('0-num_of_proteins'),
                  'first_pdb_type': step_0_prev_data.get('0-first_pdb_type'),
                  'first_pdb_id': step_0_prev_data.get('0-first_pdb_id'),
                  'first_pdb_file': step_0_prev_data.get('0-first_pdb_file'),
                  'second_pdb_type': step_0_prev_data.get('0-second_pdb_type'),
                  'second_pdb_id': step_0_prev_data.get('0-second_pdb_id'),
                  'second_pdb_file': step_0_prev_data.get('0-second_pdb_file'), }

        # SimulationForm1_DetermineRelativePosition
        step_1_prev_data = self.storage.get_step_data('1')
        step_1_prev_data = {} if step_1_prev_data is None \
            else {'x1': step_1_prev_data.get('1-x1'), 'x2': step_1_prev_data.get('1-x2'),
                  'y1': step_1_prev_data.get('1-y1'), 'y2': step_1_prev_data.get('1-y2'),
                  'z1': step_1_prev_data.get('1-z1'), 'z2': step_1_prev_data.get('1-z2')}

        # SimulationForm2_SimulationParameters
        step_2_prev_data = self.storage.get_step_data('2')
        step_2_prev_data = {} if step_2_prev_data is None \
            else {'temperature': step_2_prev_data.get('2-temperature'),
                  'production_steps': step_2_prev_data.get('2-production_steps')}

        update_data = {**step_0_prev_data, **step_1_prev_data, **step_2_prev_data, **initial_data}
        return self.initial_dict.get(step, update_data)


def show_form1(wizard: CookieWizardView):
    """
    if 'num_of_proteins'==1: return FALSE, and then navigate to step 2 (simulation parameters)
    else, if 'num_of_proteins'==2: return TRUE, and then navigate to step 1 (determine relative position)
    """
    cleaned_data = wizard.get_cleaned_data_for_step('0') or {}
    return cleaned_data.get('num_of_proteins') == '2'


################################
#   File Upload
################################

def file_upload(request):
    if request.method == "POST":
        file = UploadForm(request.POST, request.FILES)
        if file.is_valid():
            file_name, file_extension = os.path.splitext(request.FILES['file'].name)
            file_extension = file_extension.lower()
            # allowing only pdb files
            if file_extension == '.pdb':
                messages.success(request, 'File Uploaded Successfully')
                file.save()
                return HttpResponseRedirect(reverse('file_upload'))
            else:  # Should never be called, since we added FileExtensionValidator on the Upload model.
                messages.error(request, 'Upload failed: file extension has to be \'pdb\'.')

    else:
        file = UploadForm()
    return render(request, 'file_upload.html', {'form': file})

################################
#   Display Video
################################


def display_video(request, video_name=None):
    if(video_name is None):
        return HttpResponse("Please Enter Video Name")
    # video_url = os.path.join(settings.MEDIA_ROOT, 'videos', 'video_1.mp4')
    video_url = settings.MEDIA_URL + 'videos/'
    return render(request, "video_display.html", {"url": video_url})

