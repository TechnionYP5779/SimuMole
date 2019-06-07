from SimuMoleScripts.simulation_main_script import Simulation
from SimuMoleScripts.uploaded_simulation import create_animations
from SimuMoleScripts.emailSender import send_email

from .forms import UploadFiles

from formtools.wizard.views import CookieWizardView
from django.http import JsonResponse
from django.shortcuts import render
from django.core.files.storage import FileSystemStorage
from django.conf import settings
import os
import threading
import zipfile

temp = 'media/files/'  # path to temp folder


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

    video_url = settings.MEDIA_URL + 'videos/'
    context = {'simulation_status': simulation_status, 'simulation_status_during_run': simulation_status_during_run,
               'video_path': video_url}
    return JsonResponse(context)


def download__create_zip(num_of_proteins, include_pdb_file, include_dcd_file, include_animations_files, previous_page):
    files = []
    files_names = []

    if include_pdb_file:
        if previous_page == 'create_simulation':
            file_name = ''
            if num_of_proteins == '1':
                file_name = '_1___movement.pdb'
            if num_of_proteins == '2':
                file_name = 'both___1___movement__2___movement.pdb'
            files.append(os.path.join(settings.MEDIA_ROOT, 'files', file_name))
        if previous_page == 'upload_files':
            files.append(os.path.join(settings.MEDIA_ROOT, 'files', 'file_upload_pdb.pdb'))
        files_names.append('pdb.pdb')

    if include_dcd_file:
        if previous_page == 'create_simulation':
            files.append(os.path.join(settings.MEDIA_ROOT, 'files', 'trajectory.dcd'))
        if previous_page == 'upload_files':
            files.append(os.path.join(settings.MEDIA_ROOT, 'files', 'file_upload_dcd.dcd'))
        files_names.append('dcd.dcd')

    if include_animations_files:
        angels = [(0, 0, 0),
                  (90, 0, 0), (180, 0, 0), (270, 0, 0),  # X axis
                  (0, 90, 0), (0, 180, 0), (0, 270, 0),  # Y axis
                  (0, 0, 90), (0, 0, 180), (0, 0, 270), ]  # Z axis
        for i, (x, y, z) in zip(range(0, 10), angels):
            video_name = 'video_{}.mp4'.format(str(i))
            files.append(os.path.join(settings.MEDIA_ROOT, 'videos', video_name))
            video_name_at_download = 'video_{}__X{}_Y{}_Z{}.mp4'.format(str(i), str(x), str(y), str(z))
            files_names.append(video_name_at_download)

    print(files_names)
    zip_file = zipfile.ZipFile(os.path.join(settings.MEDIA_ROOT, 'files', "SimuMole_output.zip"), "w")
    for file, file_name in zip(files, files_names):
        zip_file.write(file, file_name)
    zip_file.close()


def download__zip(request):
    num_of_proteins = request.GET.get('num_of_proteins')
    include_pdb_file = (request.GET.get('pdb_file') == 'true')
    include_dcd_file = (request.GET.get('dcd_file') == 'true')
    include_animations_files = (request.GET.get('animation_files') == 'true')
    previous_page = request.GET.get('previous_page')

    download__create_zip(num_of_proteins, include_pdb_file, include_dcd_file, include_animations_files, previous_page)

    return JsonResponse({})


def download__email(request):
    num_of_proteins = request.GET.get('num_of_proteins')
    include_pdb_file = (request.GET.get('pdb_file') == 'true')
    include_dcd_file = (request.GET.get('dcd_file') == 'true')
    include_animations_files = (request.GET.get('animation_files') == 'true')
    previous_page = request.GET.get('previous_page')

    email = request.GET.get('email')

    response = {'email_success': 'true'}
    try:
        download__create_zip(num_of_proteins, include_pdb_file, include_dcd_file, include_animations_files,
                             previous_page)
        send_email(email, "SimuMole_output.zip")
    except Exception as e:
        response = {'email_success': 'false'}

    return JsonResponse(response)


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

        temperature_scale = dict_.get('temperature_scale', '')
        temperature = dict_.get('temperature', '')
        time_step_number = dict_.get('time_step_number', '')

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
        clean_dict['temperature_scale'] = temperature_scale
        clean_dict['temperature'] = temperature
        clean_dict['time_step_number'] = time_step_number

        return clean_dict

    def create_simulation_thread(self, form_dict):
        s = Simulation(form_dict['num_of_proteins'],
                       form_dict['first_pdb_type'], form_dict['first_pdb_id'],
                       form_dict['second_pdb_type'], form_dict['second_pdb_id'],
                       form_dict['x1'], form_dict['y1'], form_dict['z1'],
                       form_dict['x2'], form_dict['y2'], form_dict['z2'],
                       form_dict['degXY_1'], form_dict['degYZ_1'],
                       form_dict['degXY_2'], form_dict['degYZ_2'],
                       form_dict['temperature_scale'], form_dict['temperature'],
                       form_dict['time_step_number'])
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
        return render(self.request, 'create_simulation_result.html',
                      {'form_data': form_dict, 'num_of_proteins': form_dict['num_of_proteins'],
                       'previous_page': 'create_simulation',
                       'video_path': settings.MEDIA_URL + 'videos/'})  # todo: change "video_path"

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
                  'z1': step_1_prev_data.get('1-z1'), 'z2': step_1_prev_data.get('1-z2'),
                  'degXY_1': step_1_prev_data.get('1-degXY_1'), 'degYZ_1': step_1_prev_data.get('1-degYZ_1'),
                  'degXY_2': step_1_prev_data.get('1-degXY_2'), 'degYZ_2': step_1_prev_data.get('1-degYZ_2')}

        # SimulationForm2_SimulationParameters
        step_2_prev_data = self.storage.get_step_data('2')
        step_2_prev_data = {} if step_2_prev_data is None \
            else {'temperature_scale': step_2_prev_data.get('2-temperature_scale'),
                  'temperature': step_2_prev_data.get('2-temperature'),
                  'time_step_number': step_2_prev_data.get('2-time_step_number')}

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
#   Upload PDB & DCD
################################

def upload_files(request):
    if request.method == 'POST':
        form = UploadFiles(request.POST, request.FILES)
        if form.is_valid():
            # Create a new thread responsible for creating the animations:
            t = threading.Thread(target=create_animations, args=())
            t.setDaemon(True)
            t.start()

            # Initialize the status file:
            with open(os.path.join(settings.MEDIA_ROOT, 'files', 'simulation_status.txt'), "w+") as f:
                f.write("Processing your parameters...")
            with open(os.path.join(settings.MEDIA_ROOT, 'files', 'simulation_status_during_run.txt'), "w+") as f:
                f.write("")

            return render(request, 'create_simulation_result.html',
                          {'video_path': settings.MEDIA_URL + 'videos/',  # todo: change "video_path"
                           'previous_page': "upload_files",
                           'num_of_proteins': 0})  # num_of_proteins is irrelevant
    else:
        form = UploadFiles()
    return render(request, 'file_upload.html', {'form': form})
