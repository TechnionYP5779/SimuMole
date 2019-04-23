from django.shortcuts import render
from formtools.wizard.views import CookieWizardView
from django.http import HttpResponse
from django.http import Http404
from django.template import loader
from django.shortcuts import render

from SimuMoleScripts.simulation_main_script import Simulation
from .models import UploadForm, Upload
from django.http import HttpResponseRedirect
from django.urls import reverse
from django.contrib import messages
import os

from .forms import CreateSimulationForm


################################
#   Home
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
#   Create Simulation
################################

class SimulationWizard(CookieWizardView):
    template_name = 'create_simulation.html'

    @staticmethod
    def clean_form_dict(dict_):
        """
        get dictionary that represent the form data, and return clean dictionary data (without contradictions)
        """
        clean_dict = {}
        first_pdb_type, first_pdb_id, first_pdb_file, second_pdb_type, second_pdb_id, second_pdb_file = '', '', '', '', '', ''
        x1, y1, z1, x2, y2, z2 = '', '', '', '', '', ''

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
            x1, y1, z1 = dict_.get('x1', 0), dict_.get('y1', 0), dict_.get('z1', 0)
            x2, y2, z2 = dict_.get('x2', 0), dict_.get('y2', 0), dict_.get('z2', 0)

        temperature = dict_.get('temperature', '')

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
        clean_dict['temperature'] = temperature

        return clean_dict

    def done(self, form_list, **kwargs):
        """
        override "done": this function is called when the form is submitted
        """
        form_data = [form.cleaned_data for form in form_list]
        form_dict = {k: v for d in form_data for k, v in d.items()}  # convert list of dictionaries to one dictionary
        form_dict = self.clean_form_dict(form_dict)

        # todo 8: change parameter list
        # s = Simulation(num_of_proteins, first_pdb, second_pdb, x1, y1, z1, x2, y2, z2, temperature)
        # s.create_simulation()

        return render(self.request, 'create_simulation_result.html', {
            'form_data': form_dict,
        })

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
            else {'temperature': step_2_prev_data.get('2-temperature')}

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
    if request.method=="POST":
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
        file=UploadForm()
    return render(request,'file_upload.html',{'form':file})
