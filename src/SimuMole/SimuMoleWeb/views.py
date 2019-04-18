from django.shortcuts import render
from django.http import HttpResponse
from django.http import Http404
from django.template import loader

from formtools.wizard.views import CookieWizardView  # todo: check if need to replace with "SessionWizardView"

from SimuMoleScripts.simulation_main_script import Simulation

from SimuMoleWeb.forms import SimulationForm1_LoadPdb, \
    SimulationForm2_DetermineRelativePosition, \
    SimulationForm3_SimulationParameters


################################
#   Home
################################

def home(request):
    some_dict = {}
    return render(request, 'home.html', some_dict)


################################
#   Create Simulation
################################

class SimulationWizard(CookieWizardView):
    template_name = 'create_simulation.html'

    # override "done": this function is called when the form is submitted
    def done(self, form_list, **kwargs):
        form_data = [form.cleaned_data for form in form_list]

        form1_load_pdb = form_data[0]
        first_pdb = form1_load_pdb['first_pdb']
        second_pdb = form1_load_pdb['second_pdb']

        form2_determine_relative_position = form_data[1]
        x1 = form2_determine_relative_position['x1']
        y1 = form2_determine_relative_position['y1']
        z1 = form2_determine_relative_position['z1']
        x2 = form2_determine_relative_position['x2']
        y2 = form2_determine_relative_position['y2']
        z2 = form2_determine_relative_position['z2']

        form_3_simulation_parameters = form_data[2]
        temperature = form_3_simulation_parameters['temperature']

        s = Simulation(first_pdb, second_pdb, x1, y1, z1, temperature)
        s.create_simulation()

        return render(self.request, 'create_simulation_result.html', {
            'form_data': form_data,
        })

    # override "get_form_initial" for allow getting data from previous steps
    def get_form_initial(self, step):
        current_step = self.storage.current_step
        initial_data = self.initial_dict.get(step, {})  # initial data of current step
        prev_data = {}  # data of previous step

        if current_step == '0':  # SimulationForm1_LoadPdb
            prev_data = {}

        if current_step == '1':  # SimulationForm2_DetermineRelativePosition
            prev_data = self.storage.get_step_data('0')
            prev_data = {'first_pdb': prev_data.get('0-first_pdb', ''), 'second_pdb': prev_data.get('0-second_pdb', '')}

        if current_step == '2':  # SimulationForm3_SimulationParameters
            prev_prev_data = self.storage.get_step_data('0')
            prev_data = self.storage.get_step_data('1')
            prev_data = {
                'first_pdb': prev_prev_data.get('0-first_pdb', ''),
                'second_pdb': prev_prev_data.get('0-second_pdb', ''),
                'x1': prev_data.get('1-x1', ''), 'x2': prev_data.get('1-x2', ''),
                'y1': prev_data.get('1-y1', ''), 'y2': prev_data.get('1-y2', ''),
                'z1': prev_data.get('1-z1', ''), 'z2': prev_data.get('1-z2', ''), }

        update_data = {**prev_data, **initial_data}
        return self.initial_dict.get(step, update_data)
