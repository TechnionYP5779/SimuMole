from django.shortcuts import render, render_to_response
from django.http import HttpResponse
from django.http import Http404
from django.template import loader

from SimuMoleScripts.simulation_main_script import Simulation

from formtools.wizard.views import CookieWizardView  # todo: check if need to replace with "SessionWizardView"


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
