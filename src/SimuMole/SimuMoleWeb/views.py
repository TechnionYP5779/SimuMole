from django.shortcuts import render
from formtools.wizard.views import CookieWizardView

from SimuMoleScripts.simulation_main_script import Simulation


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
        """
        override "done": this function is called when the form is submitted
        """
        form_data = [form.cleaned_data for form in form_list]
        form_dict = {k: v for d in form_data for k, v in d.items()}  # convert list of dictionaries to one dictionary

        num_of_proteins = form_dict.get('num_of_proteins', '')
        first_pdb, second_pdb = form_dict.get('first_pdb', ''), form_dict.get('first_pdb', '')
        x1, y1, z1 = form_dict.get('x1', 0), form_dict.get('y1', 0), form_dict.get('z1', 0)
        x2, y2, z2 = form_dict.get('x2', 0), form_dict.get('y2', 0), form_dict.get('z2', 0)
        temperature = form_dict.get('temperature', '')

        s = Simulation(num_of_proteins, first_pdb, second_pdb, x1, y1, z1, x2, y2, z2, temperature)
        s.create_simulation()

        return render(self.request, 'create_simulation_result.html', {
            'form_data': form_data,
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
                  'first_pdb': step_0_prev_data.get('0-first_pdb'),
                  'second_pdb': step_0_prev_data.get('0-second_pdb')}

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
