from django import forms


# todo 1: replace CharField with DecimalField, when needed
# todo 2: add validation tests (examples below)

class SimulationForm1_LoadPdb(forms.Form):
    first_pdb = forms.CharField(
        required=True,
        label='Enter your first pdb ID',
        max_length=10)
    second_pdb = forms.CharField(
        required=True,
        label='Enter your second pdb ID',
        max_length=10)

    def clean(self):
        cleaned_data = super(SimulationForm1_LoadPdb, self).clean()
        data = {**self.initial, **cleaned_data}  # self.initial->from previous steps, cleaned_data->from current step


class SimulationForm2_DetermineRelativePosition(forms.Form):
    x1 = forms.CharField(
        required=True,
        label='Enter delta x of first pdb',
        max_length=10)
    y1 = forms.CharField(
        required=True,
        label='Enter delta y of first pdb',
        max_length=10)
    z1 = forms.CharField(
        required=True,
        label='Enter delta z of first pdb',
        max_length=10)
    x2 = forms.CharField(
        required=True,
        label='Enter delta x of second pdb',
        max_length=10)
    y2 = forms.CharField(
        required=True,
        label='Enter delta y of second pdb',
        max_length=10)
    z2 = forms.CharField(
        required=True,
        label='Enter delta z of second pdb',
        max_length=10)

    def clean(self):
        cleaned_data = super(SimulationForm2_DetermineRelativePosition, self).clean()
        data = {**self.initial, **cleaned_data}  # self.initial->from previous steps, cleaned_data->from current step


class SimulationForm3_SimulationParameters(forms.Form):
    temperature = forms.CharField(
        required=True,
        label='Enter temperature',
        max_length=10)

    def clean(self):
        cleaned_data = super(SimulationForm3_SimulationParameters, self).clean()
        data = {**self.initial, **cleaned_data}  # self.initial->from previous steps, cleaned_data->from current step
