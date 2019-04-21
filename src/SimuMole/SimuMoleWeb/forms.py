from django import forms


class SimulationForm0_NumberOfProteins(forms.Form):
    num_of_proteins = forms.ChoiceField(
        required=True,
        label='Choose whether to load one or two proteins',
        choices=[('1', 'one protein'), ('2', 'two proteins')],
        widget=forms.RadioSelect)

    def clean(self):
        cleaned_data = super(SimulationForm0_NumberOfProteins, self).clean()
        data = {**self.initial, **cleaned_data}  # self.initial->from previous steps, cleaned_data->from current step
        # print("___0___ number of pdb")
        # print("___0______" + str(data))


class SimulationForm1_LoadPdb(forms.Form):
    first_pdb = forms.CharField(
        required=True,
        label='Enter your first pdb ID',
        max_length=10)
    second_pdb = forms.CharField(
        required=True,
        label='Enter your second pdb ID',
        max_length=10, disabled=True)

    def __init__(self, *args, **kwargs):
        super(SimulationForm1_LoadPdb, self).__init__(*args, **kwargs)

        # disabled/enable the field 'second_pdb' depending on the value of 'num_of_proteins':
        num_of_proteins = self.initial['num_of_proteins']
        if num_of_proteins == '1':
            self.fields['second_pdb'] = forms.CharField(required=False, label='Enter your second pdb ID', disabled=True,
                                                        initial="Field is not relevant, you have chosen to upload a single PDB file")
        if num_of_proteins == '2':
            self.fields['second_pdb'] = forms.CharField(required=False, label='Enter your second pdb ID', max_length=10,
                                                        disabled=False)

    def clean(self):
        cleaned_data = super(SimulationForm1_LoadPdb, self).clean()
        data = {**self.initial, **cleaned_data}  # self.initial->from previous steps, cleaned_data->from current step
        # print("___1___ load pdb")
        # print("___1______" + str(data))


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
        # print("___2___ positions")
        # print("___2______" + str(data))


class SimulationForm3_SimulationParameters(forms.Form):
    temperature = forms.CharField(
        required=True,
        label='Enter temperature',
        max_length=10)

    def clean(self):
        cleaned_data = super(SimulationForm3_SimulationParameters, self).clean()
        data = {**self.initial, **cleaned_data}  # self.initial->from previous steps, cleaned_data->from current step
        # print("___3___ simulation parameters")
        # print("___3______" + str(data))
