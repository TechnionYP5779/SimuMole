from django import forms


# todo 1: change types of dields

class SimulationForm0_NumberOfProteins(forms.Form):
    num_of_proteins = forms.ChoiceField(
        required=True,
        label='Choose whether to load one or two proteins',
        choices=[('1', 'one protein'), ('2', 'two proteins')],
        widget=forms.RadioSelect)

    def clean(self):
        cleaned_data = super(SimulationForm0_NumberOfProteins, self).clean()
        data = {**self.initial, **cleaned_data}  # self.initial->from previous steps, cleaned_data->from current step
        return cleaned_data


class SimulationForm1_LoadPdb(forms.Form):
    error_css_class = 'error'
    required_css_class = 'required'

    first_pdb = forms.CharField(required=True, label='Enter your first pdb ID', max_length=10)
    second_pdb = forms.CharField(required=True, label='Enter your second pdb ID', max_length=10, disabled=True)

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

    def clean_first_pdb(self):
        first_pdb = self.cleaned_data['first_pdb']

        if not self.pdb_id_exists(first_pdb):
            raise forms.ValidationError("Invalid PDB ID")
        if not self.pdb_id_valid_by_openmm(first_pdb):
            raise forms.ValidationError("Protein not supported by OpenMM")

        return first_pdb

    def clean_second_pdb(self):
        second_pdb = self.cleaned_data['second_pdb']
        if not self.pdb_id_exists(second_pdb):
            raise forms.ValidationError("invalid PDB ID")
        if not self.pdb_id_valid_by_openmm(second_pdb):
            raise forms.ValidationError("Protein not supported by OpenMM")

        return second_pdb

    def clean(self):
        cleaned_data = super(SimulationForm1_LoadPdb, self).clean()
        data = {**self.initial, **cleaned_data}  # self.initial->from previous steps, cleaned_data->from current step
        return cleaned_data

    @staticmethod
    def pdb_id_exists(pdb_id):
        # todo 2: check that the ID is exist
        return True

    @staticmethod
    def pdb_id_valid_by_openmm(pdb_id):
        # todo 3: check that the force field of OpenMM accept this PDB
        return True


class SimulationForm2_DetermineRelativePosition(forms.Form):
    x1 = forms.CharField(required=True, label='Enter delta x of first pdb', max_length=10)
    y1 = forms.CharField(required=True, label='Enter delta y of first pdb', max_length=10)
    z1 = forms.CharField(required=True, label='Enter delta z of first pdb', max_length=10)
    x2 = forms.CharField(required=True, label='Enter delta x of second pdb', max_length=10)
    y2 = forms.CharField(required=True, label='Enter delta y of second pdb', max_length=10)
    z2 = forms.CharField(required=True, label='Enter delta z of second pdb', max_length=10)

    def clean(self):
        cleaned_data = super(SimulationForm2_DetermineRelativePosition, self).clean()
        data = {**self.initial, **cleaned_data}  # self.initial->from previous steps, cleaned_data->from current step

        if not self.position_is_valid(data['first_pdb'], data['second_pdb'],
                                      data['x1'], data['y1'], data['z1'], data['x2'], data['y2'], data['z2']):
            raise forms.ValidationError("Positions are not possible: The proteins collide with each other")

        return cleaned_data

    @staticmethod
    def position_is_valid(pdb_id_1, pdb_id_2, x1, y1, z1, x2, y2, z2):
        # todo 4: check with PyMol that the proteins do not collide with each other
        return True


class SimulationForm3_SimulationParameters(forms.Form):
    temperature = forms.CharField(required=True, label='Enter temperature', max_length=10)

    # todo 5: add field of number of time steps (and also: size of every time step)

    def clean(self):
        cleaned_data = super(SimulationForm3_SimulationParameters, self).clean()
        data = {**self.initial, **cleaned_data}  # self.initial->from previous steps, cleaned_data->from current step
        return cleaned_data
