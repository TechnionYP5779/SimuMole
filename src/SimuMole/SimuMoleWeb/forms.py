from django import forms


class LastForm(forms.Form):
    first_pdb_file = forms.FileField(required=False)
    second_pdb_file = forms.FileField(required=False)


class SimulationForm0_LoadPdb(forms.Form):
    num_of_proteins = forms.ChoiceField(
        required=True,
        label='Choose whether to load one or two proteins',
        choices=[('1', 'one protein'), ('2', 'two proteins')],
        widget=forms.RadioSelect)

    # todo 1: change types of fields "first_pdb_file" and "second_pdb_file"

    # first pdb
    first_pdb_type = forms.ChoiceField(
        required=False,
        label='Choose how to load the first pdb',
        choices=[('by_id', 'by id'), ('by_file', 'by file')],
        widget=forms.RadioSelect)
    first_pdb_id = forms.CharField(required=False, label='Enter your first pdb ID', max_length=10)
    first_pdb_file = forms.CharField(required=False, label='Upload your first pdb file', max_length=10)

    # second pdb
    second_pdb_type = forms.ChoiceField(
        required=False,
        label='Choose how to load the second pdb',
        choices=[('by_id', 'by id'), ('by_file', 'by file')],
        widget=forms.RadioSelect)
    second_pdb_id = forms.CharField(required=False, label='Enter your second pdb ID', max_length=10)
    second_pdb_file = forms.CharField(required=False, label='Upload your second pdb file', max_length=10)

    def clean(self):
        cleaned_data = super(SimulationForm0_LoadPdb, self).clean()
        data = {**self.initial, **cleaned_data}  # self.initial->from previous steps, cleaned_data->from current step
        return cleaned_data

    # first pdb validation:

    def clean_first_pdb_type(self):
        num_of_proteins = self.cleaned_data['num_of_proteins']
        first_pdb_type = self.cleaned_data['first_pdb_type']
        if num_of_proteins == '1' or num_of_proteins == '2':
            if first_pdb_type == '':
                raise forms.ValidationError("This field is required.")
        return first_pdb_type

    def clean_first_pdb_id(self):
        num_of_proteins = self.cleaned_data['num_of_proteins']
        if self.cleaned_data.get('first_pdb_type') is not None:
            first_pdb_type = self.cleaned_data['first_pdb_type']
            first_pdb_id = self.cleaned_data['first_pdb_id']
            if (num_of_proteins == '1' or num_of_proteins == '2') and first_pdb_type == 'by_id':
                if first_pdb_id == '':
                    raise forms.ValidationError("This field is required.")
                else:
                    self.pdb_id_validation(first_pdb_id)
            return first_pdb_id

    def clean_first_pdb_file(self):
        num_of_proteins = self.cleaned_data['num_of_proteins']
        if self.cleaned_data.get('first_pdb_type') is not None:
            first_pdb_type = self.cleaned_data['first_pdb_type']
            first_pdb_file = self.cleaned_data['first_pdb_file']
            if (num_of_proteins == '1' or num_of_proteins == '2') and first_pdb_type == 'by_file':
                if first_pdb_file == '':
                    raise forms.ValidationError("This field is required.")
                else:
                    self.pdb_file_validation(first_pdb_file)
            return first_pdb_file

    # second pdb validation:

    def clean_second_pdb_type(self):
        num_of_proteins = self.cleaned_data['num_of_proteins']
        second_pdb_type = self.cleaned_data['second_pdb_type']
        if num_of_proteins == '2':
            if second_pdb_type == '':
                raise forms.ValidationError("This field is required.")
        return second_pdb_type

    def clean_second_pdb_id(self):
        num_of_proteins = self.cleaned_data['num_of_proteins']
        if self.cleaned_data.get('second_pdb_type') is not None:
            second_pdb_type = self.cleaned_data['second_pdb_type']
            second_pdb_id = self.cleaned_data['second_pdb_id']
            if (num_of_proteins == '2') and second_pdb_type == 'by_id':
                if second_pdb_id == '':
                    raise forms.ValidationError("This field is required.")
                else:
                    self.pdb_id_validation(second_pdb_id)
            return second_pdb_id

    def clean_second_pdb_file(self):
        num_of_proteins = self.cleaned_data['num_of_proteins']
        if self.cleaned_data.get('second_pdb_type') is not None:
            second_pdb_type = self.cleaned_data['second_pdb_type']
            second_pdb_file = self.cleaned_data['second_pdb_file']
            if (num_of_proteins == '2') and second_pdb_type == 'by_file':
                if second_pdb_file == '':
                    raise forms.ValidationError("This field is required.")
                else:
                    self.pdb_file_validation(second_pdb_file)
            return second_pdb_file

    # pdb validation checks:

    def pdb_id_validation(self, pdb_id):
        if not self.pdb_id_exists(pdb_id):
            raise forms.ValidationError("invalid PDB id")
        if not self.pdb_id_valid_by_openmm(pdb_id):
            raise forms.ValidationError("Protein not supported by OpenMM")

    def pdb_file_validation(self, pdb_id):
        if not self.pdb_file_valid(pdb_id):
            raise forms.ValidationError("invalid PDB file")
        if not self.pdb_file_valid_by_openmm(pdb_id):
            raise forms.ValidationError("Protein not supported by OpenMM")

    @staticmethod
    def pdb_id_exists(pdb_id):
        # todo 2: check that the ID is exist
        return True

    @staticmethod
    def pdb_id_valid_by_openmm(pdb_id):
        # todo 3: check that the force field of OpenMM accept this PDB
        return True

    @staticmethod
    def pdb_file_valid(pdb_file):
        # todo 4: check that the file is valid
        return True

    @staticmethod
    def pdb_file_valid_by_openmm(pdb_file):
        # todo 5: check that the force field of OpenMM accept this PDB
        return True


class SimulationForm1_DetermineRelativePosition(forms.Form):
    x1 = forms.FloatField(required=True, label='Enter delta x of first pdb')
    y1 = forms.FloatField(required=True, label='Enter delta y of first pdb')
    z1 = forms.FloatField(required=True, label='Enter delta z of first pdb')
    x2 = forms.FloatField(required=True, label='Enter delta x of second pdb')
    y2 = forms.FloatField(required=True, label='Enter delta y of second pdb')
    z2 = forms.FloatField(required=True, label='Enter delta z of second pdb')

    def clean(self):
        cleaned_data = super(SimulationForm1_DetermineRelativePosition, self).clean()
        data = {**self.initial, **cleaned_data}  # self.initial->from previous steps, cleaned_data->from current step

        if not self.position_is_valid(data['x1'], data['y1'], data['z1'], data['x2'], data['y2'], data['z2']):
            raise forms.ValidationError("Positions are not possible: The proteins collide with each other")

        return cleaned_data

    @staticmethod
    def position_is_valid(x1, y1, z1, x2, y2, z2):
        # todo 6: check with PyMol that the proteins do not collide with each other (need to add the pdbs parameters)
        return True


class SimulationForm2_SimulationParameters(forms.Form):
    temperature = forms.FloatField(required=True, label='Enter temperature')

    # todo 7: add field of number of time steps (and also: size of every time step)

    def clean(self):
        cleaned_data = super(SimulationForm2_SimulationParameters, self).clean()
        data = {**self.initial, **cleaned_data}  # self.initial->from previous steps, cleaned_data->from current step
        return cleaned_data
