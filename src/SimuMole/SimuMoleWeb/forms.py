from django import forms

import requests

from django.core.files.uploadedfile import UploadedFile
from django.core.files.storage import FileSystemStorage
from django.conf import settings
from SimuMoleScripts.basicTrajectoryBuilder import scr
from SimuMoleScripts.fix_pdb import fix_pdb
from SimuMoleScripts.transformations import get_atoms, get_atoms_string, translate_vecs
import os
import urllib.request


class SimulationForm0_LoadPdb(forms.Form):
    num_of_proteins = forms.ChoiceField(
        required=True,
        label='Choose whether to load one or two proteins',
        choices=[('1', 'one protein'), ('2', 'two proteins')],
        widget=forms.RadioSelect)

    # first pdb
    first_pdb_type = forms.ChoiceField(
        required=False,
        label='Choose how to load the first pdb',
        choices=[('by_id', 'by id'), ('by_file', 'by file')],
        widget=forms.RadioSelect)
    first_pdb_id = forms.CharField(required=False, label='Enter your first pdb ID', max_length=10)
    first_pdb_file = forms.FileField(required=False, label='Upload your first pdb file')

    # second pdb
    second_pdb_type = forms.ChoiceField(
        required=False,
        label='Choose how to load the second pdb',
        choices=[('by_id', 'by id'), ('by_file', 'by file')],
        widget=forms.RadioSelect)
    second_pdb_id = forms.CharField(required=False, label='Enter your second pdb ID', max_length=10)
    second_pdb_file = forms.FileField(required=False, label='Upload your second pdb file')

    @staticmethod
    def save_file(file: UploadedFile, filename: str):
        file_storage = FileSystemStorage(location=os.path.join(settings.MEDIA_ROOT, 'files'))
        file_storage.delete(filename)  # delete existing file with same name (due to clean_my_file previous calls)
        file_storage.save(filename, file)  # save existing file

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
            first_pdb_file: UploadedFile = self.cleaned_data['first_pdb_file']
            if (num_of_proteins == '1' or num_of_proteins == '2') and first_pdb_type == 'by_file':
                if first_pdb_file == '':
                    raise forms.ValidationError("This field is required.")
                else:
                    self.pdb_file_validation(first_pdb_file)
                self.save_file(first_pdb_file, "_1_.pdb")
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
            second_pdb_file: UploadedFile = self.cleaned_data['second_pdb_file']
            if (num_of_proteins == '2') and second_pdb_type == 'by_file':
                if second_pdb_file == '':
                    raise forms.ValidationError("This field is required.")
                else:
                    self.pdb_file_validation(second_pdb_file)
                self.save_file(second_pdb_file, "_2_.pdb")
            return second_pdb_file

    # pdb validation checks:

    def pdb_id_validation(self, pdb_id):
        if not self.pdb_id_exists(pdb_id):
            raise forms.ValidationError("invalid PDB id")

        if not self.pdb_id_valid(pdb_id):
            raise forms.ValidationError("Protein not supported by OpenMM")

    def pdb_file_validation(self, pdb_id):
        if not self.pdb_file_valid(pdb_id):
            raise forms.ValidationError("Protein not supported by OpenMM")

    @staticmethod
    def pdb_id_exists(pdb_id):
        r = requests.get("https://files.rcsb.org/download/" + pdb_id + ".pdb").status_code
        if r == 404:
            return False
        return True

    '''
    Checks whether a pdb id can be used in an openMM simulation by downloading the relevant file and testing it.
    This function assumes the id is valid.
    '''
    @staticmethod
    def pdb_id_valid(pdb_id):
        pdb_file = "media/files/"+pdb_id+".pdb"
        urllib.request.urlretrieve("https://files.rcsb.org/view/"+pdb_id+".pdb", filename=pdb_file)
        return SimulationForm0_LoadPdb.pdb_file_valid(pdb_file)

    '''
    Checks if the given file can be used in an openMM simulation.
    If it can run without fixing then it returns true.
    If it can run with fixing it will return true AND fix the file.
    Otherwise it returns false.
    '''
    @staticmethod
    def pdb_file_valid(pdb_file):
        dcd_file = "media/files/trajectory.dcd"
        min_steps = 2000    # Steps for 1 frame
        default_temp = 293  # Room temperature in Kelvin
        fix_not_needed = True
        try:
            scr(pdb_file, min_steps, default_temp)
        except Exception:
            fix_not_needed = False
        finally:
            if os.path.exists(dcd_file):
                os.remove(dcd_file)

        if fix_not_needed:
            return True

        try:
            fix_pdb(pdb_file)
            scr(pdb_file, min_steps, default_temp)
        except Exception:
            return False
        finally:
            if os.path.exists(dcd_file):
                os.remove(dcd_file)
        return True


class SimulationForm1_DetermineRelativePosition(forms.Form):
    x1 = forms.FloatField(required=True, label='Enter delta x of first pdb')
    y1 = forms.FloatField(required=True, label='Enter delta y of first pdb')
    z1 = forms.FloatField(required=True, label='Enter delta z of first pdb')
    x2 = forms.FloatField(required=True, label='Enter delta x of second pdb')
    y2 = forms.FloatField(required=True, label='Enter delta y of second pdb')
    z2 = forms.FloatField(required=True, label='Enter delta z of second pdb')
    degXY_1 = forms.FloatField(required=True, label='Enter degrees to rotate from left to right of first pdb')
    degYZ_1 = forms.FloatField(required=True, label='Enter degrees to rotate from down to up of first pdb')
    degXY_2 = forms.FloatField(required=True, label='Enter degrees to rotate from left to right of second pdb')
    degYZ_2 = forms.FloatField(required=True, label='Enter degrees to rotate from down to up of second pdb')

    def clean(self):
        cleaned_data = super(SimulationForm1_DetermineRelativePosition, self).clean()
        data = {**self.initial, **cleaned_data}  # self.initial->from previous steps, cleaned_data->from current step

        if not self.position_is_valid(data['x1'], data['y1'], data['z1'], data['x2'], data['y2'], data['z2']
                                      , data['first_pdb_id'], data['second_pdb_id'], data['first_pdb_type'],
                                      data['second_pdb_type'], data['first_pdb_file'], data['second_pdb_file']):
            raise forms.ValidationError("Positions are not possible: The proteins collide with each other")

        return cleaned_data

    @staticmethod
    def position_is_valid(x1, y1, z1, x2, y2, z2, first_pdb_id, second_pdb_id, first_pdb_type, second_pdb_type,
                          first_pdb_file, second_pdb_file):
        # DONE 6: check with PyMol that the proteins do not collide with each other (need to add the pdbs parameters)

        # return max X,Y,Z locations from all the atoms in vecs
        def get_max_XYZ(vecs):
            return max(vecs, key=lambda v: v[0])[0], max(vecs, key=lambda v: v[1])[1], max(vecs, key=lambda v: v[2])[2]

        # return min X,Y,Z locations from all the atoms in vecs
        def get_min_XYZ(vecs):
            return min(vecs, key=lambda v: v[0])[0], min(vecs, key=lambda v: v[1])[1], min(vecs, key=lambda v: v[2])[2]

        # get the atoms of the first protein after moving it in x1,y1,z1
        if first_pdb_type == 'by_id':
            vecs1 = get_atoms_string(requests.get('https://files.rcsb.org/view/' + first_pdb_id + '.pdb').text)
        else:
            vecs1 = get_atoms('media/files/_1_.pdb')
        translate_vecs(x1, y1, z1, vecs1)

        # get the atoms of the second protein after moving it in x2,y2,z2
        if second_pdb_type == 'by_id':
            vecs2 = get_atoms_string(requests.get("https://files.rcsb.org/view/" + second_pdb_id + ".pdb").text)
        else:
            vecs2 = get_atoms('media/files/_2_.pdb')
        translate_vecs(x2, y2, z2, vecs2)

        maxX1, maxY1, maxZ1 = get_max_XYZ(vecs1)
        maxX2, maxY2, maxZ2 = get_max_XYZ(vecs2)

        minX1, minY1, minZ1 = get_min_XYZ(vecs1)
        minX2, minY2, minZ2 = get_min_XYZ(vecs2)

        dist = 1

        # check overlap in axis X, axis Y and axis Z
        resultX = (maxX1 + dist) >= minX2 and (maxX2 + dist) >= minX1
        resultY = (maxY1 + dist) >= minY2 and (maxY2 + dist) >= minY1
        resultZ = (maxZ1 + dist) >= minZ2 and (maxZ2 + dist) >= minZ1

        # check overlap of whole "boxes" of proteins
        isOverlap = resultX and resultY and resultZ

        return not isOverlap



class SimulationForm2_SimulationParameters(forms.Form):
    temperature = forms.FloatField(required=True, label='Enter temperature (Kelvin)')
    production_steps = forms.IntegerField(required=True, label='Enter number of production steps (1000 = 1 frame)')

    # todo 7: add field of number of time steps (and also: size of every time step)

    def clean(self):
        cleaned_data = super(SimulationForm2_SimulationParameters, self).clean()
        data = {**self.initial, **cleaned_data}  # self.initial->from previous steps, cleaned_data->from current step
        return cleaned_data
