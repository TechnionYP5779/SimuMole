from django import forms

import requests

from django.core.files.uploadedfile import UploadedFile
from django.core.files.storage import FileSystemStorage
from django.conf import settings
from SimuMoleScripts.basicTrajectoryBuilder import scr_for_checks
from SimuMoleScripts.fix_pdb import fix_pdb
from SimuMoleScripts.transformations import get_atoms, get_atoms_string, translate_vecs, rotate_molecular, translate_pdb
from SimuMoleScripts.uploaded_simulation import pdb_and_dcd_match

import os
import pymol


################################
#   Create Simulation
################################


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

        errors = []

        num_of_proteins = data.get('num_of_proteins', None)
        first_pdb_type = data.get('first_pdb_type', None)
        first_pdb_id, first_pdb_file = data.get('first_pdb_id', None), data.get('first_pdb_file', None)
        second_pdb_type = data.get('second_pdb_type', None)
        second_pdb_id, second_pdb_file = data.get('second_pdb_id', None), data.get('second_pdb_file', None)

        # clean_first_pdb_type
        if num_of_proteins == '1' or num_of_proteins == '2':
            if (first_pdb_type is None) or (first_pdb_type == ''):
                errors.append(forms.ValidationError("First protein: The 'PDB type' field is required."))

        # clean_first_pdb_id
        if first_pdb_type is not None:
            if (num_of_proteins == '1' or num_of_proteins == '2') and first_pdb_type == 'by_id':
                if (first_pdb_id is None) or (first_pdb_id == ''):
                    errors.append(forms.ValidationError("First protein: The 'by id' field is required."))
                else:
                    pdb_validation_result = self.pdb_id_validation(first_pdb_id, "_1_.pdb")
                    if pdb_validation_result is not None:
                        errors.append(forms.ValidationError("First protein: " + pdb_validation_result))

        # clean_first_pdb_file
        if first_pdb_type is not None:
            first_pdb_file: UploadedFile = self.cleaned_data['first_pdb_file']
            if (num_of_proteins == '1' or num_of_proteins == '2') and first_pdb_type == 'by_file':
                if (first_pdb_file is None) or (first_pdb_file == ''):
                    errors.append(forms.ValidationError("First protein: The 'by file' field is required."))
                else:
                    self.save_file(first_pdb_file, "_1_.pdb")
                    pdb_validation_result = self.pdb_file_validation("media/files/" + "_1_.pdb")
                    if pdb_validation_result is not None:
                        errors.append(forms.ValidationError("First protein: " + pdb_validation_result))

        # clean_second_pdb_type
        if num_of_proteins == '2':
            if (second_pdb_type is None) or (second_pdb_type == ''):
                errors.append(forms.ValidationError("Second protein: The 'PDB type' field is required."))

        # clean_second_pdb_id
        if second_pdb_type is not None:
            if (num_of_proteins == '2') and second_pdb_type == 'by_id':
                if (second_pdb_id is None) or (second_pdb_id == ''):
                    errors.append(forms.ValidationError("Second protein: The 'by id' field is required."))
                else:
                    pdb_validation_result = self.pdb_id_validation(second_pdb_id, "_2_.pdb")
                    if pdb_validation_result is not None:
                        errors.append(forms.ValidationError("Second protein: " + pdb_validation_result))

        # clean_second_pdb_file
        if second_pdb_type is not None:
            second_pdb_file: UploadedFile = self.cleaned_data['second_pdb_file']
            if (num_of_proteins == '2') and second_pdb_type == 'by_file':
                if (second_pdb_file is None) or (second_pdb_file == ''):
                    errors.append(forms.ValidationError("Second protein: The 'by file' field is required."))
                else:
                    self.save_file(second_pdb_file, "_2_.pdb")
                    pdb_validation_result = self.pdb_file_validation("media/files/" + "_2_.pdb")
                    if pdb_validation_result is not None:
                        errors.append(forms.ValidationError("Second protein: " + pdb_validation_result))

        if len(errors) != 0:
            raise forms.ValidationError(errors)

        return cleaned_data

    def pdb_id_validation(self, pdb_id, filename):
        if not self.pdb_id_exists(pdb_id):
            return "Invalid PDB id"
        if not self.pdb_id_valid(pdb_id, filename):
            return "Protein not supported by OpenMM"
        return None

    def pdb_file_validation(self, pdb_id):
        if not self.pdb_file_valid(pdb_id):
            return "Protein not supported by OpenMM"
        return None

    @staticmethod
    def pdb_id_exists(pdb_id):
        r = requests.get("https://files.rcsb.org/download/" + pdb_id + ".pdb").status_code
        if r == 404:
            return False
        return True

    @staticmethod
    def download_pdb(pdb_id, filename):
        response = requests.get("https://files.rcsb.org/view/" + pdb_id + ".pdb")
        pdb_file_name = "media/files/" + filename
        pdb_file = open(pdb_file_name, 'w')
        pdb_file.write(response.text)
        pdb_file.close()
        response.close()
        return pdb_file_name

    @staticmethod
    def pdb_id_valid(pdb_id, filename):
        """
        Checks whether a pdb id can be used in an openMM simulation by downloading the relevant file and testing it.
        This function assumes the id is valid.
        """
        pdb_file_name = SimulationForm0_LoadPdb.download_pdb(pdb_id, filename)
        return SimulationForm0_LoadPdb.pdb_file_valid(pdb_file_name)

    @staticmethod
    def pdb_file_valid(pdb_file_name):
        """
        Checks if the given file can be used in an openMM simulation.
        If it can run without fixing then it returns true.
        If it can run with fixing it will return true AND fix the file.
        Otherwise it returns false.
        """
        dcd_file = "media/files/scr_for_checks.dcd"

        fix_not_needed = True
        try:
            scr_for_checks(pdb_file_name)
        except Exception as e:
            fix_not_needed = False
        finally:
            if os.path.exists(dcd_file):
                os.remove(dcd_file)

        if fix_not_needed:
            return True

        try:
            fix_pdb(pdb_file_name)
            scr_for_checks(pdb_file_name)
        except Exception as e:
            return False

        return True


class SimulationForm1_DetermineRelativePosition(forms.Form):
    x1 = forms.FloatField(required=False, label='Enter delta x of first pdb')
    y1 = forms.FloatField(required=False, label='Enter delta y of first pdb')
    z1 = forms.FloatField(required=False, label='Enter delta z of first pdb')
    x2 = forms.FloatField(required=False, label='Enter delta x of second pdb')
    y2 = forms.FloatField(required=False, label='Enter delta y of second pdb')
    z2 = forms.FloatField(required=False, label='Enter delta z of second pdb')
    degXY_1 = forms.FloatField(required=False, label='Enter degrees to rotate from left to right of first pdb')
    degYZ_1 = forms.FloatField(required=False, label='Enter degrees to rotate from down to up of first pdb')
    degXY_2 = forms.FloatField(required=False, label='Enter degrees to rotate from left to right of second pdb')
    degYZ_2 = forms.FloatField(required=False, label='Enter degrees to rotate from down to up of second pdb')

    def clean(self):
        cleaned_data = super(SimulationForm1_DetermineRelativePosition, self).clean()
        data = {**self.initial, **cleaned_data}  # self.initial->from previous steps, cleaned_data->from current step

        for field in ['x1', 'y1', 'z1', 'x2', 'y2', 'z2', 'degXY_1', 'degYZ_1', 'degXY_2', 'degYZ_2']:
            if data[field] == '' or data[field] is None:
                raise forms.ValidationError("All fields are required.")

        if not self.position_is_valid(data['x1'], data['y1'], data['z1'], data['x2'], data['y2'], data['z2'],
                                      data['degXY_1'], data['degYZ_1'], data['degXY_2'], data['degYZ_2']):
            raise forms.ValidationError("Positions are not possible: The proteins collide with each other")

        self.change_relative_position(data['x1'], data['y1'], data['z1'], data['x2'], data['y2'], data['z2'],
                                      data['degXY_1'], data['degYZ_1'], data['degXY_2'], data['degYZ_2'])

        return cleaned_data

    @staticmethod
    def change_relative_position(x1, y1, z1, x2, y2, z2, degXY_1, degYZ_1, degXY_2, degYZ_2):
        """
        Save PDB file that represents the 2 PDB files after you change the positions and running pdb_fixer
        """
        # change positions
        filename_1, filename_2, pdb, temp = '_1_', '_2_', '.pdb', 'media/files/'
        filename_1_movement, filename_2_movement = filename_1 + '__movement', filename_2 + '__movement'
        translate_pdb(temp + filename_1 + pdb, temp + filename_1_movement + pdb, x1, y1, z1, degXY_1, degYZ_1)
        translate_pdb(temp + filename_2 + pdb, temp + filename_2_movement + pdb, x2, y2, z2, degXY_2, degYZ_2)

        # fix pdb
        fix_pdb(temp + filename_1_movement + pdb)
        fix_pdb(temp + filename_2_movement + pdb)

        # merge to single pdb file
        pymol.finish_launching(['pymol', '-q'])  # pymol: -q quiet launch, -c no gui, -e fullscreen
        cmd = pymol.cmd
        cmd.reinitialize()
        cmd.load(temp + filename_1_movement + pdb)
        cmd.load(temp + filename_2_movement + pdb)
        cmd.zoom()
        cmd.save(temp + "both_1_2" + pdb)

    @staticmethod
    def position_is_valid(x1, y1, z1, x2, y2, z2, degXY_1, degYZ_1, degXY_2, degYZ_2):
        """
        check with PyMol that the proteins do not collide with each other
        """

        # return max X,Y,Z locations from all the atoms in vecs
        def get_max_XYZ(vecs):
            return max(vecs, key=lambda v: v[0])[0], max(vecs, key=lambda v: v[1])[1], max(vecs, key=lambda v: v[2])[2]

        # return min X,Y,Z locations from all the atoms in vecs
        def get_min_XYZ(vecs):
            return min(vecs, key=lambda v: v[0])[0], min(vecs, key=lambda v: v[1])[1], min(vecs, key=lambda v: v[2])[2]

        # get the atoms of the first protein after moving it in x1,y1,z1
        vecs1 = get_atoms('media/files/_1_.pdb')
        translate_vecs(x1, y1, z1, vecs1)
        rotate_molecular(x1, y1, z1, degXY_1, degYZ_1, vecs1)

        # get the atoms of the second protein after moving it in x2,y2,z2
        vecs2 = get_atoms('media/files/_2_.pdb')
        translate_vecs(x2, y2, z2, vecs2)
        rotate_molecular(x2, y2, z2, degXY_2, degYZ_2, vecs2)

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
    temperature_scale = forms.ChoiceField(
        required=False,
        label='Choose a temperature scale',
        choices=[('kelvin', 'Kelvin (K)'), ('celsius', 'Celsius (°C)')],
        widget=forms.RadioSelect, initial='celsius')
    temperature = forms.FloatField(required=False, label='Enter temperature', initial=30)

    time_step_number = forms.IntegerField(required=False, label='Enter the number of time steps (frames)', initial=5)

    def clean(self):
        cleaned_data = super(SimulationForm2_SimulationParameters, self).clean()
        data = {**self.initial, **cleaned_data}  # self.initial->from previous steps, cleaned_data->from current step

        errors = []
        for field in ['temperature_scale', 'temperature', 'time_step_number']:
            if data[field] == '' or data[field] is None:
                errors.append(forms.ValidationError("All fields are required."))
            else:
                if field == 'time_step_number':  # the user type some value, therefore we can apply "int(str)"
                    if int(data['time_step_number'] < 2):
                        errors.append(forms.ValidationError("The minimum number of time steps is 2."))

        if len(errors) != 0:
            raise forms.ValidationError(errors)
        return cleaned_data


################################
#   Upload PDB & DCD
################################

class UploadFiles(forms.Form):
    pdb_file = forms.FileField(required=True, label='Upload a PDB file')
    dcd_file = forms.FileField(required=True, label='Upload a DCD file')

    def clean(self):
        errors = []

        # pdb_file
        pdb_file: UploadedFile = self.cleaned_data['pdb_file']
        file_name, file_extension = os.path.splitext(str(pdb_file))
        if file_extension.lower() != '.pdb':
            errors.append(forms.ValidationError("The PDB file format is incorrect"))

        # dcd_file
        dcd_file: UploadedFile = self.cleaned_data['dcd_file']
        file_name, file_extension = os.path.splitext(str(dcd_file))
        if file_extension.lower() != '.dcd':
            errors.append(forms.ValidationError("The DCD file format is incorrect"))

        if len(errors) != 0:
            raise forms.ValidationError(errors)

        # save files (only after checking that their format is correct)
        self.save_file(pdb_file, "file_upload_pdb.pdb")
        self.save_file(dcd_file, "file_upload_dcd.dcd")

        # check match between the files
        if not pdb_and_dcd_match("file_upload_pdb.pdb", "file_upload_dcd.dcd"):
            raise forms.ValidationError("The PDB file does not match the DCD file (the number of atoms is different)")

    @staticmethod
    def save_file(file: UploadedFile, filename: str):
        path = os.path.join(settings.MEDIA_ROOT, 'files')
        file_storage = FileSystemStorage(location=path)

        file_storage.delete(filename)  # delete existing file with same name
        file_storage.save(filename, file)  # save existing file
