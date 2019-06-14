from django import forms

import requests

from django.core.files.uploadedfile import UploadedFile
from django.core.files.storage import FileSystemStorage
from django.conf import settings
from SimuMoleScripts.basicTrajectoryBuilder import scr_for_checks
from SimuMoleScripts.fix_pdb import fix_pdb
from SimuMoleScripts.transformations import get_atoms, get_atoms_string, translate_vecs, rotate_molecular
from SimuMoleScripts.uploaded_simulation import pdb_and_dcd_match
from SimuMoleScripts.clean_status_csv import init_clean_status, read_clean_status, write_clean_status, \
    same_form__SimulationForm0_LoadPdb, complete_cleaning__SimulationForm0_LoadPdb

import os
from os import path

################################
#   Create Simulation
################################

do_checks_cnt = 0
dir_path = 'media/files/'
clean_status_path = dir_path + 'clean_status.csv'


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
        print("clean")
        cleaned_data = super(SimulationForm0_LoadPdb, self).clean()

        # we already run "clean" and the form isn't changed:
        if path.exists(clean_status_path) and same_form__SimulationForm0_LoadPdb(dir_path, cleaned_data):
            if read_clean_status(dir_path, "SimulationForm0_LoadPdb") == "pass":
                return cleaned_data
            if complete_cleaning__SimulationForm0_LoadPdb(dir_path, cleaned_data):
                write_clean_status(dir_path, "SimulationForm0_LoadPdb", "pass")
                return cleaned_data
        # this is the first time we run "clean":
        else:
            init_clean_status(dir_path, cleaned_data)

        return cleaned_data

    # ................. first pdb validation ................. #

    def clean_first_pdb_type(self):
        num_of_proteins = self.cleaned_data['num_of_proteins']
        first_pdb_type = self.cleaned_data['first_pdb_type']
        if num_of_proteins == '1' or num_of_proteins == '2':
            if first_pdb_type == '':
                raise forms.ValidationError("This field is required.")

        if not path.exists(clean_status_path) or read_clean_status(dir_path, "clean_first_pdb_type") == "pass":
            return self.cleaned_data['first_pdb_type']

        write_clean_status(dir_path, "first_pdb_type", str(first_pdb_type))
        write_clean_status(dir_path, "clean_first_pdb_type", "pass")
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
                    if not path.exists(clean_status_path) or \
                            read_clean_status(dir_path, "clean_first_pdb_id") == "pass" or \
                            (read_clean_status(dir_path, "first_pdb_id") != str(self.cleaned_data['first_pdb_id'])):
                        return first_pdb_id
                    else:
                        self.pdb_id_validation(first_pdb_id, "_1_.pdb")
                        write_clean_status(dir_path, "first_pdb_id", str(first_pdb_id))
                        write_clean_status(dir_path, "clean_first_pdb_id", "pass")
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
                    if not path.exists(clean_status_path) or \
                            read_clean_status(dir_path, "clean_first_pdb_file") == "pass" or \
                            (read_clean_status(dir_path, "first_pdb_file") != str(self.cleaned_data['first_pdb_file'])):
                        return first_pdb_file
                    else:
                        self.pdb_file_validation("media/files/" + first_pdb_file.name)
                        write_clean_status(dir_path, "first_pdb_file", str(first_pdb_file))
                        write_clean_status(dir_path, "clean_first_pdb_file", "pass")
                self.save_file(first_pdb_file, "_1_.pdb")
            return first_pdb_file

    # ................. second pdb validation ................. #

    def clean_second_pdb_type(self):
        num_of_proteins = self.cleaned_data['num_of_proteins']
        second_pdb_type = self.cleaned_data['second_pdb_type']
        if num_of_proteins == '2':
            if second_pdb_type == '':
                raise forms.ValidationError("This field is required.")

        if not path.exists(clean_status_path) or read_clean_status(dir_path, "clean_second_pdb_type") == "pass":
            return second_pdb_type if num_of_proteins == '2' else ''

        write_clean_status(dir_path, "second_pdb_type", str(second_pdb_type))
        write_clean_status(dir_path, "clean_second_pdb_type", "pass")
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
                    if not path.exists(clean_status_path) or \
                            read_clean_status(dir_path, "clean_second_pdb_id") == "pass" or \
                            (read_clean_status(dir_path, "second_pdb_id") != str(self.cleaned_data['second_pdb_id'])):
                        return second_pdb_id
                    else:
                        self.pdb_id_validation(second_pdb_id, "_2_.pdb")
                        write_clean_status(dir_path, "second_pdb_id", str(second_pdb_id))
                        write_clean_status(dir_path, "clean_second_pdb_id", "pass")
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
                    if not path.exists(clean_status_path) or \
                            read_clean_status(dir_path, "clean_second_pdb_file") == "pass" or \
                            (read_clean_status(dir_path, "first_second_file") !=
                             str(self.cleaned_data['second_pdb_file'])):
                        return second_pdb_type
                    else:
                        self.pdb_file_validation("media/files/" + second_pdb_file.name)
                        write_clean_status(dir_path, "second_pdb_file", str(second_pdb_file))
                        write_clean_status(dir_path, "clean_second_pdb_file", "pass")
                self.save_file(second_pdb_file, "_2_.pdb")
            return second_pdb_file

    # ................. pdb validation checks ................. #

    def pdb_id_validation(self, pdb_id, filename):
        if not self.pdb_id_exists(pdb_id):
            raise forms.ValidationError("invalid PDB id")

        if not self.pdb_id_valid(pdb_id, filename):
            raise forms.ValidationError("Protein not supported by OpenMM")

    def pdb_file_validation(self, pdb_id):
        if not self.pdb_file_valid(pdb_id):
            raise forms.ValidationError("Protein not supported by OpenMM")

    @staticmethod
    def pdb_id_exists(pdb_id):
        print("......... pdb_id_exists")
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
    def pdb_file_valid(pdb_file):
        """
        Checks if the given file can be used in an openMM simulation.
        If it can run without fixing then it returns true.
        If it can run with fixing it will return true AND fix the file.
        Otherwise it returns false.
        "pdb_file" == file name
        """
        print("......... pdb_file_valid by openmm")
        dcd_file = "media/files/scr_for_checks.dcd"

        fix_not_needed = True
        try:
            scr_for_checks(pdb_file)
        except Exception as e:
            print("                    ! exp 1: {}".format(str(e)))
            fix_not_needed = False
        finally:
            if os.path.exists(dcd_file):
                os.remove(dcd_file)

        if fix_not_needed:
            return True

        try:
            fix_pdb(pdb_file)
            scr_for_checks(pdb_file)
        except Exception as e:
            print("                    ! exp 2: {}".format(str(e)))
            return False
        finally:
            if os.path.exists(dcd_file):
                os.remove(dcd_file)
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

        if not self.position_is_valid(data['x1'], data['y1'], data['z1'],
                                      data['x2'], data['y2'], data['z2'],
                                      data['degXY_1'], data['degYZ_1'], data['degXY_2'], data['degYZ_2'],
                                      data['first_pdb_id'], data['second_pdb_id'], data['first_pdb_type'],
                                      data['second_pdb_type'], data['first_pdb_file'], data['second_pdb_file']):
            raise forms.ValidationError("Positions are not possible: The proteins collide with each other")

        return cleaned_data

    @staticmethod
    def position_is_valid(x1, y1, z1, x2, y2, z2, degXY_1, degYZ_1, degXY_2, degYZ_2,
                          first_pdb_id, second_pdb_id, first_pdb_type, second_pdb_type,
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
        rotate_molecular(x1, y1, z1, degXY_1, degYZ_1, vecs1)

        # get the atoms of the second protein after moving it in x2,y2,z2
        if second_pdb_type == 'by_id':
            vecs2 = get_atoms_string(requests.get("https://files.rcsb.org/view/" + second_pdb_id + ".pdb").text)
        else:
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

    # time_step_duration = forms.FloatField(required=False,
    #                                       label='Enter the duration of the time step (in fs/femto-second)',
    #                                       initial=2.0) // todo: check if we can add this field
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
