from django import forms


class CreateSimulationForm(forms.Form):
    first_pdb = forms.CharField(label='Enter your first pdb ID', max_length=10)
    second_pdb = forms.CharField(label='Enter your second pdb ID', max_length=10)
    x1 = forms.CharField(label='Enter delta x of first pdb', max_length=10)
    y1 = forms.CharField(label='Enter delta y of first pdb', max_length=10)
    z1 = forms.CharField(label='Enter delta z of first pdb', max_length=10)
    temperature = forms.CharField(label='Enter temperature', max_length=10)
