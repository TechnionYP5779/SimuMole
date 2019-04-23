from django import forms


class CreateSimulationForm(forms.Form):
    first_pdb = forms.CharField(label='Enter your first pdb ID', max_length=10)
    second_pdb = forms.CharField(label='Enter your second pdb ID', max_length=10)
    x1 = forms.FloatField(label='Enter delta x of first pdb')
    y1 = forms.FloatField(label='Enter delta y of first pdb')
    z1 = forms.FloatField(label='Enter delta z of first pdb')
    temperature = forms.FloatField(label='Enter temperature')
