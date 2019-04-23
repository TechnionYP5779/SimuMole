from django.shortcuts import render
from django.http import HttpResponse
from django.http import Http404
from django.template import loader
from django.shortcuts import render

from SimuMoleScripts.simulation_main_script import Simulation
from .models import UploadForm, Upload
from django.http import HttpResponseRedirect
from django.urls import reverse
from django.contrib import messages
import os

from .forms import CreateSimulationForm


def home(request):
    some_dict = {}
    return render(request, 'home.html', some_dict)


def create_simulation(request):
    if request.method == 'POST':
        form = CreateSimulationForm(request.POST)

        if form.is_valid():
            first_pdb = form.cleaned_data['first_pdb']
            second_pdb = form.cleaned_data['second_pdb']
            x1 = form.cleaned_data['x1']
            y1 = form.cleaned_data['y1']
            z1 = form.cleaned_data['z1']
            temperature = form.cleaned_data['temperature']

            s = Simulation(first_pdb, second_pdb, x1, y1, z1, temperature)
            s.create_simulation()

            context = {
                'first_pdb': first_pdb,
                'second_pdb': second_pdb,
                'x1': x1,
                'y1': y1,
                'z1': z1,
                'temperature': temperature,
            }

            template = loader.get_template('create_simulation_result.html')
            return HttpResponse(template.render(context, request))

        form = CreateSimulationForm()

    else:
        form = CreateSimulationForm()

    return render(request, 'create_simulation.html', {'form': form})

def file_upload(request):
    if request.method=="POST":
        file = UploadForm(request.POST, request.FILES)
        if file.is_valid():
            file_name, file_extension = os.path.splitext(request.FILES['file'].name)
            file_extension = file_extension.lower()
            # allowing only pdb files
            if file_extension == '.pdb':
                messages.success(request, 'File Uploaded Successfully')
                file.save()
                return HttpResponseRedirect(reverse('file_upload'))
            else:  # Should never be called, since we added FileExtensionValidator on the Upload model.
                messages.error(request, 'Upload failed: file extension has to be \'pdb\'.')

    else:
        file=UploadForm()
    return render(request,'file_upload.html',{'form':file})