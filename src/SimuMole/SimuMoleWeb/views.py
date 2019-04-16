from django.shortcuts import render
from django.http import HttpResponse
from django.http import Http404
from django.template import loader

from .forms import CreateSimulationForm
from SimuMoleScripts.simulation_main_script import Simulation


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
