from django.contrib import admin
from django.urls import include, path

from SimuMoleWeb import views

from SimuMoleWeb.forms import SimulationForm0_NumberOfProteins, \
    SimulationForm1_LoadPdb, \
    SimulationForm2_DetermineRelativePosition, \
    SimulationForm3_SimulationParameters
from SimuMoleWeb.views import SimulationWizard, show_form2

urlpatterns = [
    path('', views.home, name='home'),
    path('admin/', admin.site.urls),
    path('create_simulation/',
         SimulationWizard.as_view([SimulationForm0_NumberOfProteins,
                                   SimulationForm1_LoadPdb,
                                   SimulationForm2_DetermineRelativePosition,
                                   SimulationForm3_SimulationParameters],
                                  condition_dict={'2': show_form2})),
]
