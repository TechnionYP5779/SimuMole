from django.contrib import admin
from django.urls import include, path

from SimuMoleWeb import views

from SimuMoleWeb.forms import SimulationForm0_LoadPdb, \
    SimulationForm1_DetermineRelativePosition, \
    SimulationForm2_SimulationParameters
from SimuMoleWeb.views import SimulationWizard, show_form1

urlpatterns = [
    path('', views.home, name='home'),
    path('admin/', admin.site.urls),
    path('create_simulation/',
         SimulationWizard.as_view([SimulationForm0_LoadPdb,
                                   SimulationForm1_DetermineRelativePosition,
                                   SimulationForm2_SimulationParameters],
                                  condition_dict={'1': show_form1})),
]
