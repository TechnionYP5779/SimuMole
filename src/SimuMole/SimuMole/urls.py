from django.contrib import admin
from django.urls import include, path

from SimuMoleWeb import views

from SimuMoleWeb.forms import SimulationForm1_LoadPdb, \
    SimulationForm2_DetermineRelativePosition, \
    SimulationForm3_SimulationParameters
from SimuMoleWeb.views import SimulationWizard

urlpatterns = [
    path('', views.home, name='home'),
    path('admin/', admin.site.urls),
    path('create_simulation/', SimulationWizard.as_view([SimulationForm1_LoadPdb,
                                                         SimulationForm2_DetermineRelativePosition,
                                                         SimulationForm3_SimulationParameters])),

]
