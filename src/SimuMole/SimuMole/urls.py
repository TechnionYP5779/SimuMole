from django.contrib import admin
from django.urls import path
from django.conf import settings
from django.conf.urls.static import static
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

                  path('update_simulation_status/', views.update_simulation_status, name='update_simulation_status'),

                  path('download__zip/', views.download__zip, name='download__zip'),
                  path('download__email/', views.download__email, name='download__email'),

                  path('upload/', views.upload_files, name='upload_files'),

                  path('help/', views.help_page, name='help_page')

              ] + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)