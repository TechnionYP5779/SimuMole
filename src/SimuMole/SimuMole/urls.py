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


                  path('download_pdb_dcd__zip/', views.download_pdb_dcd__zip, name='download_pdb_dcd__zip'),
                  path('download_pdb_dcd__email/', views.download_pdb_dcd__email, name='download_pdb_dcd__email'),
                  path('download_animation__zip/', views.download_animation__zip, name='download_animation__zip'),
                  path('download_animation__email/', views.download_animation__email, name='download_animation__email'),
                  path('upload/', views.file_upload, name='file_upload'),

              ] + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
