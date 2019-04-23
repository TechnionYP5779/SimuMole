from django.contrib import admin
from django.urls import include, path
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
    path('create_simulation_result/', views.create_simulation),
	  path('news/', views.news),
	  path('contact/', views.contact),
	  path('about/', views.about),

    path('create_simulation/',
         SimulationWizard.as_view([SimulationForm0_LoadPdb,
                                   SimulationForm1_DetermineRelativePosition,
                                   SimulationForm2_SimulationParameters],
                                  condition_dict={'1': show_form1})),
    path('upload/', views.file_upload, name='file_upload'),
]+ static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)