from django.contrib import admin
from django.urls import include, path

from SimuMoleWeb import views

urlpatterns = [
    path('', views.home, name='home'),
    path('admin/', admin.site.urls),
    path('create_simulation/', views.create_simulation, name='create_simulation'),
    path('create_simulation_result/', views.create_simulation),
	path('news/', views.news),
	path('contact/', views.contact),
	path('about/', views.about),
]
