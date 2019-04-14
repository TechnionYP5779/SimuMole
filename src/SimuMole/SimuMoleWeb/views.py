from django.shortcuts import render
from django.http import HttpResponse
from django.http import Http404


def home(request):
    some_dict = {}
    return render(request, 'home.html', some_dict)
