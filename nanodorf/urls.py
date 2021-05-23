from django.urls import include, path
from nanodorf import views
from . import views
app_name = 'nanodorf'

urlpatterns = [
    path('', views.index, name='index'),
    path('result', views.download_file),
]
