from django.urls import include, path
from duesselpore import views
from . import views
app_name = 'duesselpore'

urlpatterns = [
    path('', views.index, name='index'),
    path('result', views.download_file),
]
