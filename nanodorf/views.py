from django.shortcuts import render
from django.http import HttpResponse
from .models import InputForm
from .apps import *
import os
from django.conf import settings
import requests
from django.contrib import messages

from django.views.generic.edit import FormView
from django import forms

import mimetypes
from time import time

def index(request):
    os.chdir(os.path.dirname(__file__))
    print(os.path.dirname(__file__))
    if request.method == 'POST':
        form = InputForm(request.POST, request.FILES)
        if form.is_valid():
            if not request.FILES =={}:
                t1 = time()
                form = form.save(commit=False)
                session_id = '%s_%s' %(form.name, str(int(time())))
                os.mkdir('users_file/%s' %session_id)
                handle_uploaded_file(request.FILES['upfile_fastq'], s_id=session_id)
                print(os.getcwd())
                t2 = time() 
                print('Upload time %i seconds' %(t2- t1))
                samples = manage_fastq_list(session_id)
                t3 = time() 
                print('Read time %i seconds' %(t3 - t2))
                create_yaml(s_id=session_id, samples=samples, ref_group= form.reference_group, readCountMinThreshold=form.readCountMinThreshold, lfcThreshold=form.lfcThreshold, adjPValueThreshold=form.adjPValueThreshold)
                write_rscript(s_id=session_id)
                
                #os.unlink('users_file/'+session_id)
                run_minimap2(s_id=session_id)
                t4 = time()
                print('Run Minimap time %i seconds' %(t4-t3))
                print('Start R Analyser')
                os.system('R < users_file/%s/RNA.R --no-save'%session_id)
                t5= time()
                print('R Runtime %i seconds'%(t5-t4))
                # shutil.copyfile('users_file/%s/Rplots.pdf', 'users_file/%s/Analysis/Results/Rplots.pdf'%session_id)
                shutil.make_archive('static/%s'%session_id, 'zip', 'users_file/%s/Analysis/Results/' %session_id)
                context = {'file_id': session_id}
                link = 'http://172.17.21.81:8000/static/%s.zip'%session_id
                send_result(str(form.name), link, form.email)
                
            return render(request, 'results_rna.html', context)

    else:
        form = InputForm()

    return render(request, 'input_RNAseq.html',
            {'form': form })


def download_file(request):
    print(os.getcwd())
    fl_path = 'nanodorf/static/test_result/test.zip'
    filename = 'downloaded_file_name.extension'

    fl = open(fl_path, 'r')
    mime_type, _ = mimetypes.guess_type(fl_path)
    response = HttpResponse(fl, content_type=mime_type)
    
    response['Content-Disposition'] = "attachment; filename=%s" % filename
    return response
