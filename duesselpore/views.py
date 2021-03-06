from django.shortcuts import render
from django.http import HttpResponse
from .models import InputForm
from .apps import *
import os
from django.conf import settings
import requests
from django.contrib import messages
import sys

from django.views.generic.edit import FormView
from django import forms

import mimetypes
from time import time
from distutils.dir_util import copy_tree
import shutil

def index(request):
    os.chdir(os.path.dirname(__file__))
    if request.method == 'POST':
        form = InputForm(request.POST, request.FILES)
        if form.is_valid():
            if not request.FILES =={}:
                t1 = time()
                form = form.save(commit=False)
                session_id = '%s_%s' %(form.name, str(int(time())))
                organism = form.reference_genes
                os.mkdir('users_file/%s' %session_id)

                print('Upload time')

                stdoutOrigin=sys.stdout 
                os.mkdir('users_file/%s/fastq'%session_id)
                os.mkdir('users_file/%s/Analysis/'%session_id)
                os.mkdir('users_file/%s/Analysis/Results'%session_id)
                os.mkdir('users_file/%s/Analysis/Results/QC'%session_id)
                sys.stdout = open('users_file/%s/Analysis/Results/log.txt' %session_id, "w")
                handle_uploaded_file(request.FILES['upfile_fastq'], s_id=session_id)

                t2 = time() 
                print('Upload time %i seconds' %(t2- t1))

                samples = manage_fastq_list(session_id)
                t3 = time() 
                print('Read time %i seconds' %(t3 - t2))
                
                create_yaml(s_id=session_id, samples=samples, NumberOfTopGene=form.NumberOfTopGene ,ref_group= form.reference_group, study_group=form.study_group, readCountMinThreshold=form.readCountMinThreshold, lfcThreshold=form.lfcThreshold, adjPValueThreshold=form.adjPValueThreshold, organism=organism, cluster_col = form.cluster_by_replica, pathway_ID=form.pathway_ID)

                #Run minimap for two first options
                if form.gene_count_method in ['Rsubread', 'HTSeq']:
                    run_minimap2(s_id=session_id, seq_method=form.seq_method, n_cpu=form.number_cpu)
                    t4 = time()
                    print('Run Minimap time %i seconds' %(t4-t3))

                    if form.gene_count_method == 'Rsubread':
                        print('Starting R Analyser Feature counts')
                        write_rscript(s_id=session_id, method = 'Rsubread', DESeq=form.Differential_expression_method)
                        os.system('R < users_file/%s/RNA.R --no-save'%session_id)
                        t5= time()
                        print('R Runtime Feature counts in %i seconds'%(t5-t4))

                    elif form.gene_count_method == 'HTSeq':
                        print('Starting HTSeq counts')
                        run_htseq_count_parallel(s_id=session_id)
                        write_rscript(s_id=session_id, method = 'HTSeq', DESeq=form.Differential_expression_method)
                        os.system('R < users_file/%s/RNA.R --no-save'%session_id)
                        t5= time()
                        print('R Runtime HTSeq in %i seconds'%(t5-t4))
                    processing_gene_counts(excel_file='users_file/%s/Analysis/Results/ExpressedGenes.xlsx' %session_id)

                elif form.gene_count_method == 'Salmon':
                    run_minimap2_transcriptome(s_id=session_id)
                    t4 = time()
                    print('Run Minimap time %i seconds' %(t4-t3))
                    print('Starting Salmon counts')
                    run_salmon_count(s_id=session_id)

            
                sys.stdout.close()
                sys.stdout=stdoutOrigin
                copy_tree('users_file/%s/Analysis/Results/'%session_id, 'static/results/%s'%session_id)
                shutil.copyfile('templates/report.html', 'static/results/%s/report.html'%session_id)
                shutil.copyfile('templates/report.html', 'users_file/%s/Analysis/Results/report.html'%session_id)
                shutil.make_archive('static/results/%s'%session_id, 'zip', 'users_file/%s/Analysis/Results/' %session_id)
                try:
                    shutil.rmtree('users_file/%s/fastq/'%session_id)
                except OSError as e:
                    print ("Error: %s - %s." % (e.filename, e.strerror))
                
                context = {'file_id': session_id}
                #link = 'http://172.17.21.81:8000/static/%s.zip'%session_id
                # send email is not implemented in local mode
                #send_result(str(form.name), link, form.email)
                
            return render(request, 'results_rna.html', context)

    else:
        form = InputForm()

    return render(request, 'input_RNAseq.html', locals())
            # {'form': form})

def processing_gene_counts(excel_file = 'ExpressedGenes.xlsx'):
    '''
    Fixed missing Genename from Bioconductor
    '''
    from pyensembl import EnsemblRelease
    import pandas as pd

    gene_list = EnsemblRelease(102)
    new_list = []
    countmat2 = pd.read_excel(excel_file)
    for id1 in countmat2.gene_id:
        try:
            new_list.append(gene_list.gene_by_id(id1).gene_name)
        except Exception:
            new_list.append('not found')
    countmat2['gene_name_full'] = new_list
    countmat2.to_excel(excel_file)

    return


def download_file(request):
    print(os.getcwd())
    fl_path = 'duesselpore/static/test_result/test.zip'
    filename = 'downloaded_file_name.extension'

    fl = open(fl_path, 'r')
    mime_type, _ = mimetypes.guess_type(fl_path)
    response = HttpResponse(fl, content_type=mime_type)
    
    response['Content-Disposition'] = "attachment; filename=%s" % filename
    return response
