'''
This is the main application functions. Read the user's upload files (fastaq format). Group by Genome barcode and align to the reference genome. 
The barcodes file contains coupling barcodes.
'''
from django.apps import AppConfig
import glob
from Bio import SeqIO
import os
from zipfile import ZipFile
import shutil
import zipfile

from django.http import HttpResponse
import numpy as np

class nanodorfConfig(AppConfig):
    name = 'nanodorf'

def handle_uploaded_file(f, s_id, f_name = 'fastq.zip'):
    with open('users_file/%s/%s' %(s_id,f_name), 'wb+') as destination:
        for chunk in f.chunks():
            destination.write(chunk)

def manage_fastq_list(s_id):
    '''
    Read user upload files in zip format and return the list of fastq entries for later process, pass it in to YAML file
    Run Quality
    '''
    zipfilename = 'users_file/%s/fastq.zip' %s_id 
    file1 = ZipFile(zipfilename)
    os.mkdir('users_file/%s/fastq'%s_id)
    os.mkdir('users_file/%s/Analysis/'%s_id)
    os.mkdir('users_file/%s/Analysis/Results'%s_id)
    os.mkdir('users_file/%s/Analysis/Results/QC'%s_id)
    path = 'users_file/%s/fastq'%s_id
    file1.extractall(path=path)
    samples_data = os.listdir('users_file/%s/fastq'%s_id)
    for i, each_group in enumerate(samples_data):
        os.system('fastqc -o users_file/%s/Analysis/Results/QC/ users_file/%s/fastq/%s/*.fastq.gz'%(s_id, s_id, each_group))
        stored_group = {}
        each_group_list = os.listdir('users_file/%s/fastq/%s'%(s_id, each_group))
        group1 = {}
        for group in each_group_list:
            group1[group.split('.')[0]] = 'fastq/%s/%s'%(each_group, group)
        stored_group[each_group] = group1
        samples_data[i] = stored_group
    print(samples_data)
    os.unlink(zipfilename)

    return samples_data    


def create_yaml(s_id, samples, yaml_file = 'config.yaml', ref_group = 0, readCountMinThreshold = 10, lfcThreshold =1,  adjPValueThreshold = 0.05, tutorialText=False, organism='human'):
    # '''
    # Create the user defined YAML from YAML template for RScript. The fastq files are separate into different groups as subfolders. Named by directory's name and file'sname.
    # '''
    import yaml   
    genome_annotation = {"human":"ftp://ftp.ensembl.org/pub/release-102/gtf/homo_sapiens/Homo_sapiens.GRCh38.102.gtf.gz",\
                        "rat": "ftp://ftp.ensembl.org/pub/release-102/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.102.gtf.gz",\
                        "mouse": "ftp://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm39.102.gtf.gz",\
                        "zebrafish": "ftp://ftp.ensembl.org/pub/release-102/gtf/danio_rerio/Danio_rerio.GRCz11.102.gtf.gz",\
                        "celegans": "ftp://ftp.ensembl.org/pub/release-102/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.102.gtf.gz"}

    with open(yaml_file) as file:
        default_config = yaml.load(file, Loader=yaml.FullLoader)
    default_config['readCountMinThreshold'] = readCountMinThreshold
    default_config['lfcThreshold'] = lfcThreshold
    default_config['adjPValueThreshold'] = adjPValueThreshold
    default_config['tutorialText'] = tutorialText
    default_config['Samples'] = samples
    default_config['referenceGroup'] = ref_group
    default_config['organism'] = organism
    default_config['genome_annotation'] = genome_annotation[organism]


    with open('users_file/%s/config.yaml'%s_id, 'w') as f:
        yaml.dump(default_config, f)
    os.mkdir('users_file/%s/Static'%s_id)
    os.mkdir('users_file/%s/Static/R'%s_id)
    shutil.copyfile('Static/R/common.R', 'users_file/%s/Static/R/common.R'%s_id)

    return


def run_minimap2(path='users_file/', s_id = 'Test_name_1618217069', organism = 'human'):
    if not os.path.exists('users_file/%s/Analysis'%s_id):
        os.mkdir('users_file/%s/Analysis'%s_id)
    path_minimap = 'users_file/%s/Analysis/Minimap'%s_id
    os.mkdir(path_minimap)
    path_flagstat = 'users_file/%s/Analysis/flagstat/'%s_id
    os.mkdir(path_flagstat)
    path1 = '%s%s'%(path,s_id)
    for group in os.listdir('users_file/%s/fastq/'%s_id):
        for fastq_file in os.listdir('users_file/%s/fastq/%s'%(s_id, group)):
            path2 = '%s/fastq/%s/%s'%(path1, group, fastq_file)
            fastq_file1 = fastq_file.split('.')[0]
            os.system('minimap2 -t16 -ax splice -k14 -secondary=no /home/ag-rossi/ReferenceData/reference_%s.mmi %s | samtools view -Sb | samtools sort - -o %s/%s.bam'%(organism, path2, path_minimap, fastq_file1))
            os.system('samtools flagstat %s/%s.bam>%s%s.txt'%(path_minimap, fastq_file1, path_flagstat, fastq_file1))
    return

def write_rscript(path='users_file/', s_id = 'Test_name_1618217069/'):
    new_R = 'setwd("/home/ag-rossi/projects/nanodorf/nanodorf/%s%s")\n'%(path, s_id)
    new_R += open('RNA.R', 'r').read()

    f = open(path+s_id+'/RNA.R', 'w')
    f.write(new_R)
    f.close()
    return
        
# 
def send_result(submission_name, link_address, recipient_email):
    from django.core.mail import send_mail
    send_mail(subject='Your RNAseq computational results',
        message= 'Dear User\nYour RNAseq analyse <%s> job is now completed\n\Your computational result will be stored in 30 days, You can download your result from this link <<%s>>\n\n\
    Thanks for using our service, if you have any problems please send an email to the developer webserver.iuf@gmail.com'%(submission_name, link_address),
        from_email=None,
        recipient_list=[recipient_email],
        fail_silently=False,
    )
    return

# Start WT and genome alignment
# Require pyensembl, minimap2
# pyensembl install --release 102 --species homo_sapiens

def get_genome_data(gene_name = 'ACE2'):
    from pyensembl import EnsemblRelease
    release = EnsemblRelease()[0]

def pos_process(xls_file = 'test11.xls'):
    import pandas as pd
    data = pd.read_excel('ExpressedGenes.xlsx')
    for col in ['barcode01', 'barcode02','barcode03','barcode04','barcode05','barcode06']:
        data[col] = data[col]/data[col].sum()
    data.to_excel('ExpressedGenes1.xlsx', index=False)

    return 

