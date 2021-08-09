'''
This is the main application functions. Read the user's upload files (fastaq format).'
'''
from django.apps import AppConfig
import os
from zipfile import ZipFile
import shutil
import yaml
import django

class duesselporeConfig(AppConfig):
    name = 'duesselpore'
    #default_auto_field = django.db.models.BigAutoField

def handle_uploaded_file(f, s_id, f_name = 'fastq1'):
    with open('users_file/%s/%s' %(s_id,f_name), 'wb+') as destination:
        for chunk in f.chunks():
            destination.write(chunk)

def manage_fastq_list(s_id):
    '''
    Read user upload files in zip format and return the list of fastq entries for later process, pass it in to YAML file
    Run Quality
    '''
    zipfilename = 'users_file/%s/fastq1' %s_id 
    #file1 = ZipFile(zipfilename)
    path = 'users_file/%s/fastq'%s_id
    #file1.extractall(path=path)
    os.system('7z x users_file/%s/fastq1 -ousers_file/%s/fastq'%(s_id,s_id))
    samples_data = os.listdir('users_file/%s/fastq'%s_id)
    for i, each_group in enumerate(samples_data):
        os.system('fastqc -o users_file/%s/Analysis/Results/QC/ users_file/%s/fastq/%s/*'%(s_id, s_id, each_group))
        stored_group = {}
        each_group_list = os.listdir('users_file/%s/fastq/%s'%(s_id, each_group))
        group1 = {}
        for group in each_group_list:
            group1[group.split('.')[0]] = 'fastq/%s/%s'%(each_group, group)
        stored_group[each_group] = group1
        samples_data[i] = stored_group
    #print(samples_data)
    os.unlink(zipfilename)

    return samples_data    

def create_yaml(s_id, samples, yaml_file = 'config.yaml', ref_group = 0, readCountMinThreshold = 10, lfcThreshold =1,  adjPValueThreshold = 0.05, tutorialText=False, organism='human', cluster_col = False):
    # '''
    # Create the user defined YAML from YAML template for RScript. The fastq files are separate into different groups as subfolders. Named by directory's name and file'sname.
    # '''   
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
    if cluster_col == 'No':
        default_config['cluster_col'] = False
    else:
        default_config['cluster_col'] = True 

    with open('users_file/%s/config.yaml'%s_id, 'w') as f:
        yaml.dump(default_config, f)
    os.mkdir('users_file/%s/Static'%s_id)
    os.mkdir('users_file/%s/Static/R'%s_id)
    shutil.copyfile('Static/R/common.R', 'users_file/%s/Static/R/common.R'%s_id)

    return

def run_minimap2(path='users_file/', s_id = 'Test_name_1618217069', organism = 'human'):
    file_org={'human':'Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz',
                'rat': 'Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz',
                'mouse':'Mus_musculus.GRCm38.dna.primary_assembly.fa.gz',
                'zebrafish':'Danio_rerio.GRCz11.dna.primary_assembly.fa.gz',
                'celegans':'Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz'}

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
            #os.system('minimap2 -t 16 -a -x map-ont --splice -k 15 -w 10 --secondary=no /home/ag-rossi/ReferenceData/reference_%s.mmi %s | samtools view -Sb | samtools sort - -o %s/%s.bam'%(organism, path2, path_minimap, fastq_file1))
            #os.system('minimap2 -t 16 -ax map-ont --splice --secondary=no /home/ag-rossi/ReferenceData/%s %s | samtools view -Sb | samtools sort - -o %s/%s.bam'%(file_org[organism], path2, path_minimap, fastq_file1))
            os.system('minimap2 -t 16 -ax splice -k14 -uf --secondary=no /home/ag-rossi/ReferenceData/%s %s | samtools view -Sb | samtools sort - -o %s/%s.bam'%(file_org[organism], path2, path_minimap, fastq_file1))
            print('minimap2 -t 16 -ax splice -k14 -uf --secondary=no /home/ag-rossi/ReferenceData/%s %s | samtools view -Sb | samtools sort - -o %s/%s.bam'%(file_org[organism], path2, path_minimap, fastq_file1))
            #os.system('minimap2 -t 16 -ax splice -uf -k14 --secondary=no /home/ag-rossi/ReferenceData/reference.mmi %s | samtools view -Sb | samtools sort - -o %s/%s.bam'%(path2, path_minimap, fastq_file1))
            os.system('samtools flagstat %s/%s.bam>%s%s.txt'%(path_minimap, fastq_file1, path_flagstat, fastq_file1))
    return

def run_minimap2_transcriptome(path='users_file/', s_id = 'Test_name_1618217069', organism = 'human'):
    file_org={'human':'Homo_sapiens.GRCh38.cdna.all.fa.gz',
                'rat': 'Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz',
                'mouse':'Mus_musculus.GRCm38.dna.primary_assembly.fa.gz',
                'zebrafish':'Danio_rerio.GRCz11.dna.primary_assembly.fa.gz',
                'celegans':'Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz'}

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
            #os.system('minimap2 -t 16 -a -x map-ont --splice -k 15 -w 10 --secondary=no /home/ag-rossi/ReferenceData/reference_%s.mmi %s | samtools view -Sb | samtools sort - -o %s/%s.bam'%(organism, path2, path_minimap, fastq_file1))
            #os.system('minimap2 -t 16 -ax map-ont --splice --secondary=no /home/ag-rossi/ReferenceData/%s %s | samtools view -Sb | samtools sort - -o %s/%s.bam'%(file_org[organism], path2, path_minimap, fastq_file1))
            os.system('minimap2 -t 16 -ax splice -k14 -uf --secondary=no /home/ag-rossi/ReferenceData/%s %s | samtools view -Sb | samtools sort - -o %s/%s.bam'%(file_org[organism], path2, path_minimap, fastq_file1))
            print('Run minimap2 for transcriptome. Organism: %s'%organism)
            print('minimap2 -t 16 -ax splice -k14 -uf --secondary=no /home/ag-rossi/ReferenceData/%s %s | samtools view -Sb | samtools sort - -o %s/%s.bam'%(file_org[organism], path2, path_minimap, fastq_file1))
            #os.system('minimap2 -t 16 -ax splice -uf -k14 --secondary=no /home/ag-rossi/ReferenceData/reference.mmi %s | samtools view -Sb | samtools sort - -o %s/%s.bam'%(path2, path_minimap, fastq_file1))
            os.system('samtools flagstat %s/%s.bam>%s%s.txt'%(path_minimap, fastq_file1, path_flagstat, fastq_file1))
    return

def run_htseq_count(path='users_file/', s_id = 'Test_name_1618217069'):
    import pandas as pd
    import glob
    #bam_path= 'users_file/%s/Analysis/Minimap/'%s_id
    #hts_out_path = 'users_file/%s/Analysis/'%s_id

    for i, bamfile in enumerate(glob.glob('users_file/%s/Analysis/Minimap/*.bam'%s_id)):
        os.system('samtools index %s'%(bamfile))
        os.system('htseq-count -s no -a 5 --nonunique=all %s ~/ReferenceData/Homo_sapiens.GRCh38.104.gtf>%s.csv'%(bamfile, bamfile[:-4]))
        if i ==0:
            df = pd.read_csv('%s.csv'%(bamfile[:-4]), delimiter='\t', header=None)
            df.columns=['gene_id', bamfile[:-4]]
        else:
            df1 = pd.read_csv('%s.csv'%(bamfile[:-4]), delimiter='\t', header=None)
            df[bamfile[:-4]] = df1[1]

    df = df[:-5]
    df.to_excel('users_file/%s/Analysis/Results/ExpressedGenes.xlsx'%s_id)
    return

def run_htseq_count_parallel(path='users_file/', s_id = 'Test_name_1618217069'):
    import pandas as pd
    import glob
    bam_files= glob.glob('users_file/%s/Analysis/Minimap/*.bam'%s_id)
    hts_out_path = 'users_file/%s/Analysis/'%s_id
    data_columns = [s.split('/')[-1][:-4] for s in glob.glob('users_file/%s/Analysis/Minimap/*.bam'%s_id)]
    data_columns.insert(0, 'gene_id')
    for i, bamfile in enumerate(glob.glob('users_file/%s/Analysis/Minimap/*.bam'%s_id)):
        os.system('samtools index %s'%(bamfile))
    os.system('htseq-count -s no -a 5 -n %i --nonunique=all %s ~/ReferenceData/Homo_sapiens.GRCh38.104.gtf>%sHTSeq_counts.csv'%(len(bam_files), " ".join(bam_files), hts_out_path))
    
    df1 = pd.read_csv('%sHTSeq_counts.csv'%hts_out_path, delimiter='\t', header=None)
    df1 = df1[:-5]    
    df1.columns = data_columns
    df1.to_excel('users_file/%s/Analysis/Results/ExpressedGenes.xlsx'%s_id)
    return


def run_salmon_count(path='users_file/', s_id = 'Test_name_1618217069'):
    import pandas as pd
    import glob
    #bam_path= 'users_file/%s/Analysis/Minimap/'%s_id
    #hts_out_path = 'users_file/%s/Analysis/'%s_id

    for i, bamfile in enumerate(glob.glob('users_file/%s/Analysis/Minimap/*.bam'%s_id)):
        os.system('samtools index %s'%(bamfile))
        os.system('htseq-count -s no -a 5 --nonunique=all %s ~/ReferenceData/Homo_sapiens.GRCh38.104.gtf>%s.csv'%(bamfile, bamfile[:-4]))
        if i ==0:
            df = pd.read_csv('%s.csv'%(bamfile[:-4]), delimiter='\t', header=None)
            df.columns=['Geneid', bamfile[:-4]]
        else:
            df1 = pd.read_csv('%s.csv'%(bamfile[:-4]), delimiter='\t', header=None)
            df[bamfile[:-4]] = df1[1]

    df = df[:-5]
    df.to_excel('users_file/%s/Analysis/Results/Expressed_gene_HTS.xlxs'%s_id)
    return


def write_rscript(path='users_file/', s_id = 'Test_name_1618217069/', method= 'Rsubread'):
    new_R = 'setwd("/home/ag-rossi/projects/duesselpore/duesselpore/%s%s")\n'%(path, s_id)
    if method == 'Rsubread':
        new_R += open('R/RNA1.R', 'r').read()
        new_R += open('R/RNA2_subread.R', 'r').read()
        new_R += open('R/RNA3.R', 'r').read()
        f = open(path+s_id+'/RNA.R', 'w')
        f.write(new_R)
        f.close()

    elif method == 'HTSeq':
        new_R += open('R/RNA1.R', 'r').read()
        new_R += open('R/RNA2_HTSeq.R', 'r').read()
        new_R += open('R/RNA3.R', 'r').read()
        f = open(path+s_id+'/RNA.R', 'w')
        f.write(new_R)
        f.close()

    elif method == 'Salmon':
        new_R += open('RNA.R', 'r').read()
        f = open(path+s_id+'/RNA.R', 'w')
        f.write(new_R)
        f.close()
    else:
        print('Method is not support')
    return
        
