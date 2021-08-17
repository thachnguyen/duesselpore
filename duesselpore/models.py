# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.db import models
from django.forms import ModelForm


class Input(models.Model):
    human = 'human'
    zebrafish = 'zebrafish'
    mouse = 'mouse'
    rat = 'rat'
    celegans = 'celegans'

    Reference_gene = [
        (human, 'Human hg38 (Homo sapiens)'),
        (mouse, 'Mouse (Mus musculus)'),
        (rat, 'Rat (Rattus norvegicus)'),
        (zebrafish, 'Zebrafish (Danio rerio)'),
        (celegans, 'Caenorhabditis elegans'),
    ]

    Rsubread = 'Rsubread'
    HTSeq = 'HTSeq'
    Salmon = 'Salmon'

    gene_count_method = [
        (Rsubread, 'Rsubread/featureCounts (Liao et. al. 2014) for gene counts'),
        (HTSeq, 'HTSeq/htseq-counts (Ander et. al. 2014) for gene counts'),
        (Salmon, 'Salmon (Patro et. al. 2017) for transcriptome counts'),
    ]
    
    DESeq2 = 'DESeq2'
    limma = 'limma'
    edgeR = 'edgeR'

    Differential_expression_analysis_method = [
        (DESeq2, 'DESeq2 (Love et. al. 2014)'),
        (limma, 'limma (Ritchie et. al. 2015)'),
        (edgeR, 'edgeR (Robinson et. al. 2010)'),
    ]

    name = models.CharField(max_length=100,  verbose_name='Your submission', default='Test_name')
    upfile_fastq = models.FileField(verbose_name= 'Upload your fastq files (all-in zip format, group by study group required)',  upload_to='users_file/', blank=True, null=True)
    gene_count_method = models.CharField(
        max_length=100,
        choices=gene_count_method,
        default=Rsubread,
    )

    Differential_expression_method = models.CharField(
        max_length=100,
        choices=Differential_expression_analysis_method,
        default=DESeq2,
    )

    # gene_templates = models.CharField(max_length=100,  verbose_name='Reference genome code', default='None')
    reference_genes = models.CharField(
        max_length=100,
        choices=Reference_gene,
        default=human,
    )

    cluster_choices= (('Yes','Yes'), ('No', 'No'))

    reference_group = models.CharField(max_length=100,  verbose_name='Reference group (reference\'s sub-directory name)', default='first group', help_text='Your reference group (must be one of your subfolder name)')
    study_group = models.CharField(max_length=100,  verbose_name='Study group (study groups\'s sub-directory name)', default='second group', help_text='Your study group (must be one of your subfolder name)')
    NumberOfTopGene = models.IntegerField(verbose_name='Number of top variance genes (For Gene Ontology)', default=30)
    readCountMinThreshold = models.IntegerField(verbose_name='readCountMinThreshold (Optional)', default=10)
    lfcThreshold = models.FloatField(verbose_name='lfcThreshold (Optional)', default=1)
    adjPValueThreshold = models.FloatField(verbose_name='adjPValueThreshold (Optional)', default=0.05)
    cluster_by_replica = models.CharField(max_length=6, choices=cluster_choices, verbose_name='Cluster replicate', default='Yes')
    pathway_ID = models.CharField(max_length=10, verbose_name='KEGG id (for pathview)', default='hsa04110')
    # gene_gene_templates_file = models.FileField(verbose_name= 'Upload your reference genome',  upload_to='users_file/', blank=False, null=True)
   
class InputForm(ModelForm):
    class Meta:
        model = Input
        fields = ['name', 'upfile_fastq', 'gene_count_method', 'Differential_expression_method', 'NumberOfTopGene','reference_group', 'study_group', 'reference_genes', 'cluster_by_replica', 'readCountMinThreshold', 'lfcThreshold' , 'adjPValueThreshold', 'pathway_ID']
        
