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
        (Salmon, 'Salmon (Patro et. al. 2017) for gene counts'),
    ]
    

    name = models.CharField(max_length=100,  verbose_name='Your submission', default='Test_name')
    upfile_fastq = models.FileField(verbose_name= 'Upload your fastq files (all-in zip format, group by study group required)',  upload_to='users_file/', blank=True, null=True)
    gene_count_method = models.CharField(
        max_length=100,
        choices=gene_count_method,
        default=Rsubread,
    )
    # gene_templates = models.CharField(max_length=100,  verbose_name='Reference genome code', default='None')
    reference_genes = models.CharField(
        max_length=100,
        choices=Reference_gene,
        default=human,
    )
    reference_group = models.CharField(max_length=100,  verbose_name='your reference group (reference\'s directory name)', default='group01')
    readCountMinThreshold = models.IntegerField(verbose_name='readCountMinThreshold (Optional)', default=10)
    lfcThreshold = models.FloatField(verbose_name='lfcThreshold (Optional)', default=1)
    adjPValueThreshold = models.FloatField(verbose_name='adjPValueThreshold (Optional)', default=0.05)
    cluster_by_replica = models.CharField(max_length=10,  verbose_name='Cluster replicate', default='FALSE')
    # gene_gene_templates_file = models.FileField(verbose_name= 'Upload your reference genome',  upload_to='users_file/', blank=False, null=True)
   
class InputForm(ModelForm):
    class Meta:
        model = Input
        fields = ['name', 'upfile_fastq', 'gene_count_method', 'reference_group','reference_genes', 'cluster_by_replica', 'readCountMinThreshold', 'lfcThreshold' , 'adjPValueThreshold']
        
