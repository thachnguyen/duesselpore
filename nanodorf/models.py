# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.db import models
from django.forms import ModelForm
from captcha.fields import ReCaptchaField
# from .forms import FileFieldForm

class Input(models.Model):
    Human_hg38 = 'hg38'
    Zebra_fish = 'zbrf'
    Mouse = 'mouse'
    Rat = 'rat'
    C_elegant = 'rat'

    Reference_gene = [
        (Human_hg38, 'Human hg38'),
        (Zebra_fish, 'Zebrafish'),
        (Mouse, 'Mouse'),
        (Rat, 'Rat'),
    ]

    name = models.CharField(max_length=100,  verbose_name='Your submission', default='Test_name')
    upfile_fastq5 = models.FileField(verbose_name= 'Upload your fast5 files (very slow, not recommend)',  upload_to='users_file/', blank=True, null=True)
    upfile_fastq = models.FileField(verbose_name= 'Upload your fastq files (all-in zip format, group by study group required)',  upload_to='users_file/', blank=True, null=True)
    # gene_templates = models.CharField(max_length=100,  verbose_name='Reference genome code', default='None')
    reference_genes = models.CharField(
        max_length=100,
        choices=Reference_gene,
        default=Human_hg38,
    )
    reference_group = models.CharField(max_length=100,  verbose_name='your reference group (reference\'s directory name)', default='group01')
    readCountMinThreshold = models.IntegerField(verbose_name='readCountMinThreshold (Optional)', default=10)
    lfcThreshold = models.FloatField(verbose_name='lfcThreshold (Optional)', default=1)
    adjPValueThreshold = models.FloatField(verbose_name='adjPValueThreshold (Optional)', default=0.05)
    # gene_gene_templates_file = models.FileField(verbose_name= 'Upload your reference genome',  upload_to='users_file/', blank=False, null=True)
    email = models.EmailField(verbose_name='Email address', default= 'thach.nguyen@iuf-duesseldorf.de')
    
class InputForm(ModelForm):
    # captcha = ReCaptchaField(public_key='6LesYJ4aAAAAAHDch9YM2xb0EeLDqAaT7zKTJZ4g',
    # private_key='6LesYJ4aAAAAAAAaarJlC14fNcBPc9QzMu1FyXpx',)
    # captcha = ReCaptchaField()
    class Meta:
        model = Input
        # fields = ['name','upfile_fastq','upfile_barcode', 'gene_templates', 'gene_gene_templates_file', 'email']
        fields = ['name', 'upfile_fastq', 'upfile_fastq5', 'reference_group','reference_genes', 'readCountMinThreshold', 'lfcThreshold' , 'adjPValueThreshold', 'email']
        
