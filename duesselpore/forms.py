from django import forms

class InputForm(forms.Form):
    files_field = forms.FileField(widget=forms.ClearableFileInput(attrs={'multiple': True}))
    barcode = forms.FileField()
    



