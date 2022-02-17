from django import forms

class InputForm(forms.Form):
    files_field = forms.FileField(widget=forms.ClearableFileInput(attrs={'multiple': True}))
    attrs={'onclick': 'numberWithCommas()',}
    attrs={'onclick': 'Start()',}
    barcode = forms.FileField()
    



