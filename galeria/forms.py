from django import forms

class UploadSDFForm(forms.Form):
    file = forms.FileField(label='Selecione um arquivo .sdf')
