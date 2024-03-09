from django.db import models
from datetime import datetime


class similaridade(models.Model):

    OPCOES_CATEGORIA = [
        ("SYNTHETIC", "Synthetic"),
        ("NATURAL", "Natural"),
        ("DRUG", "Drug"),
    ]

    id = models.BigAutoField(auto_created=True, primary_key=True) 
    nome = models.CharField(max_length=20, blank=False, default='')
    smile = models.CharField(max_length=100)
    categoria = models.CharField(max_length=100, choices=OPCOES_CATEGORIA, default='')
    sdf = models.FileField(upload_to='sdfs/')
    url = models.URLField(blank=True, null=True)
    observacoes_admin = models.TextField("Observações do Administrador", blank=True, null=True)
    publicada = models.BooleanField(default=False)
    data_mol = models.DateTimeField(default=datetime.now, blank=False)
    

    def __str__(self):
        return self.nome

 
       