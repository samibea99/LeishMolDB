# Generated by Django 4.2.7 on 2023-12-03 21:01

import datetime
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('galeria', '0011_similaridade_publicada'),
    ]

    operations = [
        migrations.AddField(
            model_name='similaridade',
            name='data_mol',
            field=models.DateTimeField(default=datetime.datetime.now),
        ),
    ]
