# Generated by Django 4.2.7 on 2023-12-03 18:04

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('galeria', '0008_similaridade_url'),
    ]

    operations = [
        migrations.AddField(
            model_name='similaridade',
            name='nome',
            field=models.CharField(default='', max_length=20),
        ),
    ]