# Generated by Django 4.2.7 on 2023-12-02 17:21

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('galeria', '0005_delete_arquivosdf_remove_similaridade_sdf_and_more'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='similaridade',
            name='log',
        ),
        migrations.RemoveField(
            model_name='similaridade',
            name='polaridade',
        ),
        migrations.AddField(
            model_name='similaridade',
            name='nome',
            field=models.CharField(max_length=100, null=True),
        ),
        migrations.AddField(
            model_name='similaridade',
            name='sdf',
            field=models.CharField(default=0.0, max_length=300),
        ),
        migrations.AlterField(
            model_name='similaridade',
            name='smile',
            field=models.CharField(max_length=200),
        ),
    ]
