# Generated by Django 4.2.7 on 2023-12-02 18:17

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('galeria', '0006_remove_similaridade_log_and_more'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='similaridade',
            name='nome',
        ),
        migrations.AlterField(
            model_name='similaridade',
            name='sdf',
            field=models.CharField(max_length=255),
        ),
        migrations.AlterField(
            model_name='similaridade',
            name='smile',
            field=models.CharField(max_length=100),
        ),
    ]