# Generated by Django 4.2.7 on 2023-12-02 14:00

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('galeria', '0004_alter_similaridade_sdf'),
    ]

    operations = [
        migrations.DeleteModel(
            name='ArquivoSDF',
        ),
        migrations.RemoveField(
            model_name='similaridade',
            name='sdf',
        ),
        migrations.AddField(
            model_name='similaridade',
            name='log',
            field=models.FloatField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name='similaridade',
            name='polaridade',
            field=models.FloatField(blank=True, null=True),
        ),
    ]
