# Generated by Django 4.2.7 on 2024-03-09 16:15

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('galeria', '0014_alter_similaridade_img'),
    ]

    operations = [
        migrations.AlterField(
            model_name='similaridade',
            name='sdf',
            field=models.FileField(upload_to='sdfs/'),
        ),
    ]
