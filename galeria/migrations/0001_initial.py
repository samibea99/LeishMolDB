# Generated by Django 4.2.7 on 2023-11-18 23:55

from django.db import migrations, models


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='similaridade',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False)),
                ('smile', models.CharField(max_length=100)),
                ('sdf', models.URLField()),
            ],
        ),
    ]
