# Generated by Django 2.2 on 2019-04-21 12:28

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('SimuMoleWeb', '0001_initial'),
    ]

    operations = [
        migrations.RenameField(
            model_name='upload',
            old_name='file_name',
            new_name='object_name',
        ),
    ]
