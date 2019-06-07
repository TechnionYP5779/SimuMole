# Create your models here.
from django.db import models
from django.forms import ModelForm
from django.core import validators


class Upload(models.Model):
    file = models.FileField(upload_to="files/", validators=[validators.FileExtensionValidator(['pdb'])])
    object_name = models.CharField(max_length=200)  # is the file object name, not the file name itself!
    upload_date = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return self.object_name


# FileUpload form class.
class UploadForm(ModelForm):
    class Meta:
        model = Upload
        fields = ('file', 'object_name',)
