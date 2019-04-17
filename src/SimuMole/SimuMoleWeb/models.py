from django.db import models

# Create your models here.
from django.db import models
from django.forms import ModelForm

class Upload(models.Model):
    file = models.FileField(upload_to="files/")
    file_name = models.CharField(max_length=200) #is the file object name, not the file name itself!
    upload_date=models.DateTimeField(auto_now_add =True)
    def __str__(self):
        return self.file_name

# FileUpload form class.
class UploadForm(ModelForm):
    class Meta:
        model = Upload
        fields = ('file','file_name',)