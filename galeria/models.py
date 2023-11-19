from django.db import models

class similaridade(models.Model):
    id = models.BigAutoField(auto_created=True, primary_key=True) 
    smile = models.CharField(max_length=100)
    sdf = models.URLField()

    def __str__(self):
        return f"similaridade [id={self.id}]"
