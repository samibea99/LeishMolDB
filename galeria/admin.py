from django.contrib import admin
from galeria.models import similaridade 

class ListandoSDF(admin.ModelAdmin):
    list_display = ("id", "nome", "smile", "sdf", "url", 'observacoes_admin', "publicada")
    list_display_links = ("id", "smile", "sdf")
    search_fields = ("smile",)
    list_filter = ("categoria",)
    list_editable = ("publicada",)
    list_per_page = 50

admin.site.register(similaridade, ListandoSDF)
