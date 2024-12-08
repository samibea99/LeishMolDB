from django.conf import settings
from django.conf.urls.static import static
from django.contrib import admin
from django.urls import path, include

urlpatterns = [
    path('LeishMolDB/admin/', admin.site.urls),
    path('LeishMolDB/', include('galeria.urls')),  # Prefixando a URL com /LeishMolDB/
]

# Servindo arquivos estáticos
urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)

# Servindo arquivos de mídia
urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
