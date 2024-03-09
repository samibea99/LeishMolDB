from django.urls import path
from django.conf import settings
from django.conf.urls.static import static
from galeria.views import index, moleculas, imagem, buscar, ordenar_moleculas_view
from . import views

urlpatterns = [
    path('', views.index, name="index"),
    path('comparar/', views.comparar, name ="comparar"),
    path('moleculas/', views.moleculas, name="moleculas"),
    path('imagem/<int:id>', views.molecula_view, name="imagem"),
    path('buscar/', views.buscar, name='buscar'),
    path('api/get_molecule/', views.get_molecule_data, name='get_molecule_data'),
    path('download/<str:filename>/', views.get_SDF_files, name='download_file'),
    path('resultado_similaridade/', views.resultado_similaridade, name='resultado_similaridade'),
    path('download_pdf/<int:id>/', views.download_pdf, name='download_pdf'),
    path('moleculas/ordenar/', ordenar_moleculas_view, name='ordenar_moleculas'),
    path('api/get_mol_data/<int:id>/', views.get_mol_data, name='api_get_mol'),
    path('teste/', views.exemplo_molecula_view, name='teste'),
] 