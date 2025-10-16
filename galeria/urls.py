from django.urls import path
from django.conf import settings
from django.conf.urls.static import static
from galeria.views import index, moleculas, imagem, buscar, ordenar_moleculas_view, comparar3d, resultado_3d, molecules_graph, molblock_from_smiles, get_alignment_data, molecule_2d_image, download_conformer, pharmacophore_search
from . import views


urlpatterns = [
    path('', views.index, name="index"),
    path('comparar/', views.comparar, name ="comparar"),
    path('moleculas/', views.moleculas, name="moleculas"),
    path('imagem/<int:id>', views.molecula_view, name="imagem"),
    path('buscar/', views.buscar, name='buscar'),
    path('api/get_molecule/', views.get_molecule_data, name='get_molecule_data'),
    path('resultado_similaridade/', views.resultado_similaridade, name='resultado_similaridade'),
    path('download_sdfs/', views.download_sdfs, name='download_sdfs'),
    path('moleculas/ordenar/', ordenar_moleculas_view, name='ordenar_moleculas'),
    path('api/get_mol_data/<int:id>/', views.get_mol_data, name='api_get_mol'),
    path('comparar3d/', views.comparar3d, name='comparar3d'),
    path('resultado_3d/', views.resultado_3d, name='resultado_3d'),
    path('api/molecules/graph/', views.molecules_graph, name='molecules_graph'),
    path('api/molecules/molblock/', views.molblock_from_smiles, name='molblock_from_smiles'),
    path('api/get_alignment/', views.get_alignment_data, name='get_alignment_data'),
    path('molecule_2d/<int:mol_id>/', views.molecule_2d_image, name='molecule_2d_image'),
    path('download_conformer/', views.download_conformer, name='download_conformer'),
    path('search/pharmacophore/', views.pharmacophore_search, name='pharmacophore_search'),
] 