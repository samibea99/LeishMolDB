# galeria/management/commands/exportar_dados.py
import os
import shutil
import pandas as pd
from django.core.management.base import BaseCommand
from django.conf import settings
from galeria.models import similaridade

class Command(BaseCommand):
    help = 'Exporta os arquivos SDF e metadados para processamento no Google Colab.'

    def handle(self, *args, **kwargs):
        self.stdout.write(self.style.SUCCESS("--- Iniciando exportação de dados para o Colab ---"))

        # Cria a pasta principal para a exportação
        export_dir = os.path.join(settings.BASE_DIR, 'dados_para_colab')
        sdfs_dir = os.path.join(export_dir, 'sdfs')
        if os.path.exists(export_dir):
            shutil.rmtree(export_dir)
        os.makedirs(sdfs_dir)
        
        self.stdout.write(f"Criando pastas de exportação em: {export_dir}")

        # Busca as moléculas no banco de dados
        moleculas = similaridade.objects.filter(publicada=True).exclude(sdf__exact='').exclude(sdf__isnull=True)
        
        metadata = []
        self.stdout.write(f"Encontradas {moleculas.count()} moléculas para exportar.")

        for mol in moleculas:
            if mol.sdf and os.path.exists(mol.sdf.path):
                # O novo nome do arquivo será apenas o ID, para criar um link único.
                novo_nome_sdf = f"{mol.id}.sdf"
                caminho_destino = os.path.join(sdfs_dir, novo_nome_sdf)
                
                # Copia o arquivo SDF para a nova pasta com o novo nome
                shutil.copy(mol.sdf.path, caminho_destino)
                
                # Adiciona os metadados à nossa lista
                metadata.append({
                    'id': mol.id,
                    'nome': mol.nome,
                    'categoria': mol.categoria,
                    'observacoes': mol.observacoes_admin,
                    'sdf_original': os.path.basename(mol.sdf.name),
                    'sdf_novo_path': os.path.join('sdfs', novo_nome_sdf) # Caminho relativo dentro do zip
                })
        
        # Cria um DataFrame do pandas com os metadados
        df = pd.DataFrame(metadata)
        metadata_path = os.path.join(export_dir, 'metadata.csv')
        df.to_csv(metadata_path, index=False)
        
        self.stdout.write(self.style.SUCCESS(f"\nExportação concluída!"))
        self.stdout.write(f"Foram exportados {len(metadata)} arquivos SDF.")
        self.stdout.write(f"Arquivo de metadados salvo em: {metadata_path}")
        self.stdout.write("\nPróximo passo: Compacte a pasta 'dados_para_colab' em um arquivo .zip e faça o upload para o Google Colab.")