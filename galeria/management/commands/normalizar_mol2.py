# galeria/management/commands/normalizar_mol2.py
import os
from django.core.management.base import BaseCommand
from django.conf import settings

class Command(BaseCommand):
    help = 'Garante que todos os arquivos .mol2 na pasta de convertidos terminem com uma nova linha.'

    def handle(self, *args, **kwargs):
        input_dir = os.path.join(settings.BASE_DIR, 'mol2_convertidos')
        
        if not os.path.exists(input_dir):
            self.stdout.write(self.style.ERROR(f"ERRO: A pasta '{input_dir}' não foi encontrada."))
            return

        self.stdout.write(self.style.SUCCESS(f"Verificando e normalizando arquivos em '{input_dir}'..."))
        
        mol2_files = [f for f in os.listdir(input_dir) if f.endswith('.mol2')]
        files_fixed = 0

        for filename in mol2_files:
            filepath = os.path.join(input_dir, filename)
            with open(filepath, 'r+') as f:
                # Lê todo o conteúdo
                content = f.read()
                # Se o conteúdo não estiver vazio e não terminar com nova linha
                if content and not content.endswith('\n'):
                    # Volta para o final do arquivo e adiciona uma nova linha
                    f.seek(0, os.SEEK_END)
                    f.write('\n')
                    files_fixed += 1
        
        self.stdout.write(self.style.SUCCESS(f"Verificação concluída. {files_fixed} arquivos foram corrigidos."))