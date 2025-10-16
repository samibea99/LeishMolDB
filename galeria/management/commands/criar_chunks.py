# galeria/management/commands/criar_chunks.py
import os
import math
import shutil
from django.core.management.base import BaseCommand
from django.conf import settings

class Command(BaseCommand):
    help = 'Cria os arquivos de chunk da biblioteca a partir dos MOL2 individuais.'

    def handle(self, *args, **kwargs):
        self.stdout.write(self.style.SUCCESS("--- Iniciando a criação dos chunks da biblioteca ---"))

        # Define os caminhos de entrada (onde estão os MOL2 individuais) e saída
        input_dir = os.path.join(settings.BASE_DIR, 'mol2_convertidos')
        output_dir = os.path.join(settings.BASE_DIR, 'data', 'library_chunks')

        if not os.path.exists(input_dir):
            self.stdout.write(self.style.ERROR(f"ERRO: A pasta de entrada '{input_dir}' não foi encontrada."))
            self.stdout.write(self.style.WARNING("Por favor, certifique-se de que a pasta 'mol2_convertidos' que você baixou do Colab está na raiz do seu projeto."))
            return

        # Limpa e cria o diretório de saída para garantir que não haja lixo antigo
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)
        os.makedirs(output_dir)
        self.stdout.write(f"Diretório de saída limpo e criado em: {output_dir}")

        # Lista todos os arquivos .mol2 que foram convertidos com sucesso
        mol2_files = sorted([f for f in os.listdir(input_dir) if f.endswith('.mol2')])
        total_files = len(mol2_files)

        if total_files == 0:
            self.stdout.write(self.style.WARNING("Nenhum arquivo .mol2 encontrado na pasta de entrada."))
            return

        self.stdout.write(f"Encontrados {total_files} arquivos .mol2 para agrupar.")

        # Define o número de chunks que queremos criar (8 é um bom número para paralelismo)
        num_chunks = 8
        chunk_size = math.ceil(total_files / num_chunks)

        # Itera de 0 a 7 para criar os 8 chunks
        for i in range(num_chunks):
            start_index = i * chunk_size
            end_index = start_index + chunk_size
            files_in_chunk = mol2_files[start_index:end_index]

            if not files_in_chunk:
                continue

            chunk_filename = os.path.join(output_dir, f'database_chunk_{i}.mol2')
            self.stdout.write(f"Criando {os.path.basename(chunk_filename)} com {len(files_in_chunk)} moléculas...")

            # Abre o arquivo de chunk para escrita
            with open(chunk_filename, 'w', encoding='utf-8') as outfile:
                # Itera sobre cada arquivo individual que pertence a este chunk
                for filename in files_in_chunk:
                    filepath = os.path.join(input_dir, filename)
                    with open(filepath, 'r', encoding='utf-8') as infile:
                        # Copia o conteúdo do arquivo individual para o arquivo de chunk
                        outfile.write(infile.read())
                    # Garante que haja uma linha em branco entre as moléculas
                    outfile.write("\n")

        self.stdout.write(self.style.SUCCESS("\nCriação dos chunks da biblioteca concluída com sucesso!"))
