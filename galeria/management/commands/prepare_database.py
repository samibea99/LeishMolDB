import os
import logging
import subprocess
import math
import shutil
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed
from django.core.management.base import BaseCommand
from django.conf import settings
from galeria.models import similaridade

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def validate_sdf_file(sdf_path):
    """
    Valida se o arquivo SDF é válido antes da conversão.
    """
    try:
        logging.info(f"[VALIDATE] Iniciando validação do arquivo: {sdf_path}")
        
        # Verifica se o arquivo existe e não está vazio
        if not os.path.exists(sdf_path):
            logging.error(f"[VALIDATE] ERRO: Arquivo não existe: {sdf_path}")
            return False, "Arquivo não existe"
        
        file_size = os.path.getsize(sdf_path)
        if file_size == 0:
            logging.error(f"[VALIDATE] ERRO: Arquivo vazio (0 bytes): {sdf_path}")
            return False, "Arquivo vazio"
        
        logging.info(f"[VALIDATE] Arquivo existe com {file_size} bytes: {sdf_path}")
        
        # Verifica se o arquivo tem pelo menos uma molécula válida
        with open(sdf_path, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read().strip()
            
        logging.info(f"[VALIDATE] Conteúdo lido: {len(content)} caracteres")
        
        # Verifica se contém marcadores básicos de SDF
        if '$$$$' not in content:
            logging.error(f"[VALIDATE] ERRO: Arquivo SDF inválido - sem marcador de fim de molécula ($$$$) em: {sdf_path}")
            logging.error(f"[VALIDATE] Primeiros 200 caracteres do arquivo: {content[:200]}")
            return False, "Arquivo SDF inválido - sem marcador de fim de molécula"
        
        # Verifica se tem pelo menos uma estrutura molecular básica
        if 'V2000' not in content and 'V3000' not in content:
            logging.error(f"[VALIDATE] ERRO: Arquivo SDF inválido - sem cabeçalho de estrutura (V2000/V3000) em: {sdf_path}")
            logging.error(f"[VALIDATE] Primeiros 200 caracteres do arquivo: {content[:200]}")
            return False, "Arquivo SDF inválido - sem cabeçalho de estrutura"
        
        # Conta quantas moléculas existem no arquivo
        mol_count = content.count('$$$$')
        logging.info(f"[VALIDATE] Arquivo SDF válido com {mol_count} molécula(s): {sdf_path}")
        return True, f"Arquivo SDF válido com {mol_count} molécula(s)"
        
    except Exception as e:
        logging.error(f"[VALIDATE] ERRO INESPERADO ao validar {sdf_path}: {e}")
        logging.error(f"[VALIDATE] Tipo do erro: {type(e).__name__}")
        return False, f"Erro ao validar arquivo: {e}"

def convert_sdf_to_mol2_antechamber(sdf_path, mol2_path):
    """
    Converte SDF para MOL2 usando AmberTools (antechamber).
    Gera coordenadas 3D perfeitas e cargas parciais.
    """
    try:
        logging.info(f"[CONVERT] Iniciando conversão: {sdf_path} -> {mol2_path}")
        
        # Primeiro valida o arquivo SDF
        is_valid, validation_msg = validate_sdf_file(sdf_path)
        if not is_valid:
            logging.error(f"[CONVERT] ERRO: Validação falhou para {sdf_path}: {validation_msg}")
            return False
        
        logging.info(f"[CONVERT] Validação passou: {validation_msg}")
        
        # Comando antechamber para conversão SDF -> MOL2 (otimizado)
        comando = [
                    "antechamber",
                    "-i", sdf_path,
                    "-fi", "sdf",
                    "-o", mol2_path,
                    "-fo", "mol2",
            "-c", "bcc",  # Método de cálculo de carga
            "-s", "1",    # Verbosidade reduzida (mais rápido)
            "-at", "gaff", # Force field GAFF
            "-j", "4"     # Usar 4 threads para antechamber
        ]
        
        logging.info(f"[CONVERT] Executando comando: {' '.join(comando)}")
        logging.info(f"[CONVERT] Diretório de trabalho: {os.getcwd()}")
        
        resultado = subprocess.run(comando, check=True, capture_output=True, text=True, timeout=120)  # Timeout reduzido para 2 minutos
        
        logging.info(f"[CONVERT] Comando executado com sucesso. Return code: {resultado.returncode}")
        if resultado.stdout:
            logging.info(f"[CONVERT] Stdout: {resultado.stdout[:500]}...")
        
        # Verifica se o arquivo foi criado e não está vazio
        if not os.path.exists(mol2_path):
            logging.error(f"[CONVERT] ERRO: Arquivo MOL2 não foi criado: {mol2_path}")
            return False
        
        mol2_size = os.path.getsize(mol2_path)
        if mol2_size == 0:
            logging.error(f"[CONVERT] ERRO: Arquivo MOL2 criado mas está vazio (0 bytes): {mol2_path}")
            return False
        
        logging.info(f"[CONVERT] Arquivo MOL2 criado com {mol2_size} bytes: {mol2_path}")
        
        # Verifica se o MOL2 gerado é válido
        with open(mol2_path, 'r', encoding='utf-8') as f:
            mol2_content = f.read()
            
        if '@<TRIPOS>MOLECULE' not in mol2_content:
            logging.error(f"[CONVERT] ERRO: MOL2 inválido - sem seção @<TRIPOS>MOLECULE em: {mol2_path}")
            logging.error(f"[CONVERT] Primeiros 300 caracteres do MOL2: {mol2_content[:300]}")
            return False
            
        if '@<TRIPOS>ATOM' not in mol2_content:
            logging.error(f"[CONVERT] ERRO: MOL2 inválido - sem seção @<TRIPOS>ATOM em: {mol2_path}")
            logging.error(f"[CONVERT] Primeiros 300 caracteres do MOL2: {mol2_content[:300]}")
            return False
        
        logging.info(f"[CONVERT] SUCESSO: Conversão completa e válida: {mol2_path}")
        return True
            
    except subprocess.CalledProcessError as e:
        stderr = e.stderr.decode('utf-8', errors='ignore') if e.stderr and isinstance(e.stderr, bytes) else str(e.stderr) if e.stderr else ""
        stdout = e.stdout.decode('utf-8', errors='ignore') if e.stdout and isinstance(e.stdout, bytes) else str(e.stdout) if e.stdout else ""
        logging.error(f"[CONVERT] ERRO: antechamber falhou para {sdf_path}")
        logging.error(f"[CONVERT] Return code: {e.returncode}")
        logging.error(f"[CONVERT] Comando que falhou: {' '.join(comando)}")
        if stderr:
            logging.error(f"[CONVERT] Stderr: {stderr}")
        if stdout:
            logging.error(f"[CONVERT] Stdout: {stdout}")
        return False
    except subprocess.TimeoutExpired as e:
        logging.error(f"[CONVERT] ERRO: Timeout (120s) na conversão de {sdf_path}")
        logging.error(f"[CONVERT] Comando que teve timeout: {' '.join(comando)}")
        return False
    except Exception as e:
        logging.error(f"[CONVERT] ERRO INESPERADO na conversão de {sdf_path}: {e}")
        logging.error(f"[CONVERT] Tipo do erro: {type(e).__name__}")
        logging.error(f"[CONVERT] Traceback: {e.__traceback__}")
        return False

def process_single_molecule(mol_data):
    """
    Processa uma única molécula. Função para processamento paralelo.
    """
    mol_obj, temp_conversion_dir, chunk_index = mol_data
    temp_mol2_path = os.path.join(temp_conversion_dir, f'temp_{mol_obj.id}.mol2')
    
    try:
        logging.info(f"[PARALLEL] Processando molécula ID: {mol_obj.id}, Nome: {mol_obj.nome}")
        
        # Verifica se a molécula tem arquivo SDF
        if not mol_obj.sdf:
            return {'success': False, 'mol_id': mol_obj.id, 'error': 'sem_sdf', 'message': 'Sem arquivo SDF associado'}
        
        # Verifica se o arquivo SDF existe fisicamente
        if not os.path.exists(mol_obj.sdf.path):
            return {'success': False, 'mol_id': mol_obj.id, 'error': 'arquivo_inexistente', 'message': f'Arquivo não encontrado: {mol_obj.sdf.path}'}
        
        # Verifica se o arquivo não está vazio
        file_size = os.path.getsize(mol_obj.sdf.path)
        if file_size == 0:
            return {'success': False, 'mol_id': mol_obj.id, 'error': 'arquivo_vazio', 'message': 'Arquivo SDF vazio'}
        
        # Converte SDF para MOL2
        if convert_sdf_to_mol2_antechamber(mol_obj.sdf.path, temp_mol2_path):
            # Lê o conteúdo do MOL2 gerado
            with open(temp_mol2_path, 'r', encoding='utf-8') as f:
                content = f.read()
                nome_seguro = mol_obj.nome.replace(" ", "_").replace("/", "_") if mol_obj.nome else f"mol_{mol_obj.id}"
                content = content.replace('@<TRIPOS>MOLECULE', f'@<TRIPOS>MOLECULE\nmol_{mol_obj.id}_{nome_seguro}')
            
            return {
                'success': True, 
                'mol_id': mol_obj.id, 
                'content': content,
                'chunk_index': chunk_index
            }
        else:
            return {'success': False, 'mol_id': mol_obj.id, 'error': 'conversao_falhou', 'message': 'Falha na conversão SDF->MOL2'}
    
    except Exception as e:
        return {'success': False, 'mol_id': mol_obj.id, 'error': 'erro_inesperado', 'message': str(e)}

class Command(BaseCommand):
    help = 'Processa a biblioteca de moléculas com uma estratégia híbrida, dividindo-a em chunks.'

    def handle(self, *args, **kwargs):
        self.stdout.write(self.style.SUCCESS("Iniciando a criação da biblioteca MOL2 (estratégia híbrida)..."))

        output_dir = os.path.join(settings.BASE_DIR, 'data', 'library_chunks')
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)
        os.makedirs(output_dir, exist_ok=True)
        
        temp_conversion_dir = os.path.join(output_dir, 'temp_conversion')
        os.makedirs(temp_conversion_dir, exist_ok=True)

        moleculas = list(similaridade.objects.filter(publicada=True))
        total_mols = len(moleculas)
        num_chunks = 8
        chunk_size = math.ceil(total_mols / num_chunks)
        
        self.stdout.write(f"Encontradas {total_mols} moléculas. A dividir em {num_chunks} chunks.")

        processadas = 0
        falhas = 0
        falhas_sem_sdf = 0
        falhas_arquivo_inexistente = 0
        falhas_arquivo_vazio = 0
        falhas_conversao = 0
        falhas_erro = 0

        # Preparar dados para processamento paralelo
        mol_data_list = []
        for i in range(num_chunks):
            start_index = i * chunk_size
            end_index = start_index + chunk_size
            moleculas_chunk = moleculas[start_index:end_index]
            for mol_obj in moleculas_chunk:
                mol_data_list.append((mol_obj, temp_conversion_dir, i))
        
        # Determinar número de workers (CPU cores - 1 para não sobrecarregar)
        num_workers = min(multiprocessing.cpu_count() - 1, 8)  # Máximo 8 workers
        self.stdout.write(f"Usando {num_workers} workers para processamento paralelo")
        self.stdout.write(f"Processando {total_mols} moléculas em paralelo...")
        
        # Processar moléculas em paralelo
        chunk_contents = {i: [] for i in range(num_chunks)}
        
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            # Submeter todas as tarefas
            future_to_mol = {executor.submit(process_single_molecule, mol_data): mol_data for mol_data in mol_data_list}
            
            # Processar resultados conforme completam
            for future in as_completed(future_to_mol):
                mol_data = future_to_mol[future]
                try:
                    result = future.result()
                    
                    if result['success']:
                        processadas += 1
                        chunk_contents[result['chunk_index']].append(result['content'])
                        logging.info(f"[PARALLEL] SUCESSO: Molécula {result['mol_id']} processada")
                    else:
                        falhas += 1
                        # Contar tipo de falha
                        if result['error'] == 'sem_sdf':
                            falhas_sem_sdf += 1
                        elif result['error'] == 'arquivo_inexistente':
                            falhas_arquivo_inexistente += 1
                        elif result['error'] == 'arquivo_vazio':
                            falhas_arquivo_vazio += 1
                        elif result['error'] == 'conversao_falhou':
                            falhas_conversao += 1
                        else:
                            falhas_erro += 1
                        
                        logging.error(f"[PARALLEL] ERRO: Molécula {result['mol_id']} - {result['message']}")
                    
                    # Atualizar progresso
                    self.stdout.write(f"Progresso: {processadas + falhas}/{total_mols}", ending='\r')
                    
                except Exception as e:
                    falhas += 1
                    falhas_erro += 1
                    logging.error(f"[PARALLEL] ERRO INESPERADO: {e}")
        
        # Escrever chunks
        for i in range(num_chunks):
            chunk_path = os.path.join(output_dir, f'database_chunk_{i}.mol2')
            logging.info(f"[CHUNK] Escrevendo chunk {i}: {chunk_path}")
            
            with open(chunk_path, 'w', encoding='utf-8') as outfile:
                for content in chunk_contents[i]:
                    outfile.write(content + "\n@<TRIPOS>MOLECULE\n")
            
            chunk_size_final = os.path.getsize(chunk_path) if os.path.exists(chunk_path) else 0
            logging.info(f"[CHUNK] Chunk {i} finalizado com {len(chunk_contents[i])} moléculas, {chunk_size_final} bytes")
        
        logging.info(f"[CLEANUP] Removendo diretório temporário: {temp_conversion_dir}")
        shutil.rmtree(temp_conversion_dir)
        
        # Relatório final
        self.stdout.write(f"\n{self.style.SUCCESS('Processamento concluído!')}")
        self.stdout.write(f"Moléculas processadas com sucesso: {processadas}")
        self.stdout.write(f"Total de falhas: {falhas}")
        self.stdout.write(f"Taxa de sucesso: {(processadas/(processadas+falhas)*100):.1f}%" if (processadas+falhas) > 0 else "Nenhuma molécula processada")
        
        if falhas > 0:
            self.stdout.write(f"\n{self.style.WARNING('Detalhes das falhas:')}")
            self.stdout.write(f"  - Sem arquivo SDF associado: {falhas_sem_sdf}")
            self.stdout.write(f"  - Arquivo SDF não encontrado: {falhas_arquivo_inexistente}")
            self.stdout.write(f"  - Arquivo SDF vazio: {falhas_arquivo_vazio}")
            self.stdout.write(f"  - Falha na conversão (SDF inválido): {falhas_conversao}")
            self.stdout.write(f"  - Erro inesperado: {falhas_erro}")
