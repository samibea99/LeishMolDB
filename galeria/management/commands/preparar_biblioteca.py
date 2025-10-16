import os
import re
import math
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

from django.core.management.base import BaseCommand
from django.conf import settings
from galeria.models import similaridade

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    raise ImportError("A biblioteca RDKit não foi encontrada. Por favor, instale-a com: conda install -c conda-forge rdkit")

def slugify(value):
    """Converte uma string em um nome de arquivo seguro."""
    value = str(value)
    value = re.sub(r'[^\w\s-]', '', value).strip().lower()
    value = re.sub(r'[-\s]+', '_', value)
    return value

def processar_e_converter_molecula(mol_obj):
    """
    Função trabalhadora que converte um arquivo SDF para o formato MOL2
    (VERSÃO FINAL - Construção manual do MOL2)
    """
    resultado = {
        "status": "sucesso", "mol_id": mol_obj.id, "mol_nome": mol_obj.nome,
        "mol2_content": None, "smile_gerado": None, "erro": "Nenhum erro."
    }

    if not mol_obj.sdf or not os.path.exists(mol_obj.sdf.path):
        resultado["status"] = "falha"
        resultado["erro"] = "Arquivo SDF não encontrado."
        return resultado

    caminho_sdf_entrada = mol_obj.sdf.path

    try:
        suppl = Chem.SDMolSupplier(caminho_sdf_entrada, strictParsing=True)
        mol = next(suppl)

        if mol is None:
            resultado["status"] = "falha"
            resultado["erro"] = "RDKit não conseguiu ler uma molécula válida do arquivo SDF."
            return resultado

        mol_com_H = Chem.AddHs(mol, addCoords=True)
        AllChem.EmbedMolecule(mol_com_H, AllChem.ETKDGv3())

        try:
            AllChem.MMFFOptimizeMolecule(mol_com_H)
        except Exception:
            AllChem.UFFOptimizeMolecule(mol_com_H)

        # --- ABORDAGEM FINAL: Construção Manual do Arquivo MOL2 ---
        
        # 1. Gerar um nome temporário seguro para a molécula
        nome_seguro = slugify(mol_obj.nome if mol_obj.nome else f"mol_{mol_obj.id}")
        mol_com_H.SetProp("_Name", nome_seguro)

        # 2. Escrever a molécula em um arquivo temporário usando um "writer"
        # Esta é a forma mais fundamental de escrita no RDKit.
        # Se isto falhar, nada mais funcionará.
        temp_file_path = f"temp_{mol_obj.id}.mol2"
        writer = Chem.MolToMol2File(mol_com_H, temp_file_path)
        
        # 3. Ler o conteúdo do arquivo que acabamos de criar
        with open(temp_file_path, 'r') as f:
            mol2_block = f.read()

        # 4. Apagar o arquivo temporário
        os.remove(temp_file_path)
        
        # 5. Adicionar nosso cabeçalho customizado (como antes)
        header = f"@<TRIPOS>MOLECULE\nmol_{mol_obj.id}_{nome_seguro}\n"
        mol2_lines = mol2_block.splitlines()
        mol_line_index = -1
        for i, line in enumerate(mol2_lines):
            if line.strip() == "@<TRIPOS>MOLECULE":
                mol_line_index = i
                break
        
        if mol_line_index != -1:
            final_mol2_content = header + "\n".join(mol2_lines[mol_line_index+2:]) # Pula duas linhas no MOL2
            resultado["mol2_content"] = final_mol2_content
        else:
             resultado["status"] = "falha"
             resultado["erro"] = "Não foi possível gerar um cabeçalho MOL2 válido."

    except Exception as e:
        resultado["status"] = "falha"
        resultado["erro"] = f"Erro inesperado durante a conversão com RDKit: {str(e)}"

    return resultado

class Command(BaseCommand):
    help = 'Prepara a biblioteca para busca 3D usando RDKit, com opção de limitar o número de moléculas.'

    def add_arguments(self, parser):
        parser.add_argument(
            '--limit', type=int, default=None,
            help='Limita o número de moléculas a serem processadas para fins de teste.'
        )

    def handle(self, *args, **kwargs):
        limit = kwargs['limit']

        self.stdout.write(self.style.SUCCESS("--- Iniciando preparação da biblioteca com RDKit (versão final) ---"))
        
        output_dir = os.path.join(settings.BASE_DIR, 'data', 'library_chunks')
        if os.path.exists(output_dir): shutil.rmtree(output_dir)
        os.makedirs(output_dir)
        
        moleculas = similaridade.objects.filter(publicada=True).exclude(sdf__exact='').exclude(sdf__isnull=True)
        
        if limit:
            moleculas = moleculas[:limit]
            self.stdout.write(self.style.WARNING(f"PROCESSANDO EM MODO DE TESTE: APENAS AS PRIMEIRAS {limit} MOLÉCULAS."))

        total_mols = len(moleculas)
        if total_mols == 0:
            self.stdout.write(self.style.WARNING("Nenhuma molécula encontrada para processar."))
            return
            
        num_chunks = 8
        chunk_size = math.ceil(total_mols / num_chunks) if total_mols > 0 else 1
        self.stdout.write(f"Encontradas {total_mols} moléculas para converter com RDKit.")

        sucessos = 0
        falhas_detalhadas = []
        chunk_contents = {i: [] for i in range(num_chunks)}
        
        with ProcessPoolExecutor(max_workers=multiprocessing.cpu_count()) as executor:
            tasks = {executor.submit(processar_e_converter_molecula, mol): math.floor(i / chunk_size) for i, mol in enumerate(moleculas)}
            for future in as_completed(tasks):
                chunk_index = tasks[future]
                resultado = future.result()
                if resultado['status'] == 'sucesso':
                    sucessos += 1
                    chunk_contents[chunk_index].append(resultado['mol2_content'])
                else:
                    falhas_detalhadas.append({"id": resultado['mol_id'], "nome": resultado['mol_nome'], "erro": resultado['erro']})
                
                progresso = sucessos + len(falhas_detalhadas)
                self.stdout.write(f"Progresso: {progresso}/{total_mols} | Válidas: {sucessos} | Inválidas (puladas): {len(falhas_detalhadas)}", ending='\r')
                self.stdout.flush()

        self.stdout.write("\n")
        self.stdout.write(self.style.SUCCESS("\n--- Processamento Concluído! ---"))
        self.stdout.write(f"Moléculas Válidas (convertidas): {sucessos}")
        self.stdout.write(f"Moléculas Inválidas (puladas): {len(falhas_detalhadas)}")

        if falhas_detalhadas:
            log_path = os.path.join(settings.BASE_DIR, 'preparacao_erros.log')
            self.stdout.write(self.style.WARNING(f"\nOs detalhes sobre as {len(falhas_detalhadas)} moléculas inválidas foram salvos em: {log_path}"))
            with open(log_path, 'w') as f:
                f.write("--- LOG DE ERROS NA PREPARAÇÃO DA BIBLIOTECA (RDKit) ---\n\n")
                for falha in sorted(falhas_detalhadas, key=lambda x: x['id']):
                    f.write(f"ID: {falha['id']}\nNome: {falha['nome']}\nErro: {falha['erro']}\n" + "-" * 20 + "\n")

        self.stdout.write("\nEscrevendo arquivos de chunk da biblioteca...")
        for i in range(num_chunks):
            if not chunk_contents[i]: continue
            chunk_path = os.path.join(output_dir, f'database_chunk_{i}.mol2')
            with open(chunk_path, 'w', encoding='utf-8') as outfile:
                outfile.write("\n".join(chunk_contents[i]))
                outfile.write("\n")
        self.stdout.write("Arquivos de chunk finalizados.")

        self.stdout.write(self.style.SUCCESS("\nPreparação da biblioteca finalizada com sucesso."))
