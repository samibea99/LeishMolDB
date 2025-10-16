from django.shortcuts import render, redirect, get_object_or_404
from .models import similaridade
from django.conf import settings
from .utils import DescritoresUtil
import subprocess
import logging
import os
import glob
import re
import shutil
import pprint 
import sys
import pickle
import billiard as multiprocessing
import uuid
from PIL import Image, ImageDraw, ImageFont
from rdkit.DataStructs import TanimotoSimilarity
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D, Generate
from django.core.files.storage import default_storage
from django.http import JsonResponse, HttpResponse, HttpResponseBadRequest, FileResponse
from django.views.decorators.http import require_GET
from django.views.decorators.cache import cache_page
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Descriptors, Lipinski, rdMolDescriptors, rdmolfiles, rdMolAlign, rdFMCS, FilterCatalog, Draw
from rdkit.Chem.FilterCatalog import FilterCatalogParams
from .forms import UploadSDFForm
from django.db.models import Q 
import pandas as pd 
from django.template.loader import render_to_string 
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter
from django.core.paginator import Paginator 
from datetime import datetime
import tempfile 
import logging   

logger = logging.getLogger(__name__)

def index(request):
    return render(request, 'galeria/index.html')

def index(request):
    total_linhas = similaridade.objects.count()
    return render(request, 'galeria/index.html', {'total_linhas': total_linhas}) 

@require_GET
@cache_page(60 * 5)  # cache por 5 minutos para evitar recomputo frequente
def molecules_graph(request):
    """
    Retorna nós e arestas de similaridade molecular usando Tanimoto.
    Ignora moléculas com SMILES inválidos.
    """
    mols = similaridade.objects.all().only('id', 'nome', 'smile')

    # monta os nós
    nodes = [
        {
            "id": str(m.id),
            "label": m.nome or f"Mol-{m.id}",
            "smiles": m.smile or "",
            "name": m.nome or "",
            "href": m.url or f"/molecule/{m.id}"
        }
        for m in mols
    ]

    fps = {}
    for m in mols:
        s = (m.smile or "").strip()
        if not s:
            print(f"[AVISO] Molécula {m.id} sem SMILES.")
            fps[m.id] = None
            continue
        try:
            mol_obj = Chem.MolFromSmiles(s)
            if mol_obj is None:
                print(f"[AVISO] Molécula {m.id} com SMILES inválido: {s}")
                fps[m.id] = None
                continue
            fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol_obj, radius=2, nBits=2048)
            fps[m.id] = fp
        except Exception as e:
            print(f"[ERRO] Falha ao processar molécula {m.id}: {e}")
            fps[m.id] = None

    try:
        threshold = float(request.GET.get("threshold", 0.2))
    except ValueError:
        threshold = 0.2

    ids = [m.id for m in mols]
    edges = []
    for i in range(len(ids)):
        for j in range(i + 1, len(ids)):
            fi, fj = fps.get(ids[i]), fps.get(ids[j])
            if fi is None or fj is None:
                continue
            try:
                sim = DataStructs.TanimotoSimilarity(fi, fj)
                if sim >= threshold:
                    edges.append({
                        "source": str(ids[i]),
                        "target": str(ids[j]),
                        "weight": float(sim)
                    })
            except Exception as e:
                print(f"[ERRO] Similaridade {ids[i]}-{ids[j]}: {e}")

    return JsonResponse({"nodes": nodes, "edges": edges})

@require_GET
def molblock_from_smiles(request):
    """
    Recebe um SMILES e retorna um MolBlock 3D (molfile) para visualização no 3Dmol.js.
    """
    smiles = request.GET.get("smiles")
    if not smiles:
        return HttpResponseBadRequest("Missing SMILES parameter")
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return HttpResponseBadRequest("Invalid SMILES")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        mol = Chem.RemoveHs(mol)
        molblock = Chem.MolToMolBlock(mol)
        return HttpResponse(molblock, content_type="chemical/x-mdl-molfile")
    except Exception as e:
        return HttpResponseBadRequest(str(e))

def moleculas(request):
    # Obter todas as moléculas publicadas
    mol_list = similaridade.objects.filter(publicada=True)

    # Paginação para exibir 96 cards por página
    paginator = Paginator(mol_list, 96)  

    # Obter o número da página a partir da requisição GET
    page_number = request.GET.get('page')
    
    # Obter a página com as moléculas
    mol = paginator.get_page(page_number)

    # Renderizar o template 'moleculas.html' com todas as moléculas misturadas
    return render(request, 'galeria/moleculas.html', {"cards": mol})

def moleculas_list_view(request):
    molecula_list = similaridade.objects.all()

    # Extração das propriedades é realizada antes da paginação, mas em um subconjunto se necessário
    for molecula in molecula_list:
        # Supondo que a função extract_molecule_properties retorna os dados necessários
        # molecula.dados = extract_molecule_properties(molecula.sdf.path)[0] if molecula.sdf.path else {}
        pass

    paginator = Paginator(molecula_list, 54)  # 54 itens por página

    page = request.GET.get('page')
    moleculas = paginator.get_page(page)

    print('Número de moléculas:', molecula_list.count())
    print('Número da página atual:', moleculas.number)
    print('Total de páginas:', moleculas.paginator.num_pages)

    return render(request, 'galeria/moleculas.html', {'moleculas': moleculas})

def ordenar_moleculas_view(request):
    # Receber parâmetros de categoria e ordenação
    ordenacao = request.GET.get('ordenacao')
    categoria = request.GET.get('categoria')

    # Log para verificar os parâmetros recebidos
    print(f"Parâmetros recebidos - Categoria: {categoria}, Ordenação: {ordenacao}")

    # Obter todas as moléculas
    moleculas = similaridade.objects.all()

    # Aplicar filtro de categoria, se especificado
    if categoria:
        moleculas = moleculas.filter(categoria=categoria)
        # Log para verificar as moléculas após filtragem
        print(f"Moléculas após filtragem por categoria '{categoria}': {moleculas}")

    # Extração de dados de logP e peso molecular
    moleculas_com_dados = []
    for mol in moleculas:
        dados = extract_molecular_weight_and_logp(mol.sdf.path)
        if dados:
            moleculas_com_dados.append((mol, dados))
            # Log para verificar dados extraídos
            print(f"Molécula ID: {mol.id}, Peso molecular: {dados['peso_molecular']}, LogP: {dados['logp']}")

    # Ordenar a lista de moléculas com base no critério de ordenação especificado
    if ordenacao:
        if ordenacao == 'peso_molecular':
            moleculas_com_dados.sort(key=lambda x: x[1]['peso_molecular'])
            print(f"Moléculas ordenadas por peso molecular: {moleculas_com_dados}")
        elif ordenacao == 'logp':
            moleculas_com_dados.sort(key=lambda x: x[1]['logp'])
            print(f"Moléculas ordenadas por logP: {moleculas_com_dados}")

    # Extração das moléculas ordenadas
    moleculas_ordenadas = [mol_dados[0] for mol_dados in moleculas_com_dados]

    # Renderizar o template com as moléculas ordenadas
    return render(request, 'galeria/moleculas.html', {'cards': moleculas_ordenadas})

def extract_data_from_sdf(sdf_path):
    # Inicializa um dicionário para armazenar os dados extraídos.
    data = {
        'Molecules': [],
        'MolecularWeights': [],
        'LogPs': [],
        'HDonors': [],
        'HAcceptors': [],
        'LipinskiRulesMet': [],
        'tpsa': []  # Adiciona uma chave para TPSA.
    }

    # Cria um objeto SDMolSupplier do RDKit que pode ler moléculas de um arquivo SDF.
    suppl = Chem.SDMolSupplier(str(sdf_path))

    # Itera sobre todas as moléculas no arquivo SDF.
    for idx, mol in enumerate(suppl):
        if mol:  # Verifica se a molécula atual é válida (não nula).
            mol_with_hydrogens = Chem.AddHs(mol)
            smiles = Chem.MolToSmiles(mol_with_hydrogens)
            mol_weight = Descriptors.MolWt(mol_with_hydrogens)
            log_p = Descriptors.MolLogP(mol_with_hydrogens)
            h_donors = Lipinski.NumHDonors(mol_with_hydrogens)
            h_acceptors = Lipinski.NumHAcceptors(mol_with_hydrogens)
            tpsa = rdMolDescriptors.CalcTPSA(mol_with_hydrogens)

            # Avalia as regras de Lipinski.
            rules_met = 0
            if mol_weight <= 500: rules_met += 1
            if log_p <= 5: rules_met += 1
            if h_donors <= 5: rules_met += 1
            if h_acceptors <= 10: rules_met += 1

            # Adiciona os dados extraídos ao dicionário 'data'.
            data['Molecules'].append(smiles)
            data['MolecularWeights'].append(mol_weight)
            data['LogPs'].append(log_p)
            data['HDonors'].append(h_donors)
            data['HAcceptors'].append(h_acceptors)
            data['LipinskiRulesMet'].append(rules_met)
            data['tpsa'].append(tpsa)  # Adiciona o valor de TPSA ao dicionário.

    return data

def get_mol_data(request, id):
    mol = get_object_or_404(similaridade, id=id)
    sdf_data = mol.sdf.read()  # Lê o conteúdo do arquivo SDF
    return HttpResponse(sdf_data, content_type="chemical/x-mdl-sdfile")

def imagem(request, mol_id):
    mol = get_object_or_404(similaridade, pk=mol_id)
    return render(request, 'galeria/imagem.html', {"mol": mol})

def download_data_view(request, pk):
    objeto = get_object_or_404(similaridade, pk=pk)
    
    # Obtem os dados do objeto e formata como desejar
    dados_para_download = objeto.dados_de_interesse()

    # Agora, cria um arquivo temporário para enviar como resposta
    filename = "dados_para_download.txt"
    with open(filename, 'w') as f:
        f.write(dados_para_download)

    # Retornar o arquivo como resposta para download
    response = FileResponse(open(filename, 'rb'))
    response['Content-Disposition'] = 'attachment; filename="dados_para_download.txt"'
    return response

def buscar(request):
    # Receber todos os parâmetros da requisição GET
    categoria = request.GET.get('categoria')
    nome_a_buscar = request.GET.get('buscar')
    ordenacao = request.GET.get('ordenacao')
    search_type = request.GET.get('search_type', 'text') # Novo parâmetro, padrão 'text'

    # Log dos parâmetros recebidos
    logging.info(f"Busca recebida: termo='{nome_a_buscar}', tipo='{search_type}', categoria='{categoria}', ordenacao='{ordenacao}'")

    # Inicializa a consulta com todas as moléculas publicadas
    mol_query = similaridade.objects.filter(publicada=True)
    
    # Aplica o filtro de categoria, se houver
    if categoria:
        mol_query = mol_query.filter(categoria__icontains=categoria)

    # Aplica o filtro de busca, dependendo do tipo
    if nome_a_buscar:
        if search_type == 'substructure':
            logging.info(f"Iniciando busca por subestrutura com SMARTS/SMILES: {nome_a_buscar}")
            try:
                query_mol = Chem.MolFromSmarts(nome_a_buscar)
                if query_mol is None:
                    query_mol = Chem.MolFromSmiles(nome_a_buscar)
                
                if query_mol is None:
                    raise ValueError("SMILES/SMARTS da subestrutura é inválido.")

                matching_ids = []
                # NOTA DE EFICIÊNCIA: Este laço itera sobre todo o banco de dados em Python.
                # Para bancos de dados muito grandes (>10.000 moléculas), isso pode se tornar lento.
                # Uma otimização futura seria usar um cartucho de banco de dados como o RDKit para PostgreSQL.
                for mol_db in similaridade.objects.all().only('id', 'smile'):
                    if not mol_db.smile: continue
                    db_mol = Chem.MolFromSmiles(mol_db.smile)
                    if db_mol and db_mol.HasSubstructMatch(query_mol):
                        matching_ids.append(mol_db.id)
                
                mol_query = mol_query.filter(id__in=matching_ids)
                logging.info(f"Encontradas {len(matching_ids)} moléculas com a subestrutura.")

            except Exception as e:
                logging.error(f"Erro na busca por subestrutura: {e}")
                mol_query = mol_query.none() # Retorna uma lista vazia em caso de erro

        else: # Busca por texto padrão (lógica que você já tinha)
            try:
                nome_a_buscar_as_int = int(nome_a_buscar)
                mol_query = mol_query.filter(id=nome_a_buscar_as_int)
            except ValueError:
                mol_query = mol_query.filter(
                    Q(nome__icontains=nome_a_buscar) | 
                    Q(smile__icontains=nome_a_buscar) | 
                    Q(categoria__icontains=nome_a_buscar) | 
                    Q(url__icontains=nome_a_buscar)
                )

    # Lógica de Ordenação (integrada com a sua função original)
    if ordenacao in ['peso_molecular', 'logp']:
        moleculas_com_dados = []
        for mol in mol_query:
            # NOTA DE EFICIÊNCIA: Esta função lê um arquivo do disco para cada molécula, em cada requisição.
            dados = extract_molecular_weight_and_logp(mol.sdf.path)
            if dados:
                moleculas_com_dados.append({
                    'mol': mol,
                    'peso_molecular': dados.get('peso_molecular', 0),
                    'logp': dados.get('logp', 0)
                })
        
        # Ordena a lista de dicionários
        reverse_sort = False # Ordenação ascendente
        moleculas_com_dados.sort(key=lambda x: x.get(ordenacao, 0), reverse=reverse_sort)
        
        # Extrai apenas os objetos 'mol' ordenados
        moleculas_ordenadas = [item['mol'] for item in moleculas_com_dados]
    else:
        # Ordenação padrão por ID (ou nome, se preferir)
        moleculas_ordenadas = mol_query.order_by('id')

    # Paginação
    paginator = Paginator(moleculas_ordenadas, 96)
    page_number = request.GET.get('page')
    cards = paginator.get_page(page_number)

    # Renderiza o template com todo o contexto necessário
    return render(request, 'galeria/moleculas.html', {
        "cards": cards, 
        "categoria": categoria, 
        "ordenacao": ordenacao,
        "buscar": nome_a_buscar,
        "search_type": search_type
    })

def extract_molecular_weight_and_logp(sdf_path):
    
    # Cria um fornecedor de moléculas a partir do arquivo SDF, que permite iterar sobre as moléculas nele contidas.
    suppl = Chem.SDMolSupplier(str(sdf_path))
    
    # Obtém a primeira molécula não nula do fornecedor. Isso é feito utilizando uma expressão geradora
    # com um 'if' condicional para filtrar moléculas nulas (inválidas).
    mol = next((m for m in suppl if m is not None), None)  # Apenas a primeira molécula válida.

    # Verifica se uma molécula válida foi encontrada; se não, retorna None.
    if mol is None:
        return None

    # Adiciona átomos de hidrogênio explicitamente à molécula, 
    # para cálculos precisos de propriedades químicas.
    mol_with_hydrogens = Chem.AddHs(mol)
    
    # Calcula o peso molecular da molécula com hidrogênios adicionados.
    mol_weight = Descriptors.MolWt(mol_with_hydrogens)
    
    # Calcula o log P (coeficiente de partição octanol-água) da molécula.
    log_p = Descriptors.MolLogP(mol_with_hydrogens)

    # Retorna um dicionário contendo o peso molecular e log P da molécula.
    return {'peso_molecular': mol_weight, 'logp': log_p} 

def molecula_view(request, id):
    # Obtem a instância do modelo com base no ID
    similaridade_instancia = get_object_or_404(similaridade, id=id)
    
    # O caminho do arquivo SDF é obtido pelo atributo 'path' do FileField
    sdf_path = similaridade_instancia.sdf.path

    # Chama a função para extrair os dados do arquivo SDF
    data = extract_data_from_sdf(sdf_path)

    # Zipa as listas: pesos moleculares, LogPs, doadores de H e aceitadores de H
    moleculas_data = [
        {
            'peso_molecular': mw,
            'logp': logp,
            'h_donors': hd,
            'h_acceptors': ha,
            'lipinski_rules_met': lr,
            'tpsa': tpsa
        }
        for mw, logp, hd, ha, lr, tpsa in zip(data['MolecularWeights'], data['LogPs'], data['HDonors'], data['HAcceptors'], data['LipinskiRulesMet'], data['tpsa'])
    ]

    # Passa os dados extraídos para o contexto do template
    return render(request, 'galeria/imagem.html', {'moleculas_data': moleculas_data, 'mol': similaridade_instancia})

def download_sdfs(request):
    # Cria um arquivo temporário
    temp_file = tempfile.TemporaryFile(mode='w+')

    # Itera sobre os objetos do modelo e escreve seu conteúdo no arquivo temporário
    for obj in similaridade.objects.all():
        sdf_path = obj.sdf.path
        with open(sdf_path, 'r') as sdf:
            temp_file.write(sdf.read() + '\n$$$$\n')

    # Move o ponteiro do arquivo para o início
    temp_file.seek(0)

    # Cria uma resposta HTTP com o conteúdo do arquivo temporário
    response = HttpResponse(temp_file.read(), content_type='chemical/x-mdl-sdfile')
    response['Content-Disposition'] = 'attachment; filename="combined_files.sdf"'

    # Fecha o arquivo temporário
    temp_file.close()

    return response

def comparar(request):
    return render(request, 'galeria/comparar.html') 

def resultado_similaridade(request):
    if request.method == 'POST':
        form = UploadSDFForm(request.POST, request.FILES)
        if form.is_valid():
            uploaded_file = request.FILES['file']
            if not uploaded_file.name.endswith('.sdf'):
                error_message = "O arquivo enviado não é um SDF."
                return render(request, 'galeria/comparar.html', {'form': form, 'error_message': error_message})

            uploaded_mol_block = uploaded_file.read()

            # Verifica se o conteúdo é binário e decodifica se necessário
            if isinstance(uploaded_mol_block, bytes):
                uploaded_mol = Chem.MolFromMolBlock(uploaded_mol_block.decode('utf-8'))
            else:
                uploaded_mol = Chem.MolFromMolBlock(uploaded_mol_block)

            # Continua o processamento apenas se uploaded_mol foi criado com sucesso
            if uploaded_mol:
                similaridades = similaridade.objects.all()
                similarities_with_ids = []

                for similaridade_obj in similaridades:
                    sdf_path = similaridade_obj.sdf.path
                    mol_supplier = Chem.SDMolSupplier(sdf_path)

                    for ref_mol in mol_supplier:
                        if ref_mol:
                            # Cálculo da similaridade
                            fps_ref = AllChem.GetMorganFingerprint(ref_mol, 2)
                            fps_upload = AllChem.GetMorganFingerprint(uploaded_mol, 2)
                            similarity = DataStructs.TanimotoSimilarity(fps_ref, fps_upload)
                            
                            # Converte a similaridade para porcentagem e arredonda para 2 casas decimais
                            similarity_percentage = f"{similarity * 100:.2f}%"
                            # Adiciona a similaridade em porcentagem e o objeto similaridade à lista
                            similarities_with_ids.append((similarity, similarity_percentage, similaridade_obj))

                # Ordena a lista pelo valor de similaridade, em ordem decrescente
                similarities_with_ids.sort(key=lambda x: x[0], reverse=True)

                # Seleciona as 10 maiores similaridades
                top_10_similarities = similarities_with_ids[:10]

                return render(request, 'galeria/resultado_similaridade.html', {'similarities_with_ids': top_10_similarities})

def get_molecule_data(request):
    molecule_id = request.GET.get('id')
    similaridade_obj = similaridade.objects.get(id=molecule_id)  # Obtém o objeto específico
    sdf_path = similaridade_obj.sdf.path  # Obtém o caminho do arquivo no sistema de arquivos

    # O RDKit espera um caminho de arquivo ou um objeto de arquivo 
    suppl = Chem.SDMolSupplier(str(sdf_path))
    mols = [mol for mol in suppl if mol is not None]

    # Assumindo que você quer trabalhar com a primeira molécula
    mol = mols[0] if mols else None

    if mol:
        # Gera coordenadas 3D se não existirem 
        if not mol.GetNumConformers():
            AllChem.Compute2DCoords(mol)
            AllChem.EmbedMolecule(mol)
            AllChem.UFFOptimizeMolecule(mol)   

    # Converte para MolBlock
    molblock = Chem.MolToMolBlock(mol)

    return JsonResponse({"molblock": molblock})

def comparar3d(request):
    form = UploadSDFForm()
    return render(request, 'galeria/comparar3d.html', {'form': form})

def has_3d_coordinates(sdf_path):
    mols = Chem.SDMolSupplier(sdf_path, removeHs=False)
    for mol in mols:
        if mol is None:
            continue
        conf = mol.GetConformer()
        positions = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
        if all(pos.z != 0.0 for pos in positions):
            return True
    return False

def extrair_metricas_completas(log_file_path):
    """
    Lê o arquivo de log do SENSAAS e extrai as métricas finais da
    seção "Selected solution".
    """
    metrics = {'gfit': 0.0, 'hfit': 0.0, 'rmse': 0.0}
    try:
        with open(log_file_path, 'r', encoding='utf-8') as f:
            content = f.read()

        # Procura pela seção da solução final no log
        selected_solution_section = re.search(r"Selected solution \(best gfit \+ hfit\):(.+)", content, re.DOTALL)

        if selected_solution_section:
            section_text = selected_solution_section.group(1)
            
            # Extrai gfit e hfit da última linha (ex: gfit= 0.453 cfit= 0.328 hfit= 0.247 ...)
            final_line_match = re.search(r"gfit=\s*([0-9.]+).*hfit=\s*([0-9.]+)", section_text)
            if final_line_match:
                metrics['gfit'] = float(final_line_match.group(1))
                metrics['hfit'] = float(final_line_match.group(2))

            # Extrai o rmse da linha que contém "rmse_all_dots"
            rmse_match = re.search(r"rmse_all_dots=\s*([0-9.]+)", section_text)
            if rmse_match:
                metrics['rmse'] = float(rmse_match.group(1))
        
        print(f"[INFO] Métricas extraídas: {metrics}")

    except Exception as e:
        print(f"[ERRO] Falha ao extrair métricas do SENSAAS: {e}")
    
    return metrics

logger = logging.getLogger(__name__)

logging.basicConfig(level=logging.INFO, format='%(asctime)s - [%(levelname)s] - %(message)s')

def robust_convert_sdf_to_mol2(sdf_path, mol2_path):
    """
    Converte SDF para MOL2 de forma robusta, com fallback e verificação de saída.
    Isto corrige o erro 'Could not load ligands' do LS-align.
    """
    logging.info(f"Iniciando conversão SDF->MOL2 para: {os.path.basename(sdf_path)}")
    
    # Tentativa 1: Com otimização
    try:
        comando = ["obabel", "-i", "sdf", sdf_path, "-o", "mol2", "--gen3d", "--minimize", "--ff", "MMFF94", "-e"]
        proc = subprocess.run(comando, capture_output=True, text=True, check=True, timeout=180)
        with open(mol2_path, 'w', encoding='utf-8') as f:
            f.write(proc.stdout)
        if os.path.exists(mol2_path) and os.path.getsize(mol2_path) > 0:
            logging.info("Conversão SDF->MOL2 com otimização bem-sucedida.")
            return True, "Sucesso com otimização."
        raise ValueError(f"Comando executou mas não gerou MOL2. Stderr: {proc.stderr.strip()}")
    except Exception as e_opt:
        logging.warning(f"Método de otimização falhou: {getattr(e_opt, 'stderr', str(e_opt))}. Tentando método simples...")

    # Tentativa 2: Método simples (fallback)
    try:
        comando_simples = ["obabel", "-i", "sdf", sdf_path, "-o", "mol2", "--gen3d", "-e"]
        proc_simple = subprocess.run(comando_simples, capture_output=True, text=True, check=True, timeout=60)
        with open(mol2_path, 'w', encoding='utf-8') as f:
            f.write(proc_simple.stdout)
        if os.path.exists(mol2_path) and os.path.getsize(mol2_path) > 0:
            logging.info("Conversão SDF->MOL2 com método simples bem-sucedida.")
            return True, "Sucesso com método simples."
        raise ValueError(f"Comando simples executou mas não gerou MOL2. Stderr: {proc_simple.stderr.strip()}")
    except Exception as e_simple:
        error_message = f"Todos os métodos de conversão SDF->MOL2 falharam. Erro final: {getattr(e_simple, 'stderr', str(e_simple))}"
        logging.error(error_message)
        return False, error_message


def convert_mol2_to_sdf_with_obabel(mol2_path, sdf_path):
    """
    Converte MOL2 para SDF de forma robusta, com um fallback inteligente
    para contornar incompatibilidades de formato.
    """
    logging.info(f"Iniciando conversão MOL2->SDF para: {os.path.basename(mol2_path)}")
    
    # --- TENTATIVA 1: MÉTODO DIRETO ---
    try:
        comando_direto = ["obabel", "-i", "mol2", mol2_path, "-o", "sdf"]
        proc = subprocess.run(comando_direto, capture_output=True, text=True, check=True, timeout=180)
        
        with open(sdf_path, 'w', encoding='utf-8') as f:
            f.write(proc.stdout)
            
        if not os.path.exists(sdf_path) or os.path.getsize(sdf_path) == 0:
            error_msg = f"Falha na conversão direta. Detalhe: {proc.stderr.strip() if proc.stderr else 'Saída vazia.'}"
            raise IOError(error_msg) # Força a ida para o bloco 'except'

        logging.info("Conversão MOL2->SDF (método direto) bem-sucedida.")
        return True, "Sucesso"

    except (subprocess.CalledProcessError, IOError, Exception) as e:
        logging.warning(f"Método direto de conversão falhou: {e}. Acionando workaround de reconstrução 2D->3D.")

        # --- TENTATIVA 2: WORKAROUND DE RECONSTRUÇÃO ---
        try:
            logging.info("Workaround: Gerando SDF 2D intermediário...")
            sdf_2d_path = sdf_path.replace('.sdf', '_2d.sdf')
            comando_2d = ["obabel", "-i", "mol2", mol2_path, "-o", "sdf", "--gen2d"]
            subprocess.run(comando_2d, capture_output=True, text=True, check=True, timeout=180)
            
            logging.info("Workaround: Gerando coordenadas 3D a partir do SDF 2D...")
            comando_3d = ["obabel", "-i", "sdf", sdf_2d_path, "-o", "sdf", "--gen3d"]
            proc_3d = subprocess.run(comando_3d, capture_output=True, text=True, check=True, timeout=180)

            with open(sdf_path, 'w', encoding='utf-8') as f:
                f.write(proc_3d.stdout)
            
            if os.path.exists(sdf_path) and os.path.getsize(sdf_path) > 0:
                logging.info("Workaround de reconstrução 2D->3D bem-sucedido!")
                # Limpa o arquivo intermediário
                if os.path.exists(sdf_2d_path):
                    os.remove(sdf_2d_path)
                return True, "Sucesso com workaround."
            else:
                raise ValueError("O workaround de reconstrução também resultou em um arquivo vazio.")

        except Exception as e_workaround:
            final_error_msg = f"TODOS os métodos de conversão falharam. Erro final do workaround: {getattr(e_workaround, 'stderr', str(e_workaround))}"
            logging.error(final_error_msg)
            return False, final_error_msg

def process_user_uploaded_file(uploaded_file, temp_dir):
    """Processa o arquivo enviado pelo usuário, salvando e convertendo para os formatos necessários."""
    logging.info(f"Processando arquivo: {uploaded_file.name}")
    filename = uploaded_file.name.lower()
    
    # Define os caminhos de saída
    user_sdf_path = os.path.join(temp_dir, 'user_query.sdf')
    user_mol2_path = os.path.join(temp_dir, 'user_query.mol2')

    if filename.endswith('.sdf'):
        logging.info("Arquivo .sdf detectado. Salvando e convertendo para .mol2...")
        with open(user_sdf_path, 'wb+') as dest:
            for chunk in uploaded_file.chunks():
                dest.write(chunk)
        success, msg = robust_convert_sdf_to_mol2(user_sdf_path, user_mol2_path)
        if not success:
            raise ValueError(f"Falha ao converter o arquivo SDF enviado: {msg}")

    elif filename.endswith('.mol2'):
        logging.info("Arquivo .mol2 detectado. Salvando e convertendo para .sdf...")
        with open(user_mol2_path, 'wb+') as dest:
            for chunk in uploaded_file.chunks():
                dest.write(chunk)
        success, msg = convert_mol2_to_sdf_with_obabel(user_mol2_path, user_sdf_path)
        if not success:
            # Este erro agora será lançado se a conversão gerar um arquivo vazio
            raise ValueError(f"Falha ao criar arquivo SDF a partir do MOL2 enviado: {msg}")
    else:
        raise ValueError("Formato de arquivo inválido. Por favor, envie um arquivo .sdf ou .mol2.")
        
    return (user_mol2_path, user_sdf_path)

# ... (As funções run_lsalign_worker e parse_lsalign_output permanecem as mesmas) ...
def run_lsalign_worker(task_info):
    """Executa uma instância do LS-align para um chunk da biblioteca."""
    lsalign_executable, user_mol2_path, chunk_path = task_info
    comando = [lsalign_executable, user_mol2_path, chunk_path]
    try:
        # logging.info(f"Executando LS-align para o chunk: {os.path.basename(chunk_path)}")
        return subprocess.run(comando, capture_output=True, text=True, timeout=300, check=True).stdout
    except subprocess.CalledProcessError as e:
        logging.error(f"Erro ao executar LS-align no chunk {os.path.basename(chunk_path)}. Stderr: {e.stderr.strip()}")
        return ""
    except Exception as e:
        logging.error(f"Exceção inesperada no worker do LS-align para {os.path.basename(chunk_path)}: {str(e)}")
        return ""

# ==============================================================================
# FUNÇÃO parse_lsalign_output FINAL (COM LIPINSKI E PAINS)
# ==============================================================================

def parse_lsalign_output(output_text):
    """
    Analisa a saída do LS-align e adiciona alertas de Lipinski e PAINS.
    """
    results = []
    pattern = re.compile(r"^(?!Query_Name)\S+\s+mol_(\d+)\S*\s+([0-9.]+)\s+([0-9.]+)", re.MULTILINE)
    
    for match in pattern.finditer(output_text):
        mol_id_str = match.group(1)
        try:
            mol_id = int(mol_id_str)
            rmsd = float(match.group(2))
            tm_score = float(match.group(3))
            
            mol_obj = similaridade.objects.get(id=mol_id)
            
            violations = 0
            is_pains = False
            try:
                if mol_obj.smile:
                    mol_rdkit = Chem.MolFromSmiles(mol_obj.smile)
                    if mol_rdkit:
                        violations = calculate_lipinski_violations(mol_rdkit)
                        is_pains = has_pains_alert(mol_rdkit)
                    else:
                        logging.warning(f"SMILES inválido no banco para mol_id {mol_id}. Não foi possível calcular alertas.")
            except Exception as e_alert:
                logging.error(f"Erro inesperado ao calcular alertas para mol_id {mol_id}: {e_alert}")
            
            results.append({
                'id': mol_obj.id, 'nome': mol_obj.nome, 'rmsd': rmsd,
                'tm_score': tm_score, 'url': mol_obj.url,
                'lipinski_violations': violations, 'is_pains': is_pains
            })
        except (ValueError, similaridade.DoesNotExist) as e:
            logging.warning(f"Ignorando resultado para mol_id '{mol_id_str}' devido a erro: {e}")
            continue
            
    return results

def perform_rdkit_alignment(user_sdf_path, db_sdf_path):
    """
    Realiza o alinhamento e calcula MCS, Shape Tanimoto e Similaridade Farmacofórica.
    """
    logging.info(f"Executando análise 3D completa entre {os.path.basename(user_sdf_path)} e {os.path.basename(db_sdf_path)}")
    try:
        suppl_user = Chem.SDMolSupplier(user_sdf_path, removeHs=False)
        suppl_db = Chem.SDMolSupplier(db_sdf_path, removeHs=False)
        mol_user = suppl_user[0]
        mol_db = suppl_db[0]
        if mol_user is None or mol_db is None: 
            raise ValueError("Molécula nula encontrada durante o carregamento do SDF.")

        # --- Alinhamento ---
        o3a = AllChem.GetO3A(mol_user, mol_db)
        rmsd = o3a.Align()

        # --- Cálculo do Shape Tanimoto ---
        shape_dist = AllChem.ShapeTanimotoDist(mol_user, mol_db)
        shape_similarity = 1 - shape_dist
        
        # =================== CORREÇÃO DO BUG AQUI ===================
        #
        # A variável 'pharm_factory' agora é definida DENTRO da função,
        # garantindo que ela sempre exista.
        #
        pharm_factory = Gobbi_Pharm2D.factory
        fp_user = Generate.Gen2DFingerprint(mol_user, pharm_factory)
        fp_db = Generate.Gen2DFingerprint(mol_db, pharm_factory)
        pharm_similarity = DataStructs.TanimotoSimilarity(fp_user, fp_db)
        # ==========================================================

        logging.info(f"Alinhamento: RMSD={rmsd:.4f}, ShapeTanimoto={shape_similarity:.4f}, PharmSimilarity={pharm_similarity:.4f}")

        # --- Cálculo do MCS ---
        mcs_result = rdFMCS.FindMCS([mol_user, mol_db], timeout=5)
        user_indices = ()
        db_indices = ()
        if mcs_result.numAtoms > 0:
            mcs_mol = Chem.MolFromSmarts(mcs_result.smartsString)
            if mcs_mol:
                user_indices = mol_user.GetSubstructMatch(mcs_mol)
                db_indices = mol_db.GetSubstructMatch(mcs_mol)

        # --- Geração dos blocos SDF ---
        aligned_user_sdf_block = Chem.MolToMolBlock(mol_user)
        with open(db_sdf_path, 'r', encoding='utf-8') as f:
            db_sdf_block = f.read()
        
        cleaned_db_block = db_sdf_block.split('$$$$')[0].strip()
        cleaned_user_block = aligned_user_sdf_block.strip()
        
        return {
            "combined_sdf": f"{cleaned_db_block}\n$$$$\n{cleaned_user_block}\n$$$$",
            "user_mcs_indices": list(user_indices),
            "db_mcs_indices": list(db_indices),
            "shape_tanimoto": shape_similarity,
            "pharm_similarity": pharm_similarity
        }

    except Exception as e:
        logging.exception("ERRO CRÍTICO dentro da função perform_rdkit_alignment.")
        raise e
# ==============================================================================
# VIEW PRINCIPAL (resultado_3d) - VERSÃO COMPLETA
# ==============================================================================

def resultado_3d(request):
    logging.info("--- NOVA REQUISIÇÃO PARA resultado_3d (Busca Conformacional) ---")
    form = UploadSDFForm(request.POST or None, request.FILES or None)
    context = {"form": form}

    if request.method != "POST":
        return render(request, "galeria/resultado_3d.html", context)
    if not form.is_valid():
        context["error"] = "Formulário inválido."
        return render(request, "galeria/resultado_3d.html", context)

    uploaded_file = request.FILES.get("file")
    if not uploaded_file:
        context["error"] = "Nenhum arquivo foi enviado."
        return render(request, "galeria/resultado_3d.html", context)
    
    session_id = str(uuid.uuid4())
    session_temp_dir = os.path.join(settings.MEDIA_ROOT, "temp_sessions", session_id)
    os.makedirs(session_temp_dir, exist_ok=True)

    try:
        user_original_sdf_path = os.path.join(session_temp_dir, "user_original.sdf")
        with open(user_original_sdf_path, 'wb+') as dest:
            for chunk in uploaded_file.chunks():
                dest.write(chunk)

        multiconformer_sdf_path = os.path.join(session_temp_dir, "user_conformers.sdf")
        generate_conformers(user_original_sdf_path, multiconformer_sdf_path)
        individual_conf_sdfs = split_multiconformer_sdf(multiconformer_sdf_path, session_temp_dir)

        all_raw_results = []
        for conf_sdf_path in individual_conf_sdfs:
            conf_mol2_path = conf_sdf_path.replace('.sdf', '.mol2')
            success, msg = robust_convert_sdf_to_mol2(conf_sdf_path, conf_mol2_path)
            if not success: continue

            chunks_dir = os.path.join(settings.BASE_DIR, "data", "library_chunks")
            lsalign_executable = os.path.join(settings.BASE_DIR, "tools", "lsalign")
            chunk_files = [os.path.join(chunks_dir, f) for f in os.listdir(chunks_dir) if f.endswith(".mol2")]
            tasks = [(lsalign_executable, conf_mol2_path, chunk) for chunk in chunk_files]
            
            with multiprocessing.Pool(processes=len(tasks)) as pool:
                for result_text in pool.imap_unordered(run_lsalign_worker, tasks):
                    parsed_results = parse_lsalign_output(result_text)
                    for res in parsed_results:
                        res['user_conformer_path'] = conf_sdf_path
                    all_raw_results.extend(parsed_results)

        if not all_raw_results:
            raise ValueError("Nenhum resultado de alinhamento válido foi encontrado.")
        
        all_results = consolidate_lsalign_results(all_raw_results)
        best_match = all_results[0]
        best_mol_obj = similaridade.objects.get(id=best_match["id"])
        
        db_sdf_path = best_mol_obj.sdf.path
        with open(user_original_sdf_path, "r", encoding='utf-8') as f:
            conteudo_original = f.read()
        with open(db_sdf_path, "r", encoding='utf-8') as f:
            conteudo_db = f.read()

        best_user_conformer_path = best_match['user_conformer_path']
        alignment_data = perform_rdkit_alignment(best_user_conformer_path, db_sdf_path)

        # --- PREPARAÇÃO DOS DADOS PARA O AVISO E DOWNLOAD ---
        # Compara os nomes dos arquivos para ver se a conformação é diferente da original
        show_conformer_alert = os.path.basename(user_original_sdf_path) != os.path.basename(best_user_conformer_path)
        
        # Gera o caminho relativo para o link de download
        best_conformer_relative_path = os.path.relpath(best_user_conformer_path, session_temp_dir)

        context.update({
            "results": all_results[:50], "best_mol_obj": best_mol_obj, "best_match": best_match,
            "conteudo_molecula_original": conteudo_original, 
            "conteudo_molecula_db": conteudo_db,
            "conteudo_molecula_alinhada": alignment_data["combined_sdf"],
            "user_mcs_indices": alignment_data["user_mcs_indices"],
            "db_mcs_indices": alignment_data["db_mcs_indices"],
            "shape_tanimoto": alignment_data.get("shape_tanimoto"),
            "pharm_similarity": alignment_data.get("pharm_similarity"),
            "user_sdf_path_for_api": user_original_sdf_path,
            "show_conformer_alert": show_conformer_alert, # <- Novo
            "best_conformer_relative_path": best_conformer_relative_path, # <- Novo
            "session_id_for_download": session_id # <- Novo
        })
        return render(request, "galeria/resultado_3d.html", context)

    except Exception as e:
        logging.exception("ERRO INESPERADO NA VIEW resultado_3d (Busca Conformacional)")
        context["error"] = f"Ocorreu um erro inesperado: {e}"
        return render(request, "galeria/resultado_3d.html", context)


# ==============================================================================
# VIEW DA API (get_alignment_data) - VERSÃO COMPLETA E CORRIGIDA
# ==============================================================================

def get_alignment_data(request):
    logging.info("--- NOVA REQUISIÇÃO PARA get_alignment_data (API COM MCS) ---")
    if request.method != 'GET':
        return JsonResponse({'error': 'Método inválido.'}, status=405)

    user_sdf_path = request.GET.get('user_sdf_path')
    db_mol_id = request.GET.get('db_mol_id')
    
    if not all([user_sdf_path, db_mol_id]) or not os.path.exists(user_sdf_path):
        return JsonResponse({'error': 'Arquivo do usuário não encontrado ou ID ausente.'}, status=400)

    try:
        db_mol_obj = similaridade.objects.get(id=db_mol_id)
        db_sdf_path = db_mol_obj.sdf.path

        # A chamada da função agora retorna um dicionário
        alignment_data = perform_rdkit_alignment(user_sdf_path, db_sdf_path)
        
        with open(db_sdf_path, 'r', encoding='utf-8') as f:
            conteudo_db = f.read()

        return JsonResponse({
            'success': True,
            'conteudo_molecula_db': conteudo_db,
            'conteudo_molecula_alinhada': alignment_data["combined_sdf"],
            'user_mcs_indices': alignment_data["user_mcs_indices"],
            'db_mcs_indices': alignment_data["db_mcs_indices"],
            'shape_tanimoto': alignment_data["shape_tanimoto"],
            'pharm_similarity': alignment_data["pharm_similarity"], # <- NOVO
            'db_mol_nome': db_mol_obj.nome
        })
        
    except Exception as e:
        logging.exception("Erro inesperado na API get_alignment_data")
        return JsonResponse({'error': str(e)}, status=500)

def generate_conformers(input_sdf_path, output_sdf_path, num_confs=20):
    """
    Gera um conjunto de conformações de baixa energia para uma molécula a partir de um arquivo SDF.
    Salva todas as conformações que foram otimizadas com sucesso em um único arquivo SDF de saída.
    """
    logging.info(f"Iniciando geração de {num_confs} conformações para {os.path.basename(input_sdf_path)}")
    suppl = Chem.SDMolSupplier(input_sdf_path, removeHs=False)
    mol_template = suppl[0]
    
    if mol_template is None:
        raise ValueError("Não foi possível carregar a molécula do arquivo SDF de entrada.")

    mol_with_hs = Chem.AddHs(mol_template)
    
    # Gera múltiplas conformações aleatórias e retorna seus IDs
    cids = AllChem.EmbedMultipleConfs(mol_with_hs, numConfs=num_confs, randomSeed=42)
    logging.info(f"Geradas {len(cids)} conformações iniciais para otimização.")
    
    # Otimiza cada conformação e retorna uma lista de tuplas (energia, id) APENAS para as que tiveram sucesso
    results = AllChem.MMFFOptimizeMoleculeConfs(mol_with_hs)
    
    writer = Chem.SDWriter(output_sdf_path)
    if not results:
        raise ValueError("Nenhuma conformação pôde ser otimizada com sucesso.")

    logging.info(f"Otimização bem-sucedida para {len(results)} de {len(cids)} conformações.")

    # =================== CORREÇÃO DO BUG AQUI ===================
    #
    # Iteramos sobre os RESULTADOS da otimização, que contém apenas as conformações válidas.
    # O resultado é uma tupla (não_convergiu, energia), mas o ID da conformação é o seu índice.
    # A maneira mais segura é iterar sobre os IDs que sabemos que foram otimizados.
    #
    successful_conf_ids = [res[0] for res in results]

    for conf_id in successful_conf_ids:
        # Pega a energia correspondente para esta conformação
        energy = next((res[1] for res in results if res[0] == conf_id), "N/A")
        
        mol_with_hs.SetProp("energy_kcal_mol", f"{energy:.4f}")
        
        # Escreve a molécula especificando qual conformação (confId) usar
        writer.write(mol_with_hs, confId=conf_id)
    # ============================================================
    
    writer.close()
    
    if not os.path.exists(output_sdf_path) or os.path.getsize(output_sdf_path) == 0:
        raise IOError("A geração de conformações resultou em um arquivo de saída vazio.")

    logging.info(f"Conformações salvas com sucesso em {os.path.basename(output_sdf_path)}")
    return output_sdf_path

def split_multiconformer_sdf(multiconformer_sdf_path, output_dir):
    """
    Divide um arquivo SDF com múltiplas moléculas/conformações em arquivos SDF individuais.
    Retorna uma lista de caminhos para os arquivos individuais.
    """
    logging.info("Dividindo arquivo multi-conformacional em arquivos individuais...")
    suppl = Chem.SDMolSupplier(multiconformer_sdf_path, removeHs=False)
    output_paths = []
    for i, mol in enumerate(suppl):
        if mol:
            file_path = os.path.join(output_dir, f"conf_{i}.sdf")
            writer = Chem.SDWriter(file_path)
            writer.write(mol)
            writer.close()
            output_paths.append(file_path)
    logging.info(f"{len(output_paths)} arquivos de conformação individuais criados.")
    return output_paths

def consolidate_lsalign_results(all_results):
    """
    Recebe uma lista de resultados de múltiplas buscas e a consolida,
    mantendo apenas o melhor score (MAIOR TM-score) para cada molécula única.
    """
    logging.info("Consolidando resultados de múltiplas buscas...")
    best_hits = {}
    for result in all_results:
        mol_id = result['id']
        # --- LÓGICA DE RANQUEAMENTO ALTERADA AQUI ---
        # Se a molécula ainda não está na lista, ou se o novo TM-score é MAIOR que o anterior,
        # substitui o resultado pelo novo, que é melhor.
        if mol_id not in best_hits or result['tm_score'] > best_hits[mol_id]['tm_score']:
            best_hits[mol_id] = result
    
    consolidated_list = list(best_hits.values())
    # Ordena a lista final pelo TM-score, do maior para o menor
    consolidated_list.sort(key=lambda x: x['tm_score'], reverse=True) 
    
    logging.info(f"Resultados consolidados para {len(consolidated_list)} moléculas únicas, ranqueado por TM-score.")
    return consolidated_list

def calculate_lipinski_violations(mol):
    """
    Calcula o número de violações da Regra de Lipinski para uma molécula RDKit.
    Retorna um inteiro (0, 1, 2, 3 ou 4).
    """
    mol_with_hs = Chem.AddHs(mol) # Cálculos são mais precisos com hidrogênios
    
    mol_weight = Descriptors.MolWt(mol_with_hs)
    log_p = Descriptors.MolLogP(mol_with_hs)
    h_donors = Lipinski.NumHDonors(mol_with_hs)
    h_acceptors = Lipinski.NumHAcceptors(mol_with_hs)

    violations = 0
    if mol_weight > 500: violations += 1
    if log_p > 5: violations += 1
    if h_donors > 5: violations += 1
    if h_acceptors > 10: violations += 1
    
    return violations

    # Inicializa o catálogo de filtros de PAINS uma vez para maior eficiência
params = FilterCatalogParams()
params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
pains_catalog = FilterCatalog.FilterCatalog(params)

def has_pains_alert(mol):
    """
    Verifica se uma molécula RDKit contém alguma subestrutura PAINS.
    Retorna True se for um PAINS, False caso contrário.
    """
    if pains_catalog.HasMatch(mol):
        return True
    return False

def calculate_lipinski_violations(mol):
    """
    Calcula o número de violações da Regra de Lipinski para uma molécula RDKit.
    """
    mol_with_hs = Chem.AddHs(mol)
    mol_weight = Descriptors.MolWt(mol_with_hs)
    log_p = Descriptors.MolLogP(mol_with_hs)
    h_donors = Lipinski.NumHDonors(mol_with_hs)
    h_acceptors = Lipinski.NumHAcceptors(mol_with_hs)

    violations = 0
    if mol_weight > 500: violations += 1
    if log_p > 5: violations += 1
    if h_donors > 5: violations += 1
    if h_acceptors > 10: violations += 1
    
    return violations

def molecule_2d_image(request, mol_id):
    """
    Gera a imagem da estrutura 2D. Se falhar (ex: SMILES inválido),
    retorna uma imagem de placeholder com uma mensagem de erro.
    """
    try:
        mol_obj = get_object_or_404(similaridade, pk=mol_id)
        if mol_obj.smile:
            mol = Chem.MolFromSmiles(mol_obj.smile)
            if mol:
                # --- CORREÇÃO AQUI: Usamos "Draw" diretamente ---
                img = Draw.MolToImage(mol, size=(200, 200))
                
                response = HttpResponse(content_type="image/png")
                img.save(response, "PNG")
                return response
        
        # Se chegamos aqui, é porque mol_obj.smile está vazio ou o SMILES é inválido
        raise ValueError("SMILES ausente ou inválido")

    except Exception as e:
        logging.error(f"Não foi possível gerar imagem 2D para mol_id {mol_id}: {e}")
        
        # --- GERAÇÃO DA IMAGEM DE PLACEHOLDER ---
        img = Image.new('RGB', (200, 200), color = (230, 230, 230))
        d = ImageDraw.Draw(img)
        
        try:
            # Tenta usar uma fonte comum. Adapte o caminho se necessário.
            font = ImageFont.truetype("DejaVuSans.ttf", 14) 
        except IOError:
            font = ImageFont.load_default()
            
        d.text((10, 80), "Erro ao gerar imagem\n(SMILES inválido)", fill=(100, 100, 100), font=font, align="center")
        
        response = HttpResponse(content_type="image/png")
        img.save(response, "PNG")
        return response

def download_conformer(request):
    """
    Serve um arquivo de conformação de uma sessão temporária para download.
    Inclui uma verificação de segurança para evitar acesso a arquivos arbitrários.
    """
    file_path_relative = request.GET.get('path')
    if not file_path_relative:
        return HttpResponseBadRequest("Caminho do arquivo não fornecido.")

    # Medida de segurança: Garante que o caminho solicitado está dentro do diretório de sessões
    base_dir = os.path.join(settings.MEDIA_ROOT, 'temp_sessions')
    full_path = os.path.abspath(os.path.join(base_dir, file_path_relative))

    if not full_path.startswith(base_dir):
        return HttpResponseBadRequest("Acesso negado. Tentativa de acesso a um caminho inválido.")

    if os.path.exists(full_path):
        return FileResponse(open(full_path, 'rb'), as_attachment=True, filename=os.path.basename(full_path))
    else:
        return HttpResponse("Arquivo não encontrado.", status=404)

def pharmacophore_search(request):
    """
    Lida com a busca farmacofórica.
    GET: Mostra a página de upload.
    POST: Recebe a molécula modelo e (por enquanto) confirma o recebimento.
    """
    if request.method == 'POST':
        form = UploadSDFForm(request.POST, request.FILES)
        if form.is_valid():
            uploaded_file = request.FILES.get('file')
            if uploaded_file:
                # --- Lógica de Teste ---
                # Cria um diretório temporário para o arquivo, apenas para confirmar que o recebemos
                session_id = str(uuid.uuid4())
                temp_dir = os.path.join(settings.MEDIA_ROOT, "temp_sessions", session_id)
                os.makedirs(temp_dir, exist_ok=True)
                
                temp_file_path = os.path.join(temp_dir, uploaded_file.name)
                with open(temp_file_path, 'wb+') as dest:
                    for chunk in uploaded_file.chunks():
                        dest.write(chunk)
                
                logging.info(f"Busca Farmacofórica: Arquivo '{uploaded_file.name}' recebido e salvo em {temp_file_path}")
                
                # Na próxima etapa, aqui entrará a lógica de busca real.
                # Por agora, retornamos uma mensagem de sucesso simples.
                return HttpResponse(f"<h1>Sucesso!</h1><p>Arquivo {uploaded_file.name} recebido. A lógica de busca será implementada aqui.</p>")
        else:
            # Se o formulário for inválido, renderiza a página novamente com os erros
            return render(request, 'galeria/pharmacophore_search.html', {'form': form, 'error': 'Formulário inválido.'})

    # Se a requisição for GET, apenas mostra a página de upload
    form = UploadSDFForm()
    return render(request, 'galeria/pharmacophore_search.html', {'form': form})

def pharmacophore_search(request):
    """
    Lida com a busca farmacofórica.
    GET: Mostra a página de upload.
    POST: Executa a busca com o limiar definido pelo usuário e mostra os resultados.
    """
    if request.method == 'POST':
        form = UploadSDFForm(request.POST, request.FILES)
        if form.is_valid():
            uploaded_file = request.FILES.get('file')
            if not uploaded_file:
                return render(request, 'galeria/pharmacophore_search.html', {'form': form, 'error': 'Nenhum arquivo enviado.'})

            temp_dir = None
            try:
                # --- PASSO 1: LER PARÂMETROS E PREPARAR A QUERY ---
                
                # Lê o limiar do formulário, com um valor padrão de 0.6
                try:
                    threshold = float(request.POST.get('threshold', 0.6))
                    if not (0.1 <= threshold <= 1.0):
                        threshold = 0.6 # Garante que o valor esteja no intervalo
                except (ValueError, TypeError):
                    threshold = 0.6
                
                temp_dir = tempfile.mkdtemp()
                temp_file_path = os.path.join(temp_dir, uploaded_file.name)
                with open(temp_file_path, 'wb+') as dest:
                    for chunk in uploaded_file.chunks():
                        dest.write(chunk)
                
                suppl = Chem.SDMolSupplier(temp_file_path, removeHs=False)
                query_mol = suppl[0]
                if not query_mol:
                    raise ValueError("Arquivo SDF da query é inválido ou não pôde ser lido.")
                
                pharm_factory = Gobbi_Pharm2D.factory
                query_fp = Generate.Gen2DFingerprint(query_mol, pharm_factory)

                # --- PASSO 2: EXECUTAR A BUSCA ROBUSTA NO BANCO DE DADOS ---
                results = []
                all_molecules = similaridade.objects.exclude(pharm_fingerprint__isnull=True)
                
                for mol_db in all_molecules:
                    try:
                        db_fp = pickle.loads(mol_db.pharm_fingerprint)
                        similarity = DataStructs.TanimotoSimilarity(query_fp, db_fp)
                        
                        if similarity >= threshold:
                            results.append({
                                'mol': mol_db,
                                'similarity': similarity
                            })
                    except Exception as e:
                        logging.warning(f"Não foi possível processar o farmacóforo para mol_id {mol_db.id}. Erro: {e}")
                        continue
                
                results.sort(key=lambda x: x['similarity'], reverse=True)
                logging.info(f"Busca Farmacofórica: Encontrados {len(results)} hits com similaridade >= {threshold}")
                
                # --- PASSO 3: RENDERIZAR A PÁGINA DE RESULTADOS ---
                return render(request, 'galeria/pharmacophore_results.html', {
                    'results': results, 
                    'query_name': uploaded_file.name,
                    'threshold': threshold
                })

            except Exception as e:
                logging.error(f"Erro durante a busca farmacofórica: {e}")
                return render(request, 'galeria/pharmacophore_search.html', {'form': form, 'error': f"Ocorreu um erro: {e}"})
            
            finally:
                if temp_dir and os.path.exists(temp_dir):
                    shutil.rmtree(temp_dir)

    # Se a requisição for GET, apenas mostra a página de upload
    form = UploadSDFForm()
    return render(request, 'galeria/pharmacophore_search.html', {'form': form})