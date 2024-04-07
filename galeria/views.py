from django.shortcuts import render, redirect, get_object_or_404
from .models import similaridade
from django.conf import settings
from django.http import HttpResponse
from .utils import DescritoresUtil
import logging
import os
from django.http import JsonResponse
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Descriptors, Lipinski, rdMolDescriptors
from .forms import UploadSDFForm
from django.db.models import Q 
import pandas as pd
from django.template.loader import render_to_string
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter
from django.core.paginator import Paginator
import tempfile

def index(request):
    return render(request, 'galeria/index.html')

def index(request):
    total_linhas = similaridade.objects.count()
    return render(request, 'galeria/index.html', {'total_linhas': total_linhas}) 

def moleculas(request):
    mol_list = similaridade.objects.order_by("id").filter(publicada=True)
    paginator = Paginator(mol_list, 96)  # 96 cards por página

    page_number = request.GET.get('page')
    mol = paginator.get_page(page_number)

    return render(request, 'galeria/moleculas.html',  {"cards": mol})

def moleculas_list_view(request):
    molecula_list = similaridade.objects.all()
    paginator = Paginator(molecula_list, 54) 

    page = request.GET.get('page')
    moleculas = paginator.get_page(page)

    print('Número de moléculas:', molecula_list.count())
    print('Número da página atual:', moleculas.number)
    print('Total de páginas:', moleculas.paginator.num_pages)

    return render(request, 'galeria/moleculas.html', {'moleculas': moleculas})

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

def ordenar_moleculas_view(request):
    ordenacao = request.GET.get('ordenacao')

    # Obtem todas as instâncias de Similaridade.
    moleculas = similaridade.objects.all()

    # Cria uma lista para armazenar as moléculas e seus dados extraídos.
    moleculas_com_dados = []
    for mol in moleculas:
        dados = extract_molecular_weight_and_logp(mol.sdf.path)
        if dados:
            moleculas_com_dados.append((mol, dados))
    
    # Se um critério de ordenação for especificado, ordene a lista.
    if ordenacao in ['peso_molecular', 'logp']:
        moleculas_com_dados.sort(key=lambda x: x[1][ordenacao])

    # Extraia as moléculas da lista de tuplas para passar para o template.
    moleculas_ordenadas = [mol_data[0] for mol_data in moleculas_com_dados]

    # Renderize o template com as moléculas ordenadas.
    return render(request, 'galeria/moleculas.html', {'cards': moleculas_ordenadas})

def get_mol_data(request, id):
    mol = get_object_or_404(Similaridade, id=id)
    sdf_data = mol.sdf.read()  # Lê o conteúdo do arquivo SDF
    return HttpResponse(sdf_data, content_type="chemical/x-mdl-sdfile")

def imagem(request, mol_id):
    mol = get_object_or_404(similaridade, pk=mol_id)
    return render(request, 'galeria/imagem.html', {"mol": mol})

def download_data_view(request, pk):
    objeto = get_object_or_404(similaridade, pk=mol_id)
    
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
    categoria = request.GET.get('categoria')
    nome_a_buscar = request.GET.get('buscar')
    mol = similaridade.objects.order_by("id").filter(publicada=True)

    if categoria:
        mol = mol.filter(categoria__icontains=categoria)

    if nome_a_buscar:
        # Verifica se o termo de busca é um número e pode ser um ID.
        try:
            nome_a_buscar_as_int = int(nome_a_buscar)
            # Primeiro tenta buscar pelo ID.
            mol = mol.filter(id=nome_a_buscar_as_int)
        except ValueError:
            # Se não for um número, faz a busca por outros campos de texto.
            mol = mol.filter(
                Q(nome__icontains=nome_a_buscar) | 
                Q(smile__icontains=nome_a_buscar) | 
                Q(categoria__icontains=nome_a_buscar) |  
                Q(url__icontains=nome_a_buscar)
            )
            
    return render(request, 'galeria/buscar.html', {"cards": mol}) 

def exemplo_molecula_view(request):
    # Criando uma molécula de exemplo (etanol)
    mol = Chem.MolFromSmiles("CCO")
    mol_with_hydrogens = Chem.AddHs(mol)

    # Preparando os dados para o template
    exemplo_data = {
        'atomos_antes': mol.GetNumAtoms(),
        'atomos_apos': mol_with_hydrogens.GetNumAtoms(),
        'peso_molecular': Descriptors.MolWt(mol_with_hydrogens),
        'logp': Descriptors.MolLogP(mol_with_hydrogens)
    }

    return render(request, 'galeria/teste.html', {'exemplo_data': exemplo_data})

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
        # Gera coordenadas 3D se não existirem (opcional)
        if not mol.GetNumConformers():
            AllChem.Compute2DCoords(mol)
            AllChem.EmbedMolecule(mol)
            AllChem.UFFOptimizeMolecule(mol)   

    # Converte para MolBlock
    molblock = Chem.MolToMolBlock(mol)

    return JsonResponse({"molblock": molblock})
 