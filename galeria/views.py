from django.shortcuts import render, redirect, get_object_or_404
from .models import similaridade
from django.conf import settings
from django.http import HttpResponse
from .utils import DescritoresUtil
import logging
import os
from django.http import JsonResponse
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Descriptors
from .forms import UploadSDFForm
from django.db.models import Q 
import pandas as pd
from django.template.loader import render_to_string
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter

def index(request):
    return render(request, 'galeria/index.html')

def index(request):
    total_linhas = similaridade.objects.count()
    return render(request, 'galeria/index.html', {'total_linhas': total_linhas}) 

def moleculas(request):
    mol = similaridade.objects.order_by("id").filter(publicada=True)
    return render(request, 'galeria/moleculas.html',  {"cards": mol})

def extract_molecular_weight_and_logp(sdf_path):
    suppl = Chem.SDMolSupplier(str(sdf_path))
    mol = next((m for m in suppl if m is not None), None)  # Apenas a primeira molécula válida.

    if mol is None:
        return None

    mol_with_hydrogens = Chem.AddHs(mol)
    mol_weight = Descriptors.MolWt(mol_with_hydrogens)
    log_p = Descriptors.MolLogP(mol_with_hydrogens)

    return {'peso_molecular': mol_weight, 'logp': log_p}

def ordenar_moleculas_view(request):
    ordenacao = request.GET.get('ordenacao')

    # Obtenha todas as instâncias de Similaridade.
    moleculas = similaridade.objects.all()

    # Crie uma lista para armazenar as moléculas e seus dados extraídos.
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
    
    # Aqui você deve obter os dados do objeto e formatá-los como desejar
    dados_para_download = objeto.dados_de_interesse()

    # Agora, vamos criar um arquivo temporário para enviar como resposta
    filename = "dados_para_download.txt"
    with open(filename, 'w') as f:
        f.write(dados_para_download)

    # Retornar o arquivo como resposta para download
    response = FileResponse(open(filename, 'rb'))
    response['Content-Disposition'] = 'attachment; filename="dados_para_download.txt"'
    return response


def buscar(request):    
    categoria = request.GET.get('categoria')
    mol = similaridade.objects.order_by("id").filter(publicada=True)

    if categoria:
        mol = mol.filter(categoria__icontains=categoria)

    if "buscar" in request.GET:
        nome_a_buscar = request.GET['buscar']
        if nome_a_buscar:
            mol = mol.filter(Q(nome__icontains=nome_a_buscar) | 
                             Q(smile__icontains=nome_a_buscar) | 
                             Q(categoria__icontains=nome_a_buscar) | 
                             Q(sdf__icontains=nome_a_buscar) | 
                             Q(url__icontains=nome_a_buscar))
            
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
    data = {'Molecules': [], 'MolecularWeights': [], 'LogPs': []}
    suppl = Chem.SDMolSupplier(str(sdf_path))

    for idx, mol in enumerate(suppl):
        if mol:
            mol_with_hydrogens = Chem.AddHs(mol)  # Adiciona hidrogênios
            smiles = Chem.MolToSmiles(mol_with_hydrogens)
            mol_weight = Descriptors.MolWt(mol_with_hydrogens)
            log_p = Descriptors.MolLogP(mol_with_hydrogens)

            data['Molecules'].append(smiles)
            data['MolecularWeights'].append(mol_weight)
            data['LogPs'].append(log_p)

    return data

def molecula_view(request, id):
    # Obter a instância do modelo com base no ID
    similaridade_instancia = get_object_or_404(similaridade, id=id)
    
    # O caminho do arquivo SDF é obtido pelo atributo 'path' do FileField
    sdf_path = similaridade_instancia.sdf.path

    # Chamar a função para extrair os dados do arquivo SDF
    data = extract_data_from_sdf(sdf_path)

    # Agora só vamos zipar as duas listas: pesos moleculares e LogPs
    moleculas_data = [
        {
            'peso_molecular': mw,
            'logp': logp
        }
        for mw, logp in zip(data['MolecularWeights'], data['LogPs'])
    ]

    # Passar os dados extraídos para o contexto do template
    return render(request, 'galeria/imagem.html', {'moleculas_data': moleculas_data, 'mol': similaridade_instancia})



def download_pdf(request, id):
    # Busca a molécula pelo ID
    molecula = similaridade.objects.get(id=id)

    # Cria uma resposta HTTP do tipo PDF
    response = HttpResponse(content_type='application/pdf')
    response['Content-Disposition'] = f'attachment; filename="mol_{molecula.id}.pdf"'
    response['X-Content-Type-Options'] = 'nosniff'

    # Cria o objeto PDF
    p = canvas.Canvas(response, pagesize=letter)

    # Define um título para o PDF
    p.setFont("Helvetica-Bold", 20)
    p.drawString(100, 750, f"Dados da Molécula: {molecula.nome}")

    # Insere os dados da molécula no PDF
    p.setFont("Helvetica", 12)
    y = 730
    attributes = [
        ("Nome", molecula.nome),
        ("Categoria", molecula.categoria),
        ("URL", molecula.url),
        ("SMILES", molecula.smile),
    ]

    for attr, value in attributes:
        p.drawString(100, y, f"{attr}: {value}")
        y -= 20

    # Adiciona um link para o arquivo SDF
    y -= 20  # Decrementa y para posicionar o link abaixo dos outros dados
    sdf_url = request.build_absolute_uri(molecula.sdf.url)
    p.drawString(100, y, "SDF: Clique aqui para baixar")
    p.linkURL(sdf_url, (100, y - 10, 300, y + 10), thickness=1)

    p.showPage()
    p.save()

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

            # Aqui é onde uploaded_mol_block deve ser definido
            uploaded_mol_block = uploaded_file.read()

            # Verifica se o conteúdo é binário e decodifica se necessário
            if isinstance(uploaded_mol_block, bytes):
                uploaded_mol = Chem.MolFromMolBlock(uploaded_mol_block.decode('utf-8'))
            else:
                uploaded_mol = Chem.MolFromMolBlock(uploaded_mol_block)

            # Continua o processamento apenas se uploaded_mol foi criado com sucesso
            if uploaded_mol:
    # Processamento e cálculo da similaridade
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

                if 0.9 <= similarity <= 1.0:
                    # Adiciona a similaridade e o objeto similaridade à lista
                    similarities_with_ids.append((similarity, similaridade_obj))

    return render(request, 'galeria/resultado_similaridade.html', {'similarities_with_ids': similarities_with_ids})


def get_molecule_data(request):
    molecule_id = request.GET.get('id')
    similaridade_obj = similaridade.objects.get(id=molecule_id)  # Obtém o objeto específico
    sdf_path = similaridade_obj.sdf.path  # Obtém o caminho do arquivo no sistema de arquivos

    # O RDKit espera um caminho de arquivo ou um objeto de arquivo, então usamos o caminho
    suppl = Chem.SDMolSupplier(str(sdf_path))
    mols = [mol for mol in suppl if mol is not None]

    # Assumindo que você quer trabalhar com a primeira molécula
    mol = mols[0] if mols else None

    if mol:
        # Gera coordenadas 3D se não existirem (opcional, se seu SDF não tiver coordenadas 3D)
        if not mol.GetNumConformers():
            AllChem.Compute2DCoords(mol)
            AllChem.EmbedMolecule(mol)
            AllChem.UFFOptimizeMolecule(mol)  # Otimização UFF; alternativa: otimização MMFF

    # Converte para MolBlock
    molblock = Chem.MolToMolBlock(mol)

    return JsonResponse({"molblock": molblock})

def get_SDF_files(request, filename):
    file_path = os.path.join(settings.BASE_DIR, filename)
    print("oiiiii", file_path)
    if os.path.exists(file_path):
        return FileResponse(open(file_path, 'rb'), as_attachment=True)
    else:
        return HttpResponseNotFound("O arquivo não foi encontrado.")