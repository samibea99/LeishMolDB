from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
import logging
import os

class DescritoresUtil:
    @staticmethod
    def calcular_descritores(sdf):
        mol = Chem.MolFromMolBlock(sdf)

        if mol:
            print("Molécula RDKit:", mol)
            logP = Descriptors.MolLogP(mol)
            massa_molecular = Descriptors.MolWt(mol)
            tpsa = Descriptors.TPSA(mol)

            print("LogP:", logP)
            print("Massa Molecular:", massa_molecular)
            print("TPSA:", tpsa)

            return logP, massa_molecular, tpsa
        else:
            return None

logging.basicConfig(level=logging.INFO, format="%(asctime)s - [%(levelname)s] - %(message)s")

def process_user_uploaded_file(uploaded_file, temp_dir):
    """
    Processa o arquivo SDF/MOL2 enviado pelo usuário.
    Gera versão SDF limpa e com coordenadas 3D (necessária para o RDKit O3A e 3Dmol.js).
    Retorna caminhos (mol2_path, sdf_path).
    """

    try:
        # Cria os caminhos temporários
        base_name = os.path.splitext(uploaded_file.name)[0]
        sdf_path = os.path.join(temp_dir, f"{base_name}.sdf")
        mol2_path = os.path.join(temp_dir, f"{base_name}.mol2")

        # Salva o arquivo original
        input_path = os.path.join(temp_dir, uploaded_file.name)
        with open(input_path, "wb") as dest:
            for chunk in uploaded_file.chunks():
                dest.write(chunk)
        logging.info(f"Arquivo salvo em: {input_path}")

        # Detecta formato (SDF ou MOL2)
        ext = os.path.splitext(uploaded_file.name)[1].lower()
        if ext not in [".sdf", ".mol2"]:
            raise ValueError(f"Formato de arquivo não suportado: {ext}")

        # Carrega molécula com RDKit
        mol = None
        if ext == ".sdf":
            suppl = Chem.SDMolSupplier(input_path, removeHs=False)
            mol = suppl[0] if suppl and len(suppl) > 0 else None
        elif ext == ".mol2":
            mol = Chem.MolFromMol2File(input_path, removeHs=False, sanitize=True)

        if mol is None:
            raise ValueError("Não foi possível carregar a molécula enviada.")

        # Gera coordenadas 3D se não existirem
        if mol.GetNumConformers() == 0:
            logging.info("Nenhum confôrmero encontrado — gerando coordenadas 3D.")
            AllChem.EmbedMolecule(mol, randomSeed=42)
            AllChem.UFFOptimizeMolecule(mol)

        # Salva versões SDF e MOL2
        Chem.MolToMolFile(mol, mol2_path)
        Chem.MolToMolFile(mol, sdf_path)
        logging.info(f"Arquivos gerados: {mol2_path}, {sdf_path}")

        return mol2_path, sdf_path

    except Exception as e:
        logging.exception("Erro ao processar o arquivo enviado pelo usuário.")
        raise e


