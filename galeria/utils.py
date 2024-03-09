from rdkit import Chem
from rdkit.Chem import Descriptors

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


