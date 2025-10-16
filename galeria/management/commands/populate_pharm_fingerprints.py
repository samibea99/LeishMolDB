# galeria/management/commands/populate_pharm_fingerprints.py

from django.core.management.base import BaseCommand
from galeria.models import similaridade
from rdkit import Chem
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D, Generate
import pickle

class Command(BaseCommand):
    help = 'Calcula e salva os fingerprints farmacofóricos (Gobbi_Pharm2D) para todas as moléculas.'

    def handle(self, *args, **options):
        self.stdout.write(self.style.SUCCESS('Iniciando o cálculo dos fingerprints farmacofóricos...'))
        
        # Inicializa a "fábrica" de farmacóforos correta
        pharm_factory = Gobbi_Pharm2D.factory
        
        molecules = similaridade.objects.all()
        total_molecules = molecules.count()
        updated_count = 0

        for i, mol_obj in enumerate(molecules):
            self.stdout.write(f'Processando molécula {i+1}/{total_molecules} (ID: {mol_obj.id})...', ending='')
            
            try:
                if not mol_obj.smile:
                    self.stdout.write(self.style.WARNING(' SMILES ausente. Pulando.'))
                    continue

                mol_rdkit = Chem.MolFromSmiles(mol_obj.smile)
                if not mol_rdkit:
                    self.stdout.write(self.style.WARNING(f' SMILES inválido. Pulando.'))
                    continue

                # Gera o fingerprint farmacofórico do tipo correto (SparseBitVect)
                fp = Generate.Gen2DFingerprint(mol_rdkit, pharm_factory)
                
                # Serializa o fingerprint para salvá-lo
                mol_obj.pharm_fingerprint = pickle.dumps(fp)
                mol_obj.save(update_fields=['pharm_fingerprint'])
                
                self.stdout.write(self.style.SUCCESS(' OK.'))
                updated_count += 1

            except Exception as e:
                self.stdout.write(self.style.ERROR(f' ERRO: {e}'))

        self.stdout.write(self.style.SUCCESS(f'\nProcesso concluído! {updated_count} de {total_molecules} fingerprints foram calculados e salvos.'))