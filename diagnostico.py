import sys
import os

print("--- INICIANDO DIAGNÓSTICO DO AMBIENTE ---")
print("-" * 40)

# 1. Verificar qual executável do Python está sendo usado
print(f"[*] Executável do Python:\n    {sys.executable}\n")

# 2. Verificar a versão do Python
print(f"[*] Versão do Python:\n    {sys.version}\n")

# 3. Verificar os caminhos de busca de módulos
print("[*] Caminhos de busca do Python (sys.path):")
for path in sys.path:
    print(f"    - {path}")
print("\n")

# 4. Tentar importar e inspecionar o RDKit
print("-" * 40)
print("[*] Investigando a biblioteca RDKit...")
try:
    import rdkit
    from rdkit import Chem

    print("[+] SUCESSO: RDKit importado.")
    
    # Verificar a versão do RDKit
    print(f"    - Versão do RDKit: {rdkit.__version__}")
    
    # Verificar DE ONDE o RDKit está sendo importado
    print(f"    - Arquivo do RDKit: {rdkit.__file__}")

    # Verificar DE ONDE o submódulo Chem está sendo importado
    print(f"    - Arquivo do rdkit.Chem: {Chem.__file__}")

    # Verificar se a função 'MolToMol2Block' realmente existe
    if hasattr(Chem, 'MolToMol2Block'):
        print("\n[+] SUCESSO: A função 'MolToMol2Block' FOI ENCONTRADA em rdkit.Chem.")
    else:
        print("\n[-] FALHA: A função 'MolToMol2Block' NÃO FOI ENCONTRADA em rdkit.Chem.")
        print("    --> Isto confirma a causa do erro 'AttributeError'.")

except ImportError as e:
    print(f"\n[-] FALHA CRÍTICA: Não foi possível importar o RDKit.")
    print(f"    Erro de importação: {e}")
except Exception as e:
    print(f"\n[-] FALHA INESPERADA: Ocorreu um erro ao inspecionar o RDKit.")
    print(f"    Erro: {e}")

print("-" * 40)
print("--- DIAGNÓSTICO CONCLUÍDO ---")
