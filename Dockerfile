# Usar uma imagem base oficial do Miniconda
FROM continuumio/miniconda3:4.12.0

# Instalar dependências de sistema necessárias (como o 'zip')
RUN apt-get update && apt-get install -y zip

# Definir o diretório de trabalho
WORKDIR /app

# Instalar Open Babel e Pandas usando Conda (o método mais robusto)
RUN conda install -c conda-forge openbabel pandas -y

# Copiar o arquivo de dependências do Django
COPY requirements.txt .

# Instalar as dependências do Django com pip
RUN pip install -r requirements.txt

# Copiar todo o resto do seu projeto
COPY . .