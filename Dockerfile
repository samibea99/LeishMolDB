# Usando a imagem base do Python
FROM python:3.9-slim

# Definir o diretório de trabalho dentro do contêiner
WORKDIR /app

# Copiar o arquivo de requisitos para dentro do contêiner
COPY requirements.txt /app/

# Instalar as dependências necessárias
RUN pip install --no-cache-dir -r requirements.txt

# Copiar o código da aplicação para dentro do contêiner
COPY . /app/

# Coletar os arquivos estáticos do Django (precisa estar configurado o STATIC_ROOT no settings.py)
RUN python manage.py collectstatic --noinput

# Expor a porta onde o Django vai rodar
EXPOSE 8000

# Comando para rodar a aplicação Django
CMD ["python", "manage.py", "runserver", "0.0.0.0:8000"]
