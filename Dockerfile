FROM tensorflow/tensorflow:latest-gpu

WORKDIR /work

COPY requirements.txt .
RUN pip install -r requirements.txt
RUN apt-get update && apt-get install graphviz -y

USER dev
