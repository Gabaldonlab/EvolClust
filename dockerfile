FROM continuumio/anaconda

MAINTAINER Marina Marcet-Houben version 1.0

RUN apt-get install mcl
RUN mkdir -p /evolclust
WORKDIR /evolclust
COPY . /evolclust
ENTRYPOINT ["python","evolclust.py"]
