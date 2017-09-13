# DOCKERFILE for trifusion-cli https://github.com/ODiogoSilva/TriFusion
FROM ubuntu:16.04
MAINTAINER Diogo N. Silva, diogosilva@medicina.ulisboa.pt

RUN apt-get update && apt-get -y install \
	python-all-dev \
	python-setuptools \
	python-configparser \
	python-matplotlib \
	python-numpy \
	python-psutil \
	python-scipy \
	python-seaborn \
	python-pandas \
	python-pip

WORKDIR /home/user

RUN pip install trifusion

RUN bash
