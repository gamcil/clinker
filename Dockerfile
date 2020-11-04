FROM ubuntu:20.04

RUN apt update && \
    apt install -y python3 python3-pip

ADD . /usr/local/clinker/
RUN pip3 install /usr/local/clinker/
