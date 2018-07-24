FROM ubuntu:xenial
LABEL maintainer "onarykov@wpi.edu"

RUN apt-get update && apt-get upgrade -y -q && apt-get install -y -q \
    software-properties-common \
	wget git python-software-properties debconf-utils
RUN export LANG=en_US.UTF-8


# Install JAVA 8
#RUN add-apt-repository ppa:webupd8team/java
#RUN apt-get update
#RUN echo oracle-java8-installer shared/accepted-oracle-license-v1-1 select true | /usr/bin/debconf-set-selections
#RUN apt-get install -y -q oracle-java8-installer
#RUN apt-get install -y -q oracle-java8-set-default
#
#ENV JAVA_HOME=/usr/lib/jvm/java-8-oracle
#ENV CLASSPATH=/usr/lib/jvm/java-8-oracle/bin


# Install InterProScan
#RUN wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.27-66.0/interproscan-5.27-66.0-64-bit.tar.gz
	#wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.27-66.0/interproscan-5.27-66.0-64-bit.tar.gz.md5 && \
	#md5sum -c interproscan-5.27-66.0-64-bit.tar.gz
#RUN tar -xvzf interproscan-5.27-66.0-64-bit.tar.gz
#ENV INTERPRO_PATH=/interproscan-5.27-66.0
ENV INTERPRO_PATH=test

#Install pip
RUN apt-get install python-pip python-dev build-essential -y

#Install dependencies
RUN apt-get install libx11-dev -y

# Install EMBOSS
RUN wget ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0.tar.gz
RUN tar -xvzf EMBOSS-6.6.0.tar.gz
WORKDIR ./EMBOSS-6.6.0
RUN ./configure
RUN make
WORKDIR ../
ENV EMBOSS_PATH=/EMBOSS-6.6.0/emboss


# Install ALT-IN Tool
#RUN git clone https://github.com/korkinlab/altintool.git
COPY . .

#FROM python:2.7-slim
#
# Installing ALT-IN Tool dependencies
RUN pip install --trusted-host pypi.python.org -r requirements.txt

# Installing Web-Runner
RUN git clone https://github.com/AlexandrNP/web-runner.git
RUN pip install --trusted-host pypi.python.org -r web-runner/requirements.txt
COPY webrunner-config/config web-runner/config
WORKDIR /web-runner


# Run server
#CMD ["./app.py"]

#Test run
#CMD ["python", "altintool.py", "test/interactors.tsv", "test/results.txt", "test/diabetes_ensembl_protein.fa", "test/string_protein.fa"]

#FROM ubuntu:xenial
#RUN python altintool.py test/interactors.tsv test/results.txt test/diabetes_ensembl_protein.fa test/string_protein.fa
#RUN echo $(ls .)
