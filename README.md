# ALT-IN
-------

ALT-IN is the instrument for predicting disruptions in protein-protein interactions induced by alternative splicing.
We recommend usage of Docker version of our tool which is available on DockerHub repository (https://hub.docker.com/r/narykov/alt-in)
and has a convenient interface (altin_docker.py)

[![DOI](https://zenodo.org/badge/150632845.svg)](https://zenodo.org/badge/latestdoi/150632845)

## INSTALLATION
-------

#### DOCKER INTERFACE

The ONLY file you would need to run Docker image is **altin_docker.py**.
The prerequisite is existing Docker installation, the image would be downloaded automatically. You would be able to run example just by executing the following commands:

 > python3 altin_docker.py -i test/interactors.tsv -o res.txt -f test/diabetes_ensembl_protein.fa test/string_protein.fa

Download would take approximately 10 min on stable internet connection.
This is an interface script for Docker container. The ONLY dependency for it is preinstalled Docker with corresponding permissions.

###### OPTIONS

* **-h, --help**

    Prints brief usage information.

* **-i, --interactors**

    A tsv file containing triplets in format "reference_isoform", "interactor", "alternative_isoform". You can use test/interactors.tsv file for the reference.

* **-d, --output-dir**

    Output directory where output would be stored. Permissions should allow for Docker to bind to this location. Default value is current directory

* **-o, --output-file**

    Name of the output file. Would be stored in directory specified by "--output-dir" option.

* **-f, --fastafiles**

    A list of fasta files, which should contain sequences for all IDs occuring in interactors file. Number of files you may specify is not limited.

#### STANDALONE

Read the following information if you wish to install ALT-IN Tool without Docker as a standalone package.
Standalone version of ALT-IN Tool depends on 2 software packages - EMBOSS and INTERPRO.
Installation was tested on Ubuntu 16.04 with Python 2.7 version of Anaconda.

###### EMBOSS 

Direct download - ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0.tar.gz

Could be downloaded from the [EMBOSS website](http://emboss.sourceforge.net/download/). In our package version 6.6.0 was used. 



###### InterPro 

Direct download - ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.27-66.0/interproscan-5.27-66.0-64-bit.tar.gz

Could be downloaded from the [InterPro website](https://www.ebi.ac.uk/interpro/download.html). In our package version 5.27-66.0 was used.



#### CONFIGURATION

You would need this step only if you want to change sequence alignment (currently EMBOSS implementation is used) or domain detection (currently performed by SUPERFAMILY).

We require two files to be accessible - pepstat from EMBOSS and ass3.pl from SUPERFAMILY.
If they are in the system path (i.e., accessible from the command line), no further actions required.

For each package that you cannot access from command line fill in corresponding line in 'config' file.
E.g., to provide path for INTERPRO fill in empty space on the second line of the file:
 INTERPRO_PATH='~/Downloads/interproscan-5.27-66.0'


#### RUNNING

To run our tool please use command in format

	python altintool.py [interactors_file.tsv] [fasta_file_1, fasta_file_2, ...] [output_file]

E.g., to run our test data you can use command

	python altintool.py test/interactors.tsv test/diabetes_ensembl_protein.fa test/string_protein.fa test/results.txt



## INPUT
-------

To make predictions you would need at least 2 files:
* Triplets consisting of the IDs of (Main isoform, Interacting partner, Alternatively spliced isoform) in tab separated format (.tsv)
* Fasta file(s) with protein sequences for the aforementioned protein IDs



## SUPPORT
-------

For questions, please email to the onarykov@wpi.edu
