# ALT-IN Tool

ALT-IN Tool is a software package that predicts alternative splining rewiring effects. This tool should be used with Python 2.7. [Anaconda](https://anaconda.org/anaconda/python) distriburion is recommended.

## INSTALLATION

For the latest version of ALT-IN Tool with all dependencies please consider using docker image. You would be able to run exaple simply by executing the following commands:

 '''bash
 docker run -it narykov/alt-in:latest bash
 cd ..
 python altintool.py test/interactors.tsv test/diabetes_ensembl_protein.fa test/string_protein.fa test/results.txt
 '''


## RUNNING

To actually make predictions you would need at least 2 files:
* Triplets consisting of the IDs of (Main isoform, Interacting partner, Alternatively spliced isoform) in tab separated format (.tsv)
* Fasta file(s) with protein sequences for the aforementioned protein IDs

To run our tool please use command in format

	python altintool.py [interactors_file.tsv] [fasta_file_1, fasta_file_2, ...] [output_file]

E.g., to run our test data you can use command

	python altintool.py test/interactors.tsv test/diabetes_ensembl_protein.fa test/string_protein.fa test/results.txt


## CONFIGURATION

You would need this step only if you want to change sequence alignment (currently EMBOSS implementation is used) or domain detection (currently performed by SUPERFAMILY).

We require two files to be accessible - pepstat from EMBOSS and ass3.pl from SUPERFAMILY.
If they are in the system path (i.e., accessible from the command line), no further actions required.

For each package that you cannot access from command line fill in corresponding line in 'config' file.
E.g., to provide path for INTERPRO fill in empty space on the second line of the file:
 INTERPRO_PATH='~/Downloads/interproscan-5.27-66.0'


## SUPPORT

For questions, please email to the onarykov@wpi.edu
