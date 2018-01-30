# AS-IN Tool

AS-IN Tool is a software package that predicts alternative splining rewiring effects.

## INSTALLATION

AS-IN Tool depends on 2 software packages - EMBOSS and INTERPRO.


##### EMBOSS 

Direct download - ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0.tar.gz

Could be downloaded from the [EMBOSS website](http://emboss.sourceforge.net/download/)
In our package version 6.6.0 was used. 



##### InterPro 

Direct download - ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.27-66.0/interproscan-5.27-66.0-64-bit.tar.gz

Could be downloaded from the [InterPro website](https://www.ebi.ac.uk/interpro/download.html)
In our package version 5.27-66.0 was used.


## CONFIGURING

We require two files to be accessible - pepstat from EMBOSS and insterproscan.sh from INTERPRO.
If they are in the system path (i.e., accessible from the command line), no further actions required.

For each package that you cannot access from command line fill in corresponding line in 'config' file.
E.g., to provide path for INTERPRO fill in empty space on the second line of the file:
 INTERPRO_PATH='~/Downloads/interproscan-5.27-66.0'


## RUNNING

To actually make predictions you would need at least 2 files:
⋅⋅*Triplets consisting of the IDs of (Main isoform, Interacting partner, Alternatively spliced isoform) in tab separated format (.tsv)
⋅⋅*Fasta file(s) with protein sequences for the aforementioned protein IDs

To run our tool pleas use command in format

	python [interactors_file.tsv] [output_file] [fasta_file_1, fasta_file_2, ...]

E.g., to run our test data you can use command

	python asintool.py test/interactors.tsv test/results.txt test/diabetes_ensembl_protein.fa test/string_protein.fa

## SUPPORT

For questions, please email to the onarykov@wpi.edu
