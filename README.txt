WEC:Weighted Edge Based Clustering to Identify Protein Complexes in PPI networks by incorporating Gene Expression Profile.

Author: Seketoulie Keretsu.

How to use:

By command line execution in windows operating system.

Step1: compile  Protein.java
Step2: compile Wec.java

Step3:  Run|     java Wec  <<ppi_input_file>> <<Gene_expresion_data_file >>  <<reference_complex_file>>   <<balanceT>>  <<weightT>>  <<filterT>>  <<enrichT>>
[eg>  java Wec collins2007.txt  gene.txt  sgd.txt  0.7   0.3   0.8   0.8  ]

[The input files should be kept in the working directory]


Parameters :

balanceT:  a value between [0-1] to balance the contibution of the similarity value and edge clustering coefficient value on the weight of the edges.

weight: The threshold weight to add a protein to a cluster to form a complex.

filterT :  The threshold value used to filter redundant complexes .

enrichT : a value to check if a highly connected protein can be added to a potential complex to enrich it.

Note:
PPI data:  The PPI network data should contain weighted interactions where the interactions are given by  Protein [space] Protein [space] weight (eg. proein1 Protein2  1.0 
Gene expression data: contains expression values of genes with time course.
sgd : a collection of complexes with each complexes containg protein names. .

