WEC: Weighted Edge Based Clustering to Identify Protein Complexes in PPI networks by incorporating Gene Expression Profile.


Author: Seketoulie Keretsu
////////////////////////////////////////
//////////////////WEC/////////////
////////////////////////////////////////
The EccProject.jar file  of WEC can be executed on system with  jdk1.8 (and other compatible jdk) installations.
The file has been successfully executed and tested on Windows 8.1 platform, running on Intel core i3 processor with 4 GB RAM.

The algorithm works by clustering Proteins based on the edge clustering coefficent and the similarity in the gene expression. For a protein to be added to a cluster, it should have high clustering coefficient and high similarity in its epression. A filtering process will filter the highly overlapping complexes based on a filtereing threshold value. The complex is finally enriched by adding protein from its neighborhood that has high neighborhood affinity.
////////////////////////////////////
/////////INPUT DATA: ///////////////
////////////////////////////////////
1. weighted PPI network data containing Protein-Protein Interactions and the weight of the    interaction.
   Example: proein1<tab>  protein2<tab>  0.5.   protein1 and protein2 are protein systematic names seperated by a tab or    space, followed by the weight of the edge or interaction.[collins2007.txt] 
For unwaighted PPI Network use 1 to indicate binary interaction:   Protein1<tab> Protein2 <tab> 1;

2.Gene Expression Profile: Table containing proteins and its corresponding expresion with    respect to time course.[gene.txt]
   Example: YNR066C  0.091940984  0.06936416  0.147559598  0.12710911  0.089786761             
3. Gold standard complexes: A reference benchmark complex set to compare the predicted       complexes with some real complexes.
///////////////////////////////////// 
/////////INPUT PARAMETERS: //////////
/////////////////////////////////////
1. Balance Factor(or Lambda): a value between [0,1] to balance the contribution of clustering coefficient   and similarity value to the weight of the Edge.

2. Edge weight Threshold: A protein is added to a cluster only if the weight of the edge       connecting to the cluster is greater than the Edge weight threshold value.
3. Enrichment threshold : To add protein with high connectivity with the complex.

4. Filtering Threshold  : To filter off high overlapping complexes.

//////////INPUT VALUES Example1:///////////
PPI network: collins2007.txt
GEP data   : gene.txt
Gold Standard File: sgd.txt
threshold parameters
Balance factor: 0.8
Edge Weight   : 0.7 
Enrichment    : 0.8
Filtering     : 0.9

//////////INPUT VALUES Example2:///////////
PPI network: HumanPPI.txt
GEP data   : Humangene.txt
Gold Standard File: HumanBenchmarkComplex.txt
threshold parameters
Balance factor: 0.9
Edge Weight   : 0.5 
Enrichment    : 0.9
Filtering     : 0.9

////////////////////////////////////////////
/////////////OUTPUT ////////////////////////
////////////////////////////////////////////
output shown on the interface:
No. of distinct proteins used.
No. of interactions that took part in the network.
No. of predicted Complexes.
No. of complexes predicted that has a match in the real benchmark complex data.

Output Files:  CHECK CURRENT WORKING DIRECTORY FOR OUPUT
The predicted complexes are shown in the file [output_PredictedComplex.txt]
The predicted complexes that has a match in the real complex are shown in [output_DetectedComplex]

NOTE: This Interface is expected to be submitted along with the Paper of WEC after adding some additional design and components.
