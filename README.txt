Files accompanying the manuscript:

Compound Functional Prediction Using Multiple Unrelated Morphological Profiling Assays

France Rose 1,*, Sreetama Basu 1,*, Elton Rexhepaj 2, 
Anne Chauchereau 3, Elaine Del Nery 2, and Auguste Genovesio 1

* These authors contributed equally to this work.
1 Institut de Biologie de l’Ecole Normale Supérieure, CNRS UMR8197,
Inserm, U1024, PSL Research University, Paris, France
2 Institut Curie, PSL Research University, Department of Translational
Research, Biophenics High-Content Screening Laboratory, F-75005,
Paris, France
3 Gustave Roussy, Prostate Cancer Group, Univ Paris-Sud, Inserm
UMR981, Villejuif, F-94805, France

Correspondence should be addressed to:
A.G. (auguste.genovesio@ens.fr)



========
CONTENTS
========

File					Description
---------------------------------------------------------------
README.txt				This file
MPM04_R1_v1.csv				Data file replicate 1 MPM04
MPM04_R2_v1.csv				Data file replicate 2 MPM04
MPM10_R1_v1.csv				Data file replicate 1 MPM10
MPM10_R2_v1.csv				Data file replicate 2 MPM10
MPM11_R1_v1.csv				Data file replicate 1 MPM11
MPM11_R2_v1.csv				Data file replicate 2 MPM11
MPM17_R1_v1.csv				Data file replicate 1 MPM17
MPM17_R2_v1.csv				Data file replicate 2 MPM17
MPM24_R1_v1.csv				Data file replicate 1 MPM24
MPM24_R2_v1.csv				Data file replicate 2 MPM24
MPM25_R1_v1.csv				Data file replicate 1 MPM25
MPM25_R2_v1.csv				Data file replicate 2 MPM25
MPM28_R1_v1.csv				Data file replicate 1 MPM28
MPM28_R2_v1.csv				Data file replicate 2 MPM28
MPM59_R1_v1.csv				Data file replicate 1 MPM59
MPM59_R2_v1.csv				Data file replicate 2 MPM59
MPM60_R1_v1.csv				Data file replicate 1 MPM60
MPM60_R2_v1.csv				Data file replicate 2 MPM60
Prestwick_R100_DoceLessN_R1_V1.csv	Data file replicate 1 Prostate assay
Prestwick_R100_DoceLessN_R2_V1.csv	Data file replicate 2 Prostate assay
Code_figure3_generatePredictions.m	Matlab script to train the model
Code_figure3_plot.m			Matlab script to plot the test results 



====
CODE
====

We use Matlab R2016b to train and run the network. 
The Machine Learning Matlab Toolbox is required. It can be found here: https://www.mathworks.com/help/stats/classificationensemble-class.html
Number of trees in each forest is set to 30 (cf Code_figure3_generatePredictions.m)



====
DATA
====

In data files, each line corresponds to a well treated with one compound (column "Contents"). Columns are both metadata - like compound, image name, targeted gene... - and averaged features. Each gene (column "TargetGene") has a unique identifier "TargetGeneNum" from 1 to 113. The "Count" column stores how many compounds in our dataset target the gene.



