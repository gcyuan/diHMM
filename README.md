diHMM v0.1 beta
===============

Overview
--------

diHMM stands for Hierarchical Hidden Markov Model. diHMM is a novel computational method for finding chromatin states at multiple scales. The model takes as input a multidimensional set of histone modifications for several cell types and classifies the genome into a preselected number of nucleosome-level and domain-level hidden states. The work is done by Eugenio Marco, Wouter Meuleman, Luca Pinello, Manolis Kellis and Guo-Cheng Yuan.

Systems Requirements
--------------------

diHMM is independent of operating systems because it is written in Matlab. Basic requirements for running diHMM include MATLAB and the Bioinformatics toolbox. diHMM includes the following packages, 'cbrewer', 'cm_and_cb_utilities', 'export_fig', 'freezeColors', 'kcenter'.


Usage
-----

Unzip the package. Change the current directory in Matlab to the 'diHMM' folder containing the code organized into subfolders. In order to run the programs, 'diHMM' and its sudirectories have to be added to the path. This can be achieved with the folowing command at the Matlab prompt

```
>> addpath(genpath(pwd))
```

The data to be analyzed with diHMM has to be placed in the folder 'data'. The data has to be preprocessed so that it is binarized. diHMM does not include a binarization function at this moment, and we have used ChromHMM to binarize the data (ChromHMM was developed by Jason Ernst and Manolis Kellis, and it is available at http://compbio.mit.edu/ChromHMM/). diHMM comes already with binarized ENCODE data for the cell lines H1-hesc, GM12878 and K562, and nine histone marks. 

To train a model use the function

```
>> runDiHMMTrain
```

By default, the model trained has 30 nucleosome-level states and 30 domain-level states. To reduce the training time a smaller number of states can be given as an input to the function. 


For example: To train the model with 20 nucleosome-level states and 10 donain-level states use the command.

```
>> runDiHMMTrain({'GM12878', 'H1hesc', 'K562'}, 'hg19', 20, 10, 'modelGHKChr17',  1)
```
See the file runDiHMMTrain.m in 'main' for a description of the input parameters. After training the model, the model and bed files to visualize the states in a genome browser are saved in the 'results' folder, subfolder 'modelGHKChr17', and the model parameters are displayed.

Nuclesome and domain states can be labeled with different functional categories, eg, Promoter, Weak Enhancer, Super-Enhancers, etc. To see available categories use the commands

```
>> plotDomainAnnotations
>> plotNucleosomeAnnotations
```

Initially, all states are labeled with the white annotations. State labels can be changed using the function changeAnnotations. For example

```
>> changeAnnotations({'GM12878', 'H1hesc', 'K562'}, 20, 10, 'modelGHKChr17', 1, 'bin',[1, 4], [1, 3])
```

will change nucleosome level (bin) states 1 and 4 to categories 1, Active Promoter, and 3, Bivalent Promoter. Analogously, for domains we can use

```
>> changeAnnotations({'GM12878', 'H1hesc', 'K562'}, 20, 10, 'modelGHKChr17', 1, 'domain',[1, 2, 3], [7, 5, 3])
```

Basic model statistics like state coverage, domain enrichments into nucleosome-level states, histone marks predominant in each domain, or state sizes can be calculated with the function changeAnnotations. For example

```
>> calculateBasicModelStatistics({'GM12878', 'H1hesc', 'K562'}, 20, 10, 'modelGHKChr17', 1)
```

Figures are save in the 'figures' folder in 'results'.

