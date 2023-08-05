# HElA: A fast and accurate tool for detection of Helitron-like elements.
Helitron-like elements including Helitron, Helentron and Helitron2 are one of the class 2 transposons. They have been found in a diverse range of species and seem to play major roles in the evolution of host genomes. Although they have benn well known for more than twenty years, they were widely admitted as one of few remaining transposons difficult to identify. Here, we propose HELA (Helitron-like elements annotator) which is a fast and accurate tool for detection of Helitron-like elements.

# Table of contents
- [Dependencies](#dependencies)
- [Installation](#installation)
  * [conda](#conda)
  * [mamba](#mamba)
  * [Manual installation](#manual-installation)
- [Usage](#usage)
- [To contact us](#to-contact-us)

# Dependencies
```
- python = 3.9.0
- r-base = 4.1
- biopython
- pybedtools = 0.9.0 =py39hd65a603_2
- r-bedtoolsr
- r-seqinr = 4.2_16 = r41h06615bd_0
- bedtools = 2.30.0
- dialign2 = 2.2.1
- mafft
- cd-hit = 4.8.1
- blast = 2.2.31
- emboss = 6.6.0
- hmmer = 3.3.2
- genometools-genometools = 1.6.2 = py39h58cc16e_6
- rnamotif = 3.1.1
```
# Installation
## conda
```
#create the hela environment
conda create -n hela
#activate the hela environment
conda activate hela
# installation 
conda install -c zhenlisme hela -c bioconda -c conda-forge
conda deactivate
```
## mamba
```
#create the hela environment
mamba create -n hela
#activate the hela environment
mamba activate hela
# install 
mamba install -c zhenlisme hela -c bioconda -c conda-forge
mamba deactivate
```
## manual installation
Before installation , you need to be sure that all dependencies have been installed in your computer and their pathes have been added into your environmental variables. All dependencies except rnamotif could be installed via conda/mamba.  
1. download the latest hela package.  
`git clone https://github.com/Zhenlisme/HELA.git`
2. switch to the source code dorectory that you cloned at the last step.  
   `cd HELA/`  
3. run configure file.  
   `bash configure.sh`
4. You can find HELA in bin/ directory.
# Usage
### 1. Activate the hela conda environment (for conda/mamba installation)  
`conda activate hela`  
### 2. Check the HELA binary  
`HELA -h`
```
usage: HELA [-h] -g GENOME [-w WINDOW] [-dm DISTANCE_DOMAIN] [-pt {0,1}] [-is1 {0,1}] [-is2 {0,1}] [-sim_tir {100,90,80}] [-p PVALUE] -o OPDIR
            [-n PROCESS]

HELA can detect and classify different variants of Helitron-like elements: Helitron, Helentron and Helitron2. email us: zhen.li3@universite-paris-
saclay.fr

optional arguments:
  -h, --help            show this help message and exit
  -g GENOME, --genome GENOME
                        The reference genome.
  -w WINDOW, --window WINDOW
                        To check terminal markers within a given window bp upstream and downstream of ORF ends, default is 10 kb.
  -dm DISTANCE_DOMAIN, --distance_domain DISTANCE_DOMAIN
                        The distance between HUH and Helicase domain, default is 2500.
  -pt {0,1}, --pair_helitron {0,1}
                        The 5' and 3' terminal signal pairs should come from the same autonomous Helitrons or not. 0: no, 1: yes. default no.
  -is1 {0,1}, --IS1 {0,1}
                        require the insertion site of Helitron seeds as AT. 0: no, 1: yes. default no.
  -is2 {0,1}, --IS2 {0,1}
                        Set the insertion site of Helentron/Helitron2 seeds as TT. 0: no, 1: yes. default no.
  -sim_tir {100,90,80}, --simtir {100,90,80}
                        Set the simarity between short inverted repeats(TIRs) of Helitron2/Helentron. Default 100.
  -p PVALUE, --pvalue PVALUE
                        The p-value for fisher test.
  -o OPDIR, --opdir OPDIR
                        The output directory.
  -n PROCESS, --process PROCESS
                        Number of threads to be used.
```
### 3. Perform a test run of HELA  
##### Here we will use the genome of Xenopus tropicalis as an example. Perform the following code. (note: the `genome.fa` is the only input which is the genome file of Xenopus tropicalis in fasta format. The main result could be found in file **RC.representative.bed** in **HELA_opt** directory.)   
`HELA -g genome.fa -pt 1 -is1 1 -is2 1 -sim_tir 90 -p 1e-5 -n 30 -o HELA_opt`  
##### Explanation for RC.representative.bed
There are 11 columns in RC.representative.bed file:  
|chrm-id|start|end|subfamily|occurence|strand|pvalue|TS_blastn_identity|variant|type|name|
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
|NC_030677.2|20246175|20247385|Helentron_left_73-Helentron_right_74|79|+|2.0432e-176|104|Helentron|nonauto|insertion_Helentron_nonauto_1|
|NC_030677.2|23640161|23640743|Helentron_left_73-Helentron_right_74|79|+|2.0432e-176|74.3|Helentron|nonauto|insertion_Helentron_nonauto_2|
......
# To contact us
For any questions, please email us: zhen.li3@universite-paris-saclay.fr or nicolas.pollet@universite-paris-saclay.fr
