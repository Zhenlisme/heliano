# HELIANO: A fast and accurate tool for detection of Helitron-like elements.
Helitron-like elements (Helitron, Helentron and Helitron2) belong the class 2 transposons. They have been found in diverse species and seem to play significant roles in the evolution of host genomes. Although known for over twenty years, Helitron sequences are still challenging to identify. Here, we propose HELIANO (Helitron-like elements annotator) as an efficient solution for detecting Helitron-like elements.

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
- rnabob = 2.2.1
```
# Installation
## mamba (Recommendation)
If you don't have a mamba, you can install it with the following commands easily.
```
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
bash Mambaforge-$(uname)-$(uname -m).sh -b
```
Then you can install HELIANO with mamba.
```
#create the HELIANO environment
mamba create -n HELIANO
#activate the HELIANO environment
mamba activate HELIANO
# install 
mamba install zhenlisme::HELIANO -c bioconda -c conda-forge
mamba deactivate
```
## conda
```
#create the HELIANO environment
conda create -n HELIANO
#activate the HELIANO environment
conda activate HELIANO
# installation 
conda install zhenlisme::HELIANO -c bioconda -c conda-forge
conda deactivate
```
## manual installation
Before installation , you need to be sure that all dependencies have been installed in your computer and their pathes have been added into your environmental variables. All dependencies could be installed via conda/mamba.  
1. download the latest HELIANO package.  
`git clone https://github.com/Zhenlisme/HELIANO.git`
2. switch to the source code dorectory that you cloned at the last step.  
   `cd HELIANO/`  
3. run configure file.  
   `bash configure.sh`
4. You can find HELIANO in bin/ directory.
# Usage
### 1. Activate the HELIANO conda environment (for conda/mamba installation)  
`conda activate HELIANO`  
### 2. Check the HELIANO binary  
`HELIANO -h`
```
usage: HELIANO [-h] -g GENOME [-w WINDOW] [-dm DISTANCE_DOMAIN] [-pt {0,1}] [-is1 {0,1}] [-is2 {0,1}] [-sim_tir {100,90,80}] [-p PVALUE] -o OPDIR [-n PROCESS] [-v]

HELIANO can detect and classify different variants of Helitron-like elements: Helitron and Helentrons. Please visit https://github.com/Zhenlisme/HELIANO/ for more information. Email us:
zhen.li3@universite-paris-saclay.fr

optional arguments:
  -h, --help            show this help message and exit
  -g GENOME, --genome GENOME
                        The genome file in fasta format.
  -w WINDOW, --window WINDOW
                        To check terminal signals within a given window bp upstream and downstream of ORF ends, default is 10 kb.
  -dm DISTANCE_DOMAIN, --distance_domain DISTANCE_DOMAIN
                        The distance between HUH and Helicase domain, default is 2500.
  -pt {0,1}, --pair_helitron {0,1}
                        For Helitron, its 5' and 3' terminal signal pairs should come from the same autonomous helitorn or not. 0: no, 1: yes. default no.
  -is1 {0,1}, --IS1 {0,1}
                        Set the insertion site of autonomous Helitron as A and T. 0: no, 1: yes. default no.
  -is2 {0,1}, --IS2 {0,1}
                        Set the insertion site of autonomous Helentron/Helitron2 as T and T. 0: no, 1: yes. default no.
  -sim_tir {100,90,80}, --simtir {100,90,80}
                        Set the simarity between short inverted repeats(TIRs) of Helitron2/Helentron. Default 100.
  -p PVALUE, --pvalue PVALUE
                        The p-value for fisher's exact test.
  -o OPDIR, --opdir OPDIR
                        The output directory.
  -n PROCESS, --process PROCESS
                        Maximum of threads to be used.
  -v, --version         show program's version number and exit
```
### 3. Perform a test run of HELIANO  
##### Here we will use the genome of Xenopus tropicalis as an example
You can download it here: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/195/GCF_000004195.4_UCB_Xtro_10.0/GCF_000004195.4_UCB_Xtro_10.0_genomic.fna.gz.

Perform the following code:

`heliano -g genome.fa -pt 1 -is1 1 -is2 1 -sim_tir 90 -p 1e-5 -n 30 -o HELIANO_opt`
### 4. HELIANO outputs
You will find two main files if program runs successfully.
1. RC.representative.bed: the predicted Helitron/Helentron coordinations.
2. RC.representative.fa: the predicted Helitron/Helentron sequences in fasta format.
Notice: the two files will be empty if no Helitron/Helentron detected by HELIANO program.
##### Explanation for RC.representative.bed
There are 11 columns in RC.representative.bed file:  
|chrm-id|start|end|subfamily|occurence|strand|pvalue|TS_blastn_identity|variant|type|name|
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
|NC_030677.2|20246175|20247385|Helentron_left_73-Helentron_right_74|79|+|2.0432e-176|104|Helentron|nonauto|insertion_Helentron_nonauto_1|
|NC_030677.2|23640161|23640743|Helentron_left_73-Helentron_right_74|79|+|2.0432e-176|74.3|Helentron|nonauto|insertion_Helentron_nonauto_2|
......
|chrm-id|chromosome id|
|start|start site of HLE|
|stop|stop site of HLE|
|subfamily|heliano classification|
|occurence|how often this subfamily occurred in genome|
|strand|strand|
|pvalue|pvalue of fisher's exact test, indicating the significance of the prediction. The lower, the more significant.|
|TS_blastn_identity|the average identity of RTS and LTS to their representative counterparts|
|variant|Helitron/Helentron|
|type|indicate the mobility of HLE, either autonomous or nonautonomous|
|name|unique identifier for each insertion|

# To contact us
For any questions, please email us: zhen.li3@universite-paris-saclay.fr or nicolas.pollet@universite-paris-saclay.fr
