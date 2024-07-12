[![Anaconda Version](https://anaconda.org/zhenlisme/heliano/badges/version.svg)](https://anaconda.org/zhenlisme/heliano/)
[![Anaconda Version](https://anaconda.org/zhenlisme/heliano/badges/latest_release_date.svg)](https://anaconda.org/zhenlisme/heliano/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](./LICENSE)

# HELIANO: A fast and accurate tool for detection of Helitron-like elements.
Helitron-like elements (HLE1 and HLE2) are DNA transposons. They have been found in diverse species and seem to play significant roles in the evolution of host genomes. Although known for over twenty years, Helitron sequences are still challenging to identify. Here, we propose HELIANO (Helitron-like elements annotator) as an efficient solution for detecting Helitron-like elements.

# Table of contents
- [Update Note](#Update-Note)
- [Dependencies](#dependencies)
- [Installation](#installation)
  * [conda](#conda)
  * [mamba](#mamba)
  * [Manual installation](#manual-installation)
- [Usage](#usage)
     - [Test run](#Perform-a-test-run-of-HELIANO)
     - [Making Consensus](#Generation-for-consensus-sequences)
     - [Dis-denovo prediction](#Dis-denovo-prediction)
- [References](#References)
- [FAQ](#Frequently-asked-questions)
- [Release history](#Release-history)
- [To contact us](#to-contact-us)

# Update Note:
1) Since version 1.1.0, HELIANO will use the term HLE1 to refer to the canonical Helitron (called Helitron in v1.0.2) and the term HLE2 to refer to the non-canonical Helitrons (called HLE2 in v1.0.2). 
See figure below:

<img src="https://github.com/Zhenlisme/heliano/assets/54898847/32b1ba39-3d4c-428f-8dda-fdde63c50003" width="660" height="420">

2) From version 1.1.0, users are allowed to input a pair file as a complementary for LTS-RTS pair information. This will help a lot to search for HLEs in close species. More information see [here](#Dis-denovo-prediction).
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
If mamba is not installed on your system, you can install it with the following commands easily.
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
mamba install zhenlisme::HELIANO -c conda-forge -c bioconda
mamba deactivate
```
## conda
```
#create the HELIANO environment
conda create -n HELIANO
#activate the HELIANO environment
conda activate HELIANO
# installation 
conda install zhenlisme::HELIANO -c conda-forge -c bioconda
conda deactivate
```
## manual installation
Before installation , you need to be sure that all dependencies have been installed in your computer and that their path are defined into your environmental variables. All dependencies could be installed via conda/mamba.  
1. download the latest HELIANO package.  
`git clone https://github.com/Zhenlisme/HELIANO.git`
2. switch to the source code dorectory that you cloned at the last step.  
   `cd HELIANO/`  
3. run configure file.  
   `bash configure.sh`
4. You can find HELIANO in the bin directory.
# Usage
### Activate the HELIANO conda environment (for conda/mamba installation)  
`conda activate HELIANO`
### Check the HELIANO binary  
`heliano -h`
```
usage: heliano [-h] -g GENOME [-w WINDOW] [-dm DISTANCE_DOMAIN] [-pt {0,1}] [-is1 {0,1}] [-is2 {0,1}] [-sim_tir {100,90,80}] [-p PVALUE]
               [-s SCORE] [-ts TERMINAL_SEQUENCE] [--dis_denovo] [--multi_ts] -o OPDIR [-n PROCESS] [-v]

heliano can detect and classify different variants of Helitron-like elements: HLE1 and HLE2. Please visit
https://github.com/Zhenlisme/heliano/ for more information. Email us: zhen.li3@universite-paris-saclay.fr

options:
  -h, --help            show this help message and exit
  -g GENOME, --genome GENOME
                        The genome file in fasta format.
  -w WINDOW, --window WINDOW
                        To check terminal signals within a given window bp upstream and downstream of ORF ends, default is 10 kb.
  -dm DISTANCE_DOMAIN, --distance_domain DISTANCE_DOMAIN
                        The distance between HUH and Helicase domain, default is 2500.
  -pt {0,1}, --pair_helitron {0,1}
                        For HLE1, its 5' and 3' terminal signal pairs should come from the same autonomous helitorn or not. 0: no, 1:
                        yes. default yes.
  -is1 {0,1}, --IS1 {0,1}
                        Set the insertion site of autonomous HLE1 as A and T. 0: no, 1: yes. default yes.
  -is2 {0,1}, --IS2 {0,1}
                        Set the insertion site of autonomous HLE2 as T and T. 0: no, 1: yes. default yes.
  -sim_tir {100,90,80}, --simtir {100,90,80}
                        Set the simarity between short inverted repeats(TIRs) of HLE2. Default 100.
  -p PVALUE, --pvalue PVALUE
                        The p-value for fisher's exact test. default is 1e-5.
  -s SCORE, --score SCORE
                        The minimum bitscore of blastn for searching for homologous sequences of terminal signals. From 30 to 55,
                        default is 32.
  -ts TERMINAL_SEQUENCE, --terminal_sequence TERMINAL_SEQUENCE
                        The terminal sequence file. You can find it in the output of previous run (named as pairlist.tbl).
  --dis_denovo          If you use this parameter, you refuse to search for LTS/RTS de novo, instead you will only use the LTS/RTS
                        information described in the terminal sequence file.
  --multi_ts            To allow an auto HLE to have multiple terminal sequences. If you enable this, you might find nonauto HLEs coming
                        from the same auto HLE have different terminal sequences.
  -o OPDIR, --opdir OPDIR
                        The output directory.
  -n PROCESS, --process PROCESS
                        Maximum of threads to be used.
  -v, --version         show program's version number and exit
```
### Perform a test run of HELIANO  
##### Here we will use the chromosome 18 of Fusarium oxysporum strain Fo5176 as an example, where you can find it in file test.fa .
Perform the following code:
`heliano -g test.fa -is1 0 -is2 0 -o test_opt -w 15000`
### HELIANO outputs
You will find two main result files when HELIANO program runs successfully.
1. RC.representative.bed: the predicted HLE1/HLE2 coordinates in bed format (available in the file test.opt.tbl in this repository).
2. RC.representative.fa: the predicted HLE1/HLE2 sequences in fasta format.
3. pairlist.tbl: The file for LTS-RTS pair information.
Other files or directories are intermediate outputs.
1. TIR_count.tbl: Table for counts of terminal inverted repeats of each HLE subfamily.
2. Boundary.tbl: Table for the conservation of flanking regions of each HLE subfamily.
3. HLE1/ or HLE2/: Directory for intermediate files when detecting HLE1/HLE2.
##### Explanation for RC.representative.bed
There are 11 columns in RC.representative.bed file:  
|chrm-id|start|end|subfamily|occurence|strand|pvalue|TS_blastn_identity|variant|type|name|
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
|CP128282.1|53617|59046|HLE2_left_18-HLE2_right_18|7|-|6.3390e-07|60|HLE2|auto|insertion_HLE2_auto_1|
|CP128282.1|83425|88824|HLE2_left_18-HLE2_right_18|7|-|6.3390e-07|60|HLE2|auto|insertion_HLE2_auto_2|
|CP128282.1|94525|99924|HLE2_left_18-HLE2_right_18|7|+|6.3390e-07|60|HLE2|auto|insertion_HLE2_auto_3|
|CP128282.1|306838|312276|HLE2_left_18-HLE2_right_18|7|+|6.3390e-07|60|HLE2|auto|insertion_HLE2_auto_4|
|CP128282.1|665681|668554|HLE2_left_33-HLE2_right_32|3|-|1.3547e-05|60|HLE2|nonauto|insertion_HLE2_nonauto_1|
|CP128282.1|668719|671599|HLE2_left_33-HLE2_right_32|3|+|1.3547e-05|60|HLE2|nonauto|insertion_HLE2_nonauto_2|
|CP128282.1|855863|858738|HLE2_left_33-HLE2_right_32|3|+|1.3547e-05|60|HLE2|nonauto|insertion_HLE2_nonauto_3|
|CP128282.1|855863|880806|HLE2_left_33-HLE2_right_32|3|+|1.3547e-05|60|HLE2|auto|insertion_HLE2_auto_5|
|CP128282.1|926221|928267|CP128282.1-926221-928267|1|-|1|0|HLE2|orfonly|insertion_HLE2_orfonly_1|
|CP128282.1|963556|968985|HLE2_left_18-HLE2_right_18|7|+|6.3390e-07|60|HLE2|auto|insertion_HLE2_auto_6|
|CP128282.1|1107206|1112635|HLE2_left_18-HLE2_right_18|7|+|6.3390e-07|60|HLE2|auto|insertion_HLE2_auto_7|
|CP128282.1|1259412|1264991|HLE2_left_18-HLE2_right_18|7|-|6.3390e-07|60|HLE2|auto|insertion_HLE2_auto_8|
##### Detailed explaination for each column.
Notice: The insertions that encode Rep/helicase are considered putative autonomous HLEs.
|Columns|Explaination|
| ---- | ---- |
|chrm-id|chromosome id|
|start|start site of HLE|
|stop|stop site of HLE|
|subfamily|heliano classification|
|occurence|how often this subfamily occurred in genome|
|strand|the insertion is on which strand|
|pvalue|pvalue of fisher's exact test, indicating the significance of the prediction. The lower, the more significant.|
|TS_blastn_identity|the average identity of RTS and LTS to their representative counterparts|
|variant|the insertion is HLE1 or HLE2|
|type|the mobility of HLE, either autonomous (auto) or nonautonomous (nonauto)|
|name|unique identifier for each insertion|
### Generation for consensus sequences
The HELIANO package also provides a program (heliano_cons) for generating consensus sequences of HLE.  
Check the usage of heliano_cons:  
`heliano_cons -h`
```
usage: heliano_cons [-h] -g GENOME -r REPSENBED -o OPDIR [-n PROCESS] [-v]

Making consensus for Helitron-like sequences. Please visit https://github.com/Zhenlisme/heliano/ for more information. Email us: zhen.li3@universite-paris-saclay.fr

optional arguments:
  -h, --help            show this help message and exit
  -g GENOME, --genome GENOME
                        The genome file in fasta format.
  -r REPSENBED, --repsenbed REPSENBED
                        The representative bed file (RC.representative.bed).
  -o OPDIR, --opdir OPDIR
                        The output directory.
  -n PROCESS, --process PROCESS
                        Maximum of threads to be used.
  -v, --version         show program's version number and exit
```
### Dis-denovo prediction
Since version 1.1.0, HELIANO enables prediction of HLEs with the help of pre-identified LTS-RTS pair file.
The `pairlist.tbl` can be either obtained from main directory of your previous run or user-defined.

You can skip denovo prediction of LTS-RTS pair process (will save a lot of time),
```
heliano -g test.fa -is1 0 -is2 0 -o test_opt -w 15000 --ts pairlist.tbl --dis_denovo
```
or not skip the denovo prediction of LTS-RTS process
```
heliano -g test.fa -is1 0 -is2 0 -o test_opt -w 15000 --ts pairlist.tbl
```
# References
### If you find HELIANO useful to you, please cite:
Li Z , Gilbert C , Peng H , Pollet N. "Systematic annotation of Helitron-like elements in eukaryote genomes using HELIANO." bioRxiv (2024): 2024-02. [doi: 10.1101/2024.02.08.579435](https://doi.org/10.1101/2024.02.08.579435)  

Li Z , Pollet N. "HELIANO: a Helitron-like element annotator." Zenodo (2024). [doi: 10.5281/zenodo.10625239](https://doi.org/10.5281/zenodo.10625239)

# Frequently asked questions
### 1. How to get fragmented copies of HLEs?
HELIANO is designed to predict complete insertions of Helitron-like elements (HLE), with the limitation that fragmented insertions will not be reported. To identify fragmented insertions, we recommend running RepeatMasker or BLASTN using HELIANO predictions as the query. Before you run RepeatMasker or BLASTN, we suggest mask the HLE query with a trusted non-HLE TE database because other non-HLE TEs might insert into long HLEs which would inflates sequence length and result in misannotation.
### 2. How to choose parameters properly?
For a precise and quick search, you can use the strigent parameter '-pt 1 -is1 1 -is2 1 -p 1e-5 -s 30 -pt 1 -sim_tir 100' that considered the preferred insertion sites of HLE. For big or complex genomes (e.g., maize genome), I just recommed you use the strigent parameter set. But not all HLEs obey their regular preferring insertion sites. If you want to explore more in your interested genome, you can use the loose parameter set, e.g., '-is1 0 -is2 0 -sim_tir 90', and you will have more predictions and longer execution time. Note that the parameters 'is2' and '-sim_tir' are only for HLE2s, and 'is1' and '-pt' are only for Helitrons.
# Release history
### v1.0.1
Initial version
### v1.1.0
1. Replace term Helitron as HLE1 and Helentron as HLE2.
2. Enable to predict HLEs based on a pre-identified LTS-RTS pair file. (see -ts and --dis_denovo parameters)
3. Add a new parameter that allows an auto HLE to have multiple terminal sequences. (see '--multi_ts' parameter)
### v1.2.0
1. Add parameter '--nearest' that allow users to find terminal pairs whose LTS and RTS are closest with each other. By default, HELIANO will try to find the furthest pairs.
2. Add parameter '-dn' that allow user to define the length of nonautonomous HLEs. By default (dn 0), HELIANO will deduce it automatically.
# To contact us
For any questions, please open an issue in the 'issues' section or email us: zhen.li3@universite-paris-saclay.fr or nicolas.pollet@universite-paris-saclay.fr
