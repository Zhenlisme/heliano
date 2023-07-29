# HElA: A rapid and accurate tool for detection of Helitron-like elements.
Helitron-like elements including Helitron, Helentron and Helitron2 are one of the class 2 transposons. They have been found in a diverse range of species and seem to play major roles in the evolution of host genomes. Although they have benn well known for more than twenty years, they were widely admitted as one of few remaining transposons difficult to identify. Here, we propose HELA (Helitron-like elements annotator) which is a rapid and accurate tool for detection of Helitron-like elements.
# Dependencies:
```
- python =3.9.0
- r-base =4.1
- biopython
- pybedtools =0.9.0=py39hd65a603_2
- r-bedtoolsr
- r-seqinr =4.2_16=r41h06615bd_0
- bedtools =2.30.0
- dialign2 =2.2.1
- mafft
- cd-hit =4.8.1
- blast =2.2.31
- emboss =6.6.0
- hmmer =3.3.2
- genometools-genometools =1.6.2=py39h58cc16e_6
```
# Install
## conda
```
#creat enviroment
conda create -n hela
#activate enviroment
conda activate hela
# install 
conda install -c zhenlisme hela -c bioconda -c conda-forge
conda deactivate
```
## mamba
```
#creat enviroment
mamba create -n hela
#activate enviroment
mamba activate hela
# install 
mamba install -c zhenlisme hela -c bioconda -c conda-forge
mamba deactivate
```
# Usage
```
HELA -h
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

# To contact us
For any questions, please email us:
### Zhen Li
zhenliseme@gmail.com or li.zhen3@universite-paris-saclay.fr
### Nicolas Pollet
nicolas.pollet@universite-paris-saclay.fr
