usage: HELA [-h] -g GENOME [-w WINDOW] [-dm DISTANCE_DOMAIN] [-pt {0,1}] [-is1 {0,1}] [-is2 {0,1}]
                  [-sim_tir {100,90,80}] [-p PVALUE] -o OPDIR [-n PROCESS]

Helitron detection.

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
