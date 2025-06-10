library(bedtoolsr)
library(parallel)

###set argument#####
args=commandArgs(T)
bedtools_path = args[1]
bed_fisher_dir=args[2]
boundary_dir = args[3]
Genome_path=args[4]
Genome_size = args[5]
MAX_CPU_num = detectCores() -1
CPU_num=as.numeric(args[6])
CPU_num = ifelse(CPU_num < MAX_CPU_num, CPU_num, MAX_CPU_num)
CPU_num = ifelse(CPU_num<1, 1, CPU_num)
options(bedtools.path = bedtools_path)

##################### To output sequences ####################
Genome_size_df=read.csv2(Genome_size, stringsAsFactors = F, header = F, sep = '\t')
colnames(Genome_size_df)=c('chr', 'length')
fisher_file_list = list.files(bed_fisher_dir, full.names = T, pattern = '.bed')

cl = makeCluster(CPU_num)
clusterExport(cl, c('bedtools_path', 'Genome_path', 'Genome_size_df'), envir = .GlobalEnv)
clusterEvalQ(cl, list(library(bedtoolsr), options(bedtools.path = bedtools_path)))

parSapply(cl, fisher_file_list, function(x){
  pairname=gsub('\\.bed', '', basename(x))
  candidate_bed=bt.sort(x)
  merged_candidate = bt.merge(d=100, c = '4,5,6', o='first,mean,first')
  merged_candidate = merged_candidate[order(merged_candidate$V5, decreacing=T), ]
  row_count = nrow(merged_candidate)
  if(row_count>=2){
    row_count = ifelse(row_count>20, 20, row_count)
    merged_candidate = merged_candidate[1:row_count, c(1:6)]
    colnames(merged_candidate)=c('chr', 'start', 'stop', 'family', 'copy', 'strand')
    ### To output 100 bp terminal flanking
    terminal_flank = merged_candidate
    terminal_flank$start=terminal_flank$start-101
    terminal_flank$stop=terminal_flank$stop+99
    terminal_flank = merge(terminal_flank, Genome_size_df)
    terminal_flank$start = ifelse(terminal_flank$start<0, 0, terminal_flank$start)
    terminal_flank$stop = ifelse(terminal_flank$stop>terminal_flank$length, terminal_flank$length, terminal_flank$stop)
    terminal_flank = terminal_flank[which(terminal_flank$stop - terminal_flank$start >200), ]
  }
  
})

stopCluster (cl)