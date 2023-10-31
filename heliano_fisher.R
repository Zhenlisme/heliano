library(bedtoolsr)
library(parallel)

###set argument#####
args=commandArgs(T)
bedtools_path = args[1]
Genome_size = args[2]
joint_bedpath_df=args[3]
bed_opt=args[4]
MAX_CPU_num = detectCores() -1
CPU_num=as.numeric(args[5])
CPU_num = ifelse(CPU_num < MAX_CPU_num, CPU_num, MAX_CPU_num)
CPU_num = ifelse(CPU_num<1, 1, CPU_num)
PVALUE=as.numeric(args[6])

options(bedtools.path = bedtools_path)
#####################

Genome_size_df=read.csv2(Genome_size, stringsAsFactors = F, header = F, sep = '\t')
joint_df = read.csv2(joint_bedpath_df, stringsAsFactors = F, header = F, sep = '\t')
colnames(joint_df)=c('left', 'right')

cl = makeCluster(CPU_num)
clusterExport(cl, c('Genome_size_df', 'bedtools_path'), envir = .GlobalEnv)
clusterEvalQ(cl, list(library(bedtoolsr), options(bedtools.path = bedtools_path)))

significance_df = parApply(cl, joint_df, MARGIN = 1, function(x){
    left_bed = read.csv2(x[1], stringsAsFactors = F, header = F, sep = '\t')
    right_bed = read.csv2(x[2], stringsAsFactors = F, header = F, sep = '\t')
    fisher_result = bt.fisher(a=left_bed, b=right_bed, g=Genome_size_df, nonamecheck = FALSE, m = TRUE, s = TRUE)
    pvalue=fisher_result$right
    return(c(x[1], x[2], pvalue))
  })

stopCluster (cl)
significance_df=t(significance_df) 
colnames(significance_df)=c('leftpath', 'rightpath', 'pvalue')
significance_df=as.data.frame(significance_df, stringsAsFactors = F)
significance_df$pvalue=as.numeric(significance_df$pvalue)
write.table(significance_df, file = 'RawPvalue.txt', quote = F, sep="\t", col.names = F, row.names = F)

significance_df=significance_df[which(significance_df$pvalue<=PVALUE), ]

if(length(rownames(significance_df))){
  ## make joint and merge pvalue file
  cl = makeCluster(CPU_num)
  clusterExport(cl, c('bedtools_path', 'bed_opt'), envir = .GlobalEnv)
  clusterEvalQ(cl, list(library(bedtoolsr), options(bedtools.path = bedtools_path)))
  
  parApply(cl, significance_df, MARGIN = 1, function(x){
    pvalue = x[3]
    left_bed = read.csv2(x[1], stringsAsFactors = F, header = F, sep = '\t')
    right_bed = read.csv2(x[2], stringsAsFactors = F, header = F, sep = '\t')
    window_df = bt.intersect(a=left_bed, b=right_bed, nonamecheck = FALSE, s=TRUE, wo=TRUE)
    
    if(length(rownames(window_df))>0){
      Start = ifelse(window_df$V12=='+', window_df$V2, window_df$V8)
      Stop = ifelse(window_df$V12=='+', window_df$V9, window_df$V3)
      name=paste(window_df$V4, window_df$V10, sep = '-')
      Bitscore=(as.numeric(window_df$V5)+as.numeric(window_df$V11))/2
      window_df=data.frame(chrm=window_df$V1, start=Start, stop=Stop, 
                           combiname=name, count=1, strand=window_df$V12, 
                           pvalue=1, Bscore=Bitscore, stringsAsFactors=F)
      left_name = gsub('.bed', '', basename(x[1]))
      right_name = gsub('.bed', '', basename(x[2]))
      
      window_df = window_df[order(window_df$chrm, as.numeric(window_df$start)), ]
      cluster_df = bt.cluster(window_df, d=-10)
      cluster_count = length(unique(cluster_df[, 9]))
      window_df$pvalue=pvalue
      window_df$count=cluster_count
      filename=paste(bed_opt, '/',left_name, '-', right_name, '.bed', sep = '')
      write.table(window_df, file = filename, quote = F, row.names = F, col.names = F, sep = '\t')
    }
  })

  stopCluster (cl)
  }



