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
ORF_BED=args[7]
Strategy=args[8]
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

## make joint and merge pvalue file
cl = makeCluster(CPU_num)
clusterExport(cl, c('bedtools_path', 'bed_opt', 'ORF_BED'), envir = .GlobalEnv)
clusterEvalQ(cl, list(library(bedtoolsr), options(bedtools.path = bedtools_path)))

if(Strategy=='1'){
  if(length(rownames(significance_df))){
    parApply(cl, significance_df, MARGIN = 1, function(x){
      #apply(significance_df, MARGIN = 1, function(x){
      #print(x)
      pvalue = x[3]
      window_df = bt.intersect(a=x[1], b=x[2], nonamecheck = FALSE, s=TRUE, wo=TRUE)
      if(length(rownames(window_df))>0){
        Start = ifelse(window_df$V6=='+', window_df$V7, window_df$V14)
        Stop = ifelse(window_df$V6=='+', window_df$V14, window_df$V7)
        name=paste(window_df$V4, window_df$V11, sep = '-')
        Bitscore=(as.numeric(window_df$V5)+as.numeric(window_df$V12))/2
        window_df=data.frame(chrm=window_df$V1, start=Start, stop=Stop, 
                             combiname=name, count=1, strand=window_df$V6, 
                             pvalue=1, Bscore=Bitscore, stringsAsFactors=F)
        left_name = gsub('.bed', '', basename(x[1]))
        right_name = gsub('.bed', '', basename(x[2]))
        
        window_df = window_df[order(window_df$chrm, as.numeric(window_df$start)), ]
        cluster_df = bt.cluster(window_df, d=-10)
        cluster_count = length(unique(cluster_df[, 9]))
        window_df$pvalue=pvalue
        window_df$count=cluster_count
        window_df$repalt='rep'
        #print(head(window_df))
        ##################### dedup function ############################
        dedup_func = function(inpt_fisher_df){
          intersection_df = bt.intersect(inpt_fisher_df, inpt_fisher_df, F=1, s=T, wo=T)
          
          ## To select the one who cover any other lines
          intersection_df = intersection_df[which(!(intersection_df$V2==intersection_df$V11 & intersection_df$V3==intersection_df$V12)),
                                            c(1:9)]
          intersection_df = unique(intersection_df)
          non_overlapped_bed = bt.intersect(inpt_fisher_df, intersection_df, wa=T, f=1, F=1, v=T, s=T)
          
          if(nrow(intersection_df)>0){
            intersection_df$V9='alt'
          }
          final_opt = rbind(non_overlapped_bed, intersection_df)
          return(final_opt)
        }
        ################## Function end ################################
        if(file.exists(ORF_BED)){
          ORF_fisherpart = bt.intersect(window_df, ORF_BED, F=1, wa=T)
          nonauto_fisherpart = bt.intersect(window_df, ORF_BED, F=.5, wa=T, v = T)
          Auto_dedup = data.frame()
          Nonauto_dedup = data.frame()
          if(nrow(ORF_fisherpart)>0){
            Auto_dedup = dedup_func(ORF_fisherpart)
          }
          if(nrow(nonauto_fisherpart>0)){
            Nonauto_dedup = dedup_func(nonauto_fisherpart)
          }
          Dedup_df = rbind(Auto_dedup, Nonauto_dedup)
          
        }else{
          Dedup_df = dedup_func(window_df)
        }
        filename=paste(bed_opt, '/',left_name, '-', right_name, '.bed', sep = '')
        write.table(Dedup_df, file = filename, quote = F, row.names = F, col.names = F, sep = '\t')
      }
    })
    
  }
}else{
  if(length(rownames(significance_df))){
    parApply(cl, significance_df, MARGIN = 1, function(x){
      pvalue = x[3]
      window_df = bt.intersect(a=x[1], b=x[2], nonamecheck = FALSE, s=TRUE, wo=TRUE)
      if(length(rownames(window_df))>0){
        Start = ifelse(window_df$V6=='+', window_df$V7, window_df$V14)
        Stop = ifelse(window_df$V6=='+', window_df$V14, window_df$V7)
        name=paste(window_df$V4, window_df$V11, sep = '-')
        Bitscore=(as.numeric(window_df$V5)+as.numeric(window_df$V12))/2
        window_df=data.frame(chrm=window_df$V1, start=Start, stop=Stop, 
                             combiname=name, count=1, strand=window_df$V6, 
                             pvalue=1, Bscore=Bitscore, stringsAsFactors=F)
        left_name = gsub('.bed', '', basename(x[1]))
        right_name = gsub('.bed', '', basename(x[2]))
        
        window_df = window_df[order(window_df$chrm, as.numeric(window_df$start)), ]
        cluster_df = bt.cluster(window_df, d=-10)
        cluster_count = length(unique(cluster_df[, 9]))
        window_df$pvalue=pvalue
        window_df$count=cluster_count
        window_df$repalt='rep'
        filename=paste(bed_opt, '/',left_name, '-', right_name, '.bed', sep = '')
        write.table(window_df, file = filename, quote = F, row.names = F, col.names = F, sep = '\t')
      }
    })
    
  }
}


stopCluster (cl)

## To filter out pairs whose sequences are similar.
fisher_bedfilelist = list.files(bed_opt)

cl = makeCluster(CPU_num)
clusterExport(cl, c('bedtools_path', 'bed_opt'), envir = .GlobalEnv)
clusterEvalQ(cl, list(library(bedtoolsr), options(bedtools.path = bedtools_path)))

parSapply(cl, fisher_bedfilelist, function(x){
  fisherbed_file = paste(bed_opt, '/', x, sep = '')
  split_list = strsplit(gsub('.bed', '', x),split = '-')
  leftname = split_list[[1]][1]
  rightname = split_list[[1]][2]
  classname = strsplit(leftname, split = '_left_')[[1]][1]
  left_blasnbed = paste('./SubBlastnBed/', classname, '_left/', leftname, '.bed', sep = '')
  right_blasnbed = paste('./SubBlastnBed/', classname, '_right/', rightname, '.bed', sep = '')
  left_bed = read.csv2(left_blasnbed, stringsAsFactors = F, header = F, sep = '\t')
  intersection_df = bt.intersect(a=left_blasnbed, b=right_blasnbed, nonamecheck = FALSE, s=FALSE, f=0.8)
  proportion = nrow(unique(intersection_df))/nrow(left_bed)
  if(proportion>=0.3){
    file.remove(fisherbed_file)
  }
})
stopCluster (cl)
