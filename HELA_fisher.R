library(bedtoolsr)
library(parallel)
###set argument#####
args=commandArgs(T)
bedtools_path = args[1]
Genome_size = args[2]
joint_bedpath_df=args[3]
bed_opt=args[4]
MAX_CPU_num = detectCores()
CPU_num=as.numeric(args[5])
CPU_num = ifelse(CPU_num < MAX_CPU_num, CPU_num, MAX_CPU_num)
PVALUE=as.numeric(args[6])
JOINT_BED=args[7]

#setwd('~/remote2/Helitron/Xenopus/XL_test_1e-3/')
#bedtools_path = '/home/zhenli/bioinfo/bedtools2/bin/'
#Genome_size = 'Genome.size'
#joint_bedpath_df = 'Helitron_SubBed/left_right.path.join'
#bed_opt='bed_opt.tbl'
#CPU_num = 4
#PVALUE=1e-5
#JOINT_BED='Helitron.joint.filtered.bed'
## bedtools path
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
    fisher_result = bt.fisher(a=left_bed, b=right_bed, g=Genome_size_df, m = TRUE, s = TRUE)
    pvalue=fisher_result$right
    left_name = gsub('.bed', '', basename(x[1]))
    right_name = gsub('.bed', '', basename(x[2]))
    name_marker=paste(left_name, right_name, sep='-')
    return(c(name_marker, pvalue))
  })

stopCluster (cl)
significance_df=t(significance_df) 
colnames(significance_df)=c('combiname', 'pvalue')
significance_df=as.data.frame(significance_df, stringsAsFactors = F)
significance_df$pvalue=as.numeric(significance_df$pvalue)
significance_df=significance_df[which(significance_df$pvalue<=PVALUE), ]

## read JOINT_BED
joint_bed = read.csv2(JOINT_BED, stringsAsFactors = F, header = F, sep = '\t')
colnames(joint_bed)=c('chrm', 'start', 'stop', 'combiname', 'count', 'strand', 'fake_p', 'Bscore')

significance_df = merge(joint_bed, significance_df)
significance_df = significance_df[,c(2,3,4,1,5,6,9,8)]


##caculate cluster count
significance_df=significance_df[order(significance_df$chrm, significance_df$start),]
significance_df=bt.cluster(significance_df, d=-10)
colnames(significance_df)=c('chrm', 'start', 'stop', 'combiname', 'fake_count', 'strand', 'pvalue', 'Bscore', 'cluster')
significance_df$combiname=as.character(significance_df$combiname)

combiname_list=unique(significance_df$combiname)

Count_list = lapply(combiname_list, function(x){
  count=length(unique(significance_df[which(significance_df$combiname==x), 9]))
  return(c(x, count))
})
Count_df = as.data.frame(do.call('rbind', Count_list), stringsAsFactors = F)
colnames(Count_df)=c('combiname', 'count')

significance_df=merge(significance_df, Count_df)

significance_df=significance_df[,c(2, 3, 4, 1, 10, 6, 7, 8)]
significance_df=significance_df[order(significance_df$chrm, significance_df$start),]
write.table(significance_df, file = bed_opt, quote = F, sep="\t", col.names = F, row.names = F)

pairname_info_df=unique(significance_df[,c(4,5,7)])
write.table(pairname_info_df, file = "../Pairname.siginfo.tbl", quote = F, sep="\t", col.names = F, row.names = F, append=TRUE)

