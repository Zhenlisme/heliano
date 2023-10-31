## This script was used to split big bed file into small bed files.
library(bedtoolsr)
library(parallel)

###set argument
options(scipen = 200)

args=commandArgs(T)
left_beddir=args[1]
right_beddir=args[2]
combinid_file=args[3]
subed_dir=args[4]
Genome_size_file = args[5]
HALF_DIST = as.numeric(args[6])
bedtools_path = args[7]
# Set bedtools option
options(bedtools.path = bedtools_path)

## multiple processing
MAX_CPU_num = detectCores() -1
CPU_num=as.numeric(args[8])
CPU_num = ifelse(CPU_num < MAX_CPU_num, CPU_num, MAX_CPU_num)
CPU_num = ifelse(CPU_num<1, 1, CPU_num)

## To make directory
if(dir.exists(subed_dir)){
  unlink(subed_dir,recursive = TRUE)
}

left_flank_dir=paste(subed_dir, 'leftflank', sep = '/')
dir.create(left_flank_dir, recursive = TRUE)
right_flank_dir=paste(subed_dir, 'rightflank', sep = '/')
dir.create(right_flank_dir, recursive = TRUE)

## read genome size
genome_size_df = read.csv2(Genome_size_file, sep = '\t', header = F, stringsAsFactors = F)
colnames(genome_size_df)=c('chrmid', 'length')

#### Output left bed file ####
left_bedfile_list = list.files(left_beddir, full.names = T)

cl = makeCluster(CPU_num)
clusterExport(cl, c('left_flank_dir', 'genome_size_df', 'bedtools_path', 'HALF_DIST'), envir = .GlobalEnv)
clusterEvalQ(cl, list(library(bedtoolsr), options(bedtools.path = bedtools_path, scipen = 200)))

parSapply(cl, left_bedfile_list, function(x){
  sub_leftbed = read.csv2(x, stringsAsFactors = F, header = F, sep = '\t')
  sub_leftbed$V2=as.numeric(sub_leftbed$V2)
  sub_leftbed$V3=as.numeric(sub_leftbed$V3)
  
  ## To merge tandem repeats
  sub_leftbed = sub_leftbed[order(sub_leftbed$V1, sub_leftbed$V2), ]
  sub_leftbed = bt.merge(sub_leftbed, s=TRUE, d = 100, c="4,5,6", o="first,max,first")
  
  ## To extend bed files
  sub_leftbed$V3=ifelse(sub_leftbed$V6=='+', sub_leftbed$V2+HALF_DIST, sub_leftbed$V3)
  sub_leftbed$V2=ifelse(sub_leftbed$V6=='+', sub_leftbed$V2, sub_leftbed$V3-HALF_DIST)
  sub_leftbed$V2=ifelse(sub_leftbed$V2<=0, 1, sub_leftbed$V2)
  colnames(sub_leftbed)=c('chrmid', 'start', 'stop', 'name', 'score', 'strand')
  sub_leftbed = merge(sub_leftbed, genome_size_df)
  sub_leftbed$V3 = ifelse(sub_leftbed$stop <= sub_leftbed$length, sub_leftbed$stop, sub_leftbed$length)
  sub_leftbed = sub_leftbed[, c('chrmid', 'start', 'stop', 'name', 'score', 'strand')]
  ## To output
  bnm=basename(x)
  sub_leftfile = paste(left_flank_dir, '/', bnm, sep = '')
  bt.sort(sub_leftbed, output = sub_leftfile)
})
stopCluster (cl)

#### Output right bed file####
right_bedfile_list = list.files(right_beddir, full.names = T)
## multiple processing
cl = makeCluster(CPU_num)
clusterExport(cl, c('right_flank_dir', 'genome_size_df', 'bedtools_path', 'HALF_DIST'), envir = .GlobalEnv)
clusterEvalQ(cl, list(library(bedtoolsr), options(bedtools.path = bedtools_path, scipen = 200)))

parSapply(cl, right_bedfile_list, function(x){
  sub_rightbed = read.csv2(x, stringsAsFactors = F, header = F, sep = '\t')
  sub_rightbed$V2=as.numeric(sub_rightbed$V2)
  sub_rightbed$V3=as.numeric(sub_rightbed$V3)
  
  ## To merge tandem repeats
  sub_rightbed = sub_rightbed[order(sub_rightbed$V1, sub_rightbed$V2), ]
  sub_rightbed = bt.merge(sub_rightbed, s=TRUE, d = 100, c="4,5,6", o="first,max,first")

  ## To extend bed files
  sub_rightbed$V3=ifelse(sub_rightbed$V6=='+', sub_rightbed$V3, sub_rightbed$V2+HALF_DIST)
  sub_rightbed$V2=ifelse(sub_rightbed$V6=='+', sub_rightbed$V3-HALF_DIST, sub_rightbed$V2)
  sub_rightbed$V2=ifelse(sub_rightbed$V2<=0, 1, sub_rightbed$V2)
  colnames(sub_rightbed)=c('chrmid', 'start', 'stop', 'name', 'score', 'strand')
  sub_rightbed = merge(sub_rightbed, genome_size_df, by='chrmid')
  sub_rightbed$stop = ifelse(sub_rightbed$stop <= sub_rightbed$length, sub_rightbed$stop, sub_rightbed$length)
  
  ## To output 
  sub_rightbed = sub_rightbed[, c('chrmid', 'start', 'stop', 'name', 'score', 'strand')]
  bnm=basename(x)
  sub_rightfile = paste(right_flank_dir, '/', bnm, sep = '')
  bt.sort(sub_rightbed, output = sub_rightfile)
})
stopCluster (cl)
## To create file path file

## read combined id file
leftname_list = sapply(left_bedfile_list, function(x){strsplit(basename(x), split = '\\.')[[1]][1]})
rightname_list = sapply(right_bedfile_list, function(x){strsplit(basename(x), split = '\\.')[[1]][1]})

combine_id=read.csv2(combinid_file, stringsAsFactors = F, header = F, sep = '\t')
colnames(combine_id)=c('leftname', 'rightname')
combine_id=combine_id[which(combine_id$leftname %in% leftname_list & combine_id$rightname %in% rightname_list), ]

combine_id$leftname=paste(subed_dir, '/leftflank/', combine_id$leftname, '.bed',sep = '')
combine_id$rightname=paste(subed_dir, '/rightflank/', combine_id$rightname, '.bed', sep = '')
write.table(combine_id, file = paste(subed_dir, '/', 'left_right.path.join', sep = ''), 
            quote = F, sep="\t", col.names = F, row.names = F)
