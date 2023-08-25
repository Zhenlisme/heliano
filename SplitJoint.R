## This script was used to split big bed file into small bed files.
library(parallel)
library(bedtoolsr)

###set argument
options(scipen = 200)

args=commandArgs(T)
left_bedfile=args[1]
right_bedfile=args[2]
combinid_file=args[3]
subed_dir=args[4]
MAX_CPU_num = detectCores() -1
CPU_num=as.numeric(args[5])
CPU_num = ifelse(CPU_num < MAX_CPU_num, CPU_num, MAX_CPU_num)
CPU_num = ifelse(CPU_num<1, 1, CPU_num)
Genome_size_file = args[6]
HALF_DIST = as.numeric(args[7])
bedtools_path = args[8]

# Set bedtools option
options(bedtools.path = bedtools_path)

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

## read combined id file
if(file.info(combinid_file)$size != 0){
  combine_id=read.csv2(combinid_file, stringsAsFactors = F, header = F, sep = '\t')
  colnames(combine_id)=c('leftname', 'rightname')
}

#### Output left bed file ####

## loading left bed file####
left_bed=read.csv2(left_bedfile, stringsAsFactors = F, header = F, sep = '\t')
colnames(left_bed)=c('chrmid', 'start', 'stop', 'name', 'score', 'strand')

if(file.info(combinid_file)$size == 0){
  leftname_list = unique(left_bed$name)
}else{
  leftname_list = unique(combine_id$leftname)
  left_drop_list = setdiff(combine_id$leftname, unique(left_bed$name))
  leftname_list = leftname_list[which(!leftname_list %in% left_drop_list)]
}

cl = makeCluster(CPU_num)
clusterExport(cl, c('left_bed', 'left_flank_dir', 'genome_size_df', 'bedtools_path', 'HALF_DIST'), envir = .GlobalEnv)
clusterEvalQ(cl, list(library(bedtoolsr), options(bedtools.path = bedtools_path, scipen = 200)))

parSapply(cl, leftname_list, function(x){
  sub_leftbed = left_bed[which(left_bed$name==x), ]
  
  ## To merge tandem repeats
  sub_leftbed = sub_leftbed[order(sub_leftbed$chrmid, sub_leftbed$start), ]
  sub_leftbed = bt.merge(sub_leftbed, s=TRUE, d = 100, c="4,5,6", o="first,max,first")
  colnames(sub_leftbed)=c('chrmid', 'start', 'stop', 'name', 'score', 'strand')
  
  ## To extend bed files
  sub_leftbed$start=as.numeric(sub_leftbed$start)
  sub_leftbed$stop=as.numeric(sub_leftbed$stop)
  sub_leftbed$stop=ifelse(sub_leftbed$strand=='+', sub_leftbed$start+HALF_DIST, sub_leftbed$stop)
  sub_leftbed$start=ifelse(sub_leftbed$strand=='+', sub_leftbed$start, sub_leftbed$stop-HALF_DIST)
  sub_leftbed$start=ifelse(sub_leftbed$start<=0, 1, sub_leftbed$start)
  sub_leftbed = merge(sub_leftbed, genome_size_df, by='chrmid')
  sub_leftbed$stop = ifelse(sub_leftbed$stop <= sub_leftbed$length, sub_leftbed$stop, sub_leftbed$length)
  
  ## To output 
  sub_leftbed = sub_leftbed[, c('chrmid', 'start', 'stop', 'name', 'score', 'strand')]
  sub_leftfile = paste(left_flank_dir, '/', x, '.bed', sep = '')
  bt.sort(sub_leftbed, output = sub_leftfile)
})
stopCluster (cl)

rm(left_bed)
gc()
#### Output right bed file####

## loading right bed file
right_bed=read.csv2(right_bedfile, stringsAsFactors = F, header = F, sep = '\t')
colnames(right_bed)=c('chrmid', 'start', 'stop', 'name', 'score', 'strand')

if(file.info(combinid_file)$size == 0){
  rightname_list = unique(right_bed$name)
}else{
  rightname_list = unique(combine_id$rightname)
  right_drop_list = setdiff(combine_id$rightname, unique(right_bed$name))
  rightname_list = rightname_list[which(!rightname_list %in% right_drop_list)]
}

cl = makeCluster(CPU_num)
clusterExport(cl, c('right_bed', 'right_flank_dir', 'genome_size_df', 'bedtools_path', 'HALF_DIST'), envir = .GlobalEnv)
clusterEvalQ(cl, list(library(bedtoolsr), options(bedtools.path = bedtools_path, scipen = 200)))

parSapply(cl, rightname_list, function(x){
  sub_rightbed = right_bed[which(right_bed$name==x), ]
  
  ## To merge tandem repeats
  sub_rightbed = sub_rightbed[order(sub_rightbed$chrmid, sub_rightbed$start), ]
  sub_rightbed = bt.merge(sub_rightbed, s=TRUE, d = 100, c="4,5,6", o="first,max,first")
  colnames(sub_rightbed)=c('chrmid', 'start', 'stop', 'name', 'score', 'strand')
  
  ## To extend bed files
  right_bed$start=as.numeric(right_bed$start)
  right_bed$stop=as.numeric(right_bed$stop)
  right_bed$stop=ifelse(right_bed$strand=='+', right_bed$stop, right_bed$start+HALF_DIST)
  right_bed$start=ifelse(right_bed$strand=='+', right_bed$stop-HALF_DIST, right_bed$start)
  right_bed$start=ifelse(right_bed$start<=0, 1, right_bed$start)
  sub_rightbed = merge(sub_rightbed, genome_size_df, by='chrmid')
  sub_rightbed$stop = ifelse(sub_rightbed$stop <= sub_rightbed$length, sub_rightbed$stop, sub_rightbed$length)
  
  ## To output 
  sub_rightbed = sub_rightbed[, c('chrmid', 'start', 'stop', 'name', 'score', 'strand')]
  sub_rightfile = paste(right_flank_dir, '/', x, '.bed', sep = '')
  bt.sort(sub_rightbed, output = sub_rightfile)
})
stopCluster (cl)

rm(right_bed)
gc()

## To create file path file

## if file is empty, then organize all possibilities
if(file.info(combinid_file)$size == 0){
  combine_id = t(do.call('cbind', lapply(leftname_list,function(x){
    sapply(rightname_list, function(y){c(x, y)})})))
  combine_id = as.data.frame(combine_id, stringsAsFactors=F)
  colnames(combine_id)=c('leftname', 'rightname')
}else{
  combine_id=combine_id[which(combine_id$leftname %in% leftname_list & combine_id$rightname %in% rightname_list), ]
}

combine_id$leftname=paste(subed_dir, '/leftflank/', combine_id$leftname, '.bed',sep = '')
combine_id$rightname=paste(subed_dir, '/rightflank/', combine_id$rightname, '.bed', sep = '')
write.table(combine_id, file = paste(subed_dir, '/', 'left_right.path.join', sep = ''), 
            quote = F, sep="\t", col.names = F, row.names = F)
