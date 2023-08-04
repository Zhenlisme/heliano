library(seqinr)

args=commandArgs(T)

wkdir = args[1]
OPTFILE = args[2]

setwd(wkdir)

file_list = list.files('./')

aln_list = file_list[sapply(file_list, function(x){endsWith(x, '.fa.aln')})]

calculate_distance = function(x){
  aln = read.alignment(x, format = 'fasta')
  distance = as.vector(dist.alignment(aln))
  ## means only one sequence
  if(length(distance)==0){
    distance = 1
  }
  ## na means the two sequences have very low identity.
  distance[is.na(distance)] = 1
  identity =1 - distance**2
  string_vector = strsplit(x, split = '\\.')[[1]]
  name = string_vector[1]
  direction = string_vector[2]
  return(c(name, direction, mean(identity)))
}


if(length(aln_list)>=2){
  identity_df = as.data.frame(do.call('rbind', lapply(aln_list, function(x){calculate_distance(x)})), stringsAsFactors = F)
  colnames(identity_df)=c('pairname', 'direction', 'identity')
#  identity_df = tidyr::spread(identity_df, key = direction, value = identity)
#  identity_df[which(is.na(identity_df$left)), 'left'] = 0
#  identity_df[which(is.na(identity_df$right)), 'right'] = 0

  write.table(identity_df, file = OPTFILE, sep = '\t', quote = F, col.names = T, row.names = F)
}

