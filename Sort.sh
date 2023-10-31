inputfile=${1}
CPU_num=${2}

sort -k1,1 -k2,2n --parallel=${CPU_num} ${inputfile} > ${inputfile}.tmp

mv ${inputfile}.tmp ${inputfile}
