#!/bin/bash
msmutect=/storage/bfe_maruvka/avrahamk/MSMuTect_v3.1

######
# Input
# $1 Pair_ID
# $2 Tumor  BAM. Full PATH
# $3 Normal BAM. Full PATH
# $4 Loci file list

#Calling Tumor alleles
python3 $msmutect/main.py -I $2 -l $4 -O $1.Tumor.hist
python3 $msmutect/reformat_histogram.py $1.Tumor.hist
python3 $msmutect/calculate_alleles.py  $1.Tumor.hist.mot $msmutect/data/probability_table.csv 


#Calling Normal alleles
python3 $msmutect/main.py -I $3 -l $4 -O $1.Normal.hist
python3 $msmutect/reformat_histogram.py $1.Normal.hist
python3 $msmutect/calculate_alleles.py  $1.Normal.hist.mot $msmutect/data/probability_table.csv 

#Taking the shared loci
sh $msmutect/shell/Shared_loci_v3.sh $1.Tumor.hist.mot.all  $1.Normal.hist.mot.all A


#Calling mutations
python3 $msmutect/Find_mutations2.py  $1.Tumor.hist.mot.all.tmp.par.reg  $1.Normal.hist.mot.all.tmp.par.reg $msmutect/data/probability_table.csv 8 0.3 0.031 > $1.mut

#Simple output format 

tr '\n' ' ' < $1.mut | awk '{gsub("@","\n");print $0}' > $1.mut.cln
        
python3 reformat_output.py $1.mut.cln $1.maf_like



