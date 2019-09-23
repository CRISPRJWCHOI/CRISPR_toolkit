#!/bin/bash

####################
## User parameter ##
###################################

user=JaeWoo
project=JaeWoo_test_samples4
target_window=16-48
indel_check_pos=39-40
target_ref_alt=A,G
PAM_seq=NGG
PAM_pos=43-45
Guide_pos=23-42

thread=15

gap_open=-10 ## default
gap_extend=1 ## default

###################################

python=/media/hkim/Pipeline/Indel_searcher_2/miniconda2/bin/python

[ ! -d ./Output/${user} ] && { `mkdir ./Output/${user}`; }
[ ! -d ./Output/${user}/${project} ] && { `mkdir ./Output/${user}/${project}`; }
[ ! -d ./Output/${user}/${project}/Log ] && { `mkdir ./Output/${user}/${project}/Log`; }

$python ./Run_BaseEdit_freq.py --python $python --user $user --project $project --ednafull $EDNAFULL -w $target_window --indel_check_pos $indel_check_pos\
                               --target_ref_alt $target_ref_alt --PAM_seq $PAM_seq --PAM_pos $PAM_pos --Guide_pos $Guide_pos \
                               --gap_open $gap_open --gap_extend $gap_extend -t $thread > ./Output/${user}/${project}/Log/log.txt 2>&1 &
