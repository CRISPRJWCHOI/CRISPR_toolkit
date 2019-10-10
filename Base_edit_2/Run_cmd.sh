#!/bin/bash

####################
## User parameter ##
###################################

user=JaeWoo
project=JaeWoo_test_samples
target_window=20-59
indel_check_pos=50-51
target_ref_alt=A,G
PAM_seq=NGG
PAM_pos=54-56
Guide_pos=23-53

thread=15

gap_open=-10 ## default
gap_extend=1 ## default

###################################

while read python_path;do
    python=$python_path
done < ../PythonPath.txt

[ ! -d ./Output/${user} ] && { `mkdir ./Output/${user}`; }
[ ! -d ./Output/${user}/${project} ] && { `mkdir ./Output/${user}/${project}`; }
[ ! -d ./Output/${user}/${project}/Log ] && { `mkdir ./Output/${user}/${project}/Log`; }

nohup $python ./Run_BaseEdit_freq.py --python $python --user $user --project $project -w $target_window --indel_check_pos $indel_check_pos \
                               --target_ref_alt $target_ref_alt --PAM_seq $PAM_seq --PAM_pos $PAM_pos --Guide_pos $Guide_pos \
                               --gap_open $gap_open --gap_extend $gap_extend -t $thread > ./Output/${user}/${project}/Log/log.txt 2>&1 &
