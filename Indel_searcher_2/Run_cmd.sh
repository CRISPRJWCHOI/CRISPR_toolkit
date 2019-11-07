#!/bin/bash

####################
## User parameter ##
###################################

user=JaeWoo
project=JaeWoo_test_samples
pam_type=Cas9
pam_pos=Forward
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

nohup $python ./Run_indel_searcher.py --python $python --user $user --project $project --pam_type $pam_type --pam_pos $pam_pos -t $thread > ./Output/${user}/${project}/Log/log.txt 2>&1 &
