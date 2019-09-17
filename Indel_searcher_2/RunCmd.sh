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


python=/media/hkim/Pipeline/Indel_searcher_2/miniconda2/bin/python
EDNAFULL=/media/hkim/Pipeline/Indel_searcher_2/miniconda2/pkgs/crispresso2-2.0.30-py27h14c3975_0/lib/python2.7/site-packages/CRISPResso2/EDNAFULL

[ ! -d ./Output/${user} ] && { `mkdir ./Output/${user}`; }
[ ! -d ./Output/${user}/${project} ] && { `mkdir ./Output/${user}/${project}`; }
[ ! -d ./Output/${user}/${project}/Log ] && { `mkdir ./Output/${user}/${project}/Log`; }

$python ./Run_indel_searcher.py --python $python --user $user --project $project --ednafull $EDNAFULL --pam_type $pam_type --pam_pos $pam_pos -t $thread #> ./Output/${user}/${project}/Log/log.txt 2>&1 &
