#!/bin/bash

####################
## User parameter ##
###################################

user=SH
project=24K_screening

###################################



while read python_path;do
    python=$python_path
done < ../PythonPath.txt

nohup $python ./Each_base_summary.py $user $project > ./Output/${user}/${project}/Log/Each_base_summary_log.txt 2>&1 & 
