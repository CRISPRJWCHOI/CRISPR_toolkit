#!/bin/bash


####################
## User parameter ##
###################################

user=SH
project=24K_screening
window=25-34
thread=4

###################################



while read python_path;do
    python=$python_path
done < ../PythonPath.txt

nohup $python ./Sequence_freq.py $user $project $window $thread > ./Output/${user}/${project}/Log/Sequence_freq_log.txt 2>&1 &
