#!/bin/bash


####################
## User parameter ##
####################################

user=SH
project=p53_screening
thread=2


####################################





while read python_path;do
    python=$python_path
done < ../PythonPath.txt


nohup $python ./Summary_Random_barcode.py -u $user -p $project -t $thread > ./Output/${user}/${project}/Log/Random_barcode_log.txt 2>&1 &
