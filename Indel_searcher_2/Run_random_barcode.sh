#!/bin/bash


####################
## User parameter ##
####################################

user=JaeWoo
project=JaeWoo_test_samples
thread=2


####################################





while read python_path;do
    python=$python_path
done < ../PythonPath.txt


nohup $python ./Summary_Random_barcode.py -u $user -p $project -t $thread > Random_barcode_log.txt 2>&1 &
