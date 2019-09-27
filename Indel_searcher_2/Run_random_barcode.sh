#!/bin/bash


####################
## User parameter ##
####################################

user=JaeWoo
project=JaeWoo_test_samples
thread=15


####################################





while read python_path;do
    python=$python_path
done < ../PythonPath.txt


$python ./Summary_Random_barcode.py -u $user -p $project -t 15
