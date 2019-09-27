#!/bin/bash


####################
## User parameter ##
####################################

user=JaeWoo
project=JaeWoo_test_samples


####################################





while read python_path;do
    python=$python_path
done < ../PythonPath.txt

$python ./BaseEdit_input_converter.py $user $project
