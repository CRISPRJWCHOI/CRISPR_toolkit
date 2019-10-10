#!/bin/bash


####################
## User parameter ##
####################################

user=SH
project=p53_screening
flash=FLASH-1.2.11-Linux-x86_64
thread=4

####################################


while read python_path;do
    python=$python_path
done < ../PythonPath.txt

nohup $python ./Flash_pair_read_merge.py $user $project $flash $thread > ./Output/${user}/${project}/Log/flash_log.txt 2>&1 &
