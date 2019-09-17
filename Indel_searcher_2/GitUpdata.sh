#!/bin/bash

git init
git add *.py RunCmd.sh MakeUserFolder.sh Kill_jobs.sh  
git commit -m "update"
git config --global user.email "abmrcbioinformatics@gmail.com"
git config --global user.name "CRISPRJWCHOI"
git remote add origin https://github.com/CRISPRJWCHOI/NGS_DNA_pipeline.git
git pull origin master
git push -u origin master
