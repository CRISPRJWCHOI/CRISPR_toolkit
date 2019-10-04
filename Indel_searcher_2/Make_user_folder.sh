#!/bin/bash

user=JaeWoo
project=JaeWoo_test_samples

[ ! -d ./Input ] && { `mkdir ./Input`; }
[ ! -d ./User ] && { `mkdir ./User`; }
[ ! -d ./Output ] && { `mkdir ./Output`; }

[ ! -d ./Input/${user} ] && { `mkdir ./Input/${user}`; }
[ ! -d ./Input/${user}/FASTQ ] && { `mkdir ./Input/${user}/FASTQ`; }
[ ! -d ./Input/${user}/FASTQ/${project} ] && { `mkdir ./Input/${user}/FASTQ/${project}`; }
[ ! -d ./Input/${user}/Reference ] && { `mkdir ./Input/${user}/Reference`; }
[ ! -d ./Input/${user}/Reference/${project} ] && { `mkdir ./Input/${user}/Reference/${project}`; }

[ ! -d ./User/${user} ] && { `mkdir ./User/${user}`; }
> ./User/${user}/${project}.txt
