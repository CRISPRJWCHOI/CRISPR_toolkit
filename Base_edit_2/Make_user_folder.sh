#!/bin/bash

user=JaeWoo
project=JaeWoo_test_samples4


[ ! -d ./Input ] && { `mkdir ./Input`; }
[ ! -d ./User ] && { `mkdir ./User`; }
[ ! -d ./Output ] && { `mkdir ./Output`; }

[ ! -d ./Input/${user} ] && { `mkdir ./Input/${user}`; }
[ ! -d ./Input/${user}/Query ] && { `mkdir ./Input/${user}/Query`; }
[ ! -d ./Input/${user}/Query/${project} ] && { `mkdir ./Input/${user}/Query/${project}`; }
[ ! -d ./Input/${user}/Reference ] && { `mkdir ./Input/${user}/Reference`; }
[ ! -d ./Input/${user}/Reference/${project} ] && { `mkdir ./Input/${user}/Reference/${project}`; }

[ ! -d ./User/${user} ] && { `mkdir ./User/${user}`; }
> ./User/${user}/${project}.txt
