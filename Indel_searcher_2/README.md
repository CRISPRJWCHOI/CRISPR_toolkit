# Indel_searcher_2
Fast CRISPR indel search tool

### Prerequisites to run
```
 install the miniconda2
 https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh

 Run the conda package manager
 conda config --add channels defaults
 conda config --add channels bioconda
 conda config --add channels conda-forge
 conda install CRISPResso2
 
 vi ~/.bashrc
 export PATH=$PATH:/path/to/minicodna2/bin
 
 vi RunCmd.sh
 python=path/to/miniconda2/bin/python2.7
 EDNAFULL=path/to/miniconda2/pkgs/crispresso2-2.0.30-py27h14c3975_0/lib/python2.7/site-packages/CRISPResso2/EDNAFULL
```
