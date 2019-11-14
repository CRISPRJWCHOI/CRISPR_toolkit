# Indel_searcher_2
Fast CRISPR indel search tool

### Prerequisites to run
```
 # install the miniconda2.
 https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh

 # Run the conda package manager.
 conda config --add channels defaults
 conda config --add channels bioconda
 conda config --add channels conda-forge
 conda install CRISPResso2
 
 vi ~/.bashrc
 export PATH=$PATH:/path/to/minicodna2/bin

 vi Make_user_folder.sh
 # Modify the user name and project name.
 user=JaeWoo
 project=JaeWoo_test_samples
 ./Make_user_folder.sh

 vi Run_cmd.sh
 # Modify the parameters. The user and project name must be the same as that used in the 'Make_user_folder.sh'.
 user=JaeWoo
 project=JaeWoo_test_samples
 pam_type=Cas9
 pam_pos=Forward
 thread=15
 ./Run_cmd.sh
```
