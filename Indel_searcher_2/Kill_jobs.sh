#!/bin/bash

# Confirm the jobs.
# ps aux | grep hkim | grep BaseEdit_freq_ver1.0.py | less

kill -9 $(ps aux | grep hkim | grep Run_indel_searcher | awk '{print$2}')
kill -9 $(ps aux | grep hkim | grep Indel_searcher_crispresso_hash | awk '{print$2}')
