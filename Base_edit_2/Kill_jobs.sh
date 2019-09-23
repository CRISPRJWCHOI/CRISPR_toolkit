#!/bin/bash

# Confirm the jobs.
# ps aux | grep hkim | grep BaseEdit_freq_ver1.0.py | less

kill -9 $(ps aux | grep hkim | grep Run_BaseEdit_freq.py | awk '{print$2}')
kill -9 $(ps aux | grep hkim | grep BaseEdit_freq_crispresso.py | awk '{print$2}')
