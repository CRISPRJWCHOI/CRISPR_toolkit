#!/bin/bash

nohup ./Run_BaseEdit_freq_ver1.0.py -t 30 -w 16-48 --indel_check_pos 39-40 --target_ref_alt A,T --PAM_seq NGG --PAM_pos 43-45 --Guide_pos 23-42 --gap_open 20 --gap_extend 1 --end_open 20 --end_extend 1 &
