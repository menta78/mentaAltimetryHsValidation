#!/bin/bash
#BSUB -q p_short
#BSUB -P 0419
#BSUB -n 4
#BSUB -J menta_exp0
#BSUB -R "span[ptile=36]"
#BSUB -o menta_exp0_%J.out
#BSUB -e menta_exp0_%J.err


py=/users_home/opa/lm09621/usr/miniconda3/bin/python

$py validateWwmUOST_ICE.py > run_validateWwmUOST_ICE.log 2>&1

