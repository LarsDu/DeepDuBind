#!/usr/bin/bash
#nohup python -um run_dubind --json_params=EGR1_train.json --mode='eval' > EGR1_infC_eval.log 2>&1 & tail -f EGR1_infC_eval.log
nohup python -um run_dubind --json_params=CTCF.json --mode='eval' > CTCF_infC_eval.log 2>&1 & tail -f CTCF_infC_eval.log
nohup python -um run_dubind --json_params=GABPA.json --mode='eval' > GABPA_infC_eval.log 2>&1 & tail -f GABPA_infC_eval.log
