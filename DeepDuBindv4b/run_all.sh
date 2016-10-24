#!/usr/bin/bash 
                                                                                      
#nohup python -um run_dubind --json_params=ATF2.json --mode='train' > ATF2_infC_train.log 2>&1 & tail -f ATF2_infC_train.log
#nohup python -um run_dubind --json_params=CTCF.json --mode='train' > CTCF_infC_train.log 2>&1 & tail -f CTCF_infC_train.log
#nohup python -um run_dubind --json_params=E2F1.json --mode='train' > E2F1_infC_train.log 2>&1 & tail -f E2F1_infC_train.log
#nohup python -um run_dubind --json_params=EGR1.json --mode='train' > EGR1_infC_train.log 2>&1 & tail -f EGR1_infC_train.log
#nohup python -um run_dubind --json_params=FoxA1.json --mode='train' > FoxA1_infC_train.log 2>&1 & tail -f FoxA1_infC_train.log
#nohup python -um run_dubind --json_params=FoxA2.json --mode='train' > FoxA2_infC_train.log 2>&1 & tail -f FoxA2_infC_train.log
#nohup python -um run_dubind --json_params=GABPA.json --mode='train' > GABPA_infC_train.log 2>&1 & tail -f GABPA_infC_train.log
#nohup python -um run_dubind --json_params=HNF4A.json --mode='train' > HNF4A_infC_train.log 2>&1 & tail -f HNF4A_infC_train.log
#nohup python -um run_dubind --json_params=JUND.json --mode='train' > JUND_infC_train.log 2>&1 & tail -f JUND_infC_train.log
#nohup python -um run_dubind --json_params=MAX.json --mode='train' > MAX_infC_train.log 2>&1 & tail -f MAX_infC_train.log
#nohup python -um run_dubind --json_params=NANOG.json --mode='train' > NANOG_infC_train.log 2>&1 & tail -f NANOG_infC_train.log
nohup python -um run_dubind --json_params=REST.json --mode='train' > REST_infC_train.log 2>&1 & tail -f REST_infC_train.log
#nohup python -um run_dubind --json_params=TAF1.json --mode='train' > TAF1_infC_train.log 2>&1 & tail -f TAF1_infC_train.log


