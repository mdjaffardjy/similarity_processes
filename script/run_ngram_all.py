#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 15:34:09 2022

@author: marinedjaffardjy
"""


### this script computes levenshtein for nf proc with tools, nf proc without tools, snk proc with tools
### it also computes levenshtein for workflows -- using batches (?)
import json
import jellyfish
import numpy
import pandas as pd
import math
import time
import ngram

################
# DEFINING VARIABLES
################

################
# IMPORTING THE FILES
################

#import snakemake rules
#import snakemake rule with tools
with open('json/snk_rule_info_tool.json') as f:
    snk_rules_tools = json.load(f)
    
with open('json/nf_proc_tool.json') as f:
    nf_rules_tools = json.load(f)    

################
# DEFINING THE FUNCTIONS
################    
    
#compute levenshtein for snakemake
def ngram_proc(rules, resume=0, output_file = "/home/marinedjaffardjy/Documents/Code/Similarite_process/json/ngram_snk_tool",n_size=3):
    #input : 
    #    list of snakemake rules in form of a json
    #    resume : int at which the computation will be restarted
    #    output_file : string model of the names for the output files
    #    output_file_resume : string file for the scores already computed in the case of a resume
    #output : table of dict with scores. also saves the scores for all pairs of processes in json (path )
    v1 = rules.copy()[resume:-1]
    v2 = rules.copy()[resume+1:]
    i = 0
    scores = []
    #if(resume>0):
        #with open(outputfile_resume) as f:
            #scores = json.load(f)
            #f.close
    #else :
        #scores = []
    time1= time.time()
    for rule1 in v1:
        print(str(i)+"/"+str(len(v1)))
        i+=1
        for rule2 in v2:
            score = ngram.NGram.compare(rule1['code'],rule2['code'],N=n_size)
            scores.append({"rule1":rule1["wf_orig"]+"/"+rule1["name"],
                           "rule2":rule2["wf_orig"]+"/"+rule2["name"],
                           "ngram":score})
        if(len(v2)>=2):
            v2 = v2[1:]
        if(i%50==0 or i ==len(v1) ):
            with open(output_file+"_"+str(i+resume)+".json","w") as f:
                print(output_file+"_"+str(i+resume)+".json")
                json.dump(scores,f)
                f.close
                scores = []
            
            print("time elapsed : "+str(time.time()-time1)+" sec ")
    return scores


#### launch script

snk_proc_path = "/home/marinedjaffardjy/Documents/Code/Similarite_process/json/ngram_snk_tools/ngram_3_snk_tools"
nf_proc_path = "/home/marinedjaffardjy/Documents/Code/Similarite_process/json/ngram_nf_tools/ngram_3_nf_tools"

#test ngram

#snk_scores = ngram_proc(snk_rules_tools[:60],0,snk_proc_path)


#snk tools proc

#snk_scores = ngram_proc(snk_rules_tools,0,snk_proc_path)

#nf tools proc

nf_scores = ngram_proc(nf_rules_tools,0,nf_proc_path)



#snk_no_tools_scores = levenshtein_proc(snk_rules_no_tools,0,snk_all_path)

#nf_no_tools_scores = levenshtein_proc(nf_rules_no_tools,0,nf_all_path)




















