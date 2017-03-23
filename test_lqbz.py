# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 17:16:16 2017

@author: lu
"""

import math 
import numpy as np
import pandas as pd
from scipy.stats import rankdata

'''
a = np.array([[1,2,2,3,3,1],
              [7,9,6,8,3,2],
              [6,1,1,1,1,5],
              [4,3,2,6,7,8]])
a = np.array([[1,1,1,1,1,1],[1,1,1,1,1,1],[1,1,1,1,1,1]])
print a.shape 


ranks = [rankdata(a[:,j]/4.0,'dense') for j in range(6)]
for i in ranks :
    print i 
q1 = np.percentile(ranks,95, axis=1) 
q2 = np.percentile(ranks,85, axis=1) 
for i in q1 :
    print i
for i in q2 :
    print i
'''
fr = open(r't.txt','r')
set1 = set()
for line in fr :
    line_str = line.strip().split('\t')
    set1.add(line_str[4])
for i in set1 :
    print i
    