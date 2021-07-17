#!/usr/bin/env python3


from numpy import *
import pandas as pd
import os
from andromedaPy import galaxy
if os.getenv('LAUNCH')==None :
    print('warning,  bootup Jupyter after ignition.sh')
from scipy import linalg
import itertools

m = [-1,0,1]
n = [0.5]
sg = pd.DataFrame(data = list(itertools.product(m,n,m,n,m,n,m,n,m,n,m,n)))
sg['index'] = sg.index
sg.to_csv('sg.csv',index = False)
