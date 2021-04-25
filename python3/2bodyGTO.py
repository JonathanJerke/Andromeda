#!/usr/bin/env python3


from numpy import *
import pandas as pd
import os
from andromedaPy import galaxy
if os.getenv('LAUNCH')==None :
    print('warning,  bootup Jupyter after ignition.sh')
from scipy import linalg
import itertools

l = [0]
b = [1.2,2,2.8]
x = [0]
gto = pd.DataFrame(data = list(itertools.product(l,b,x,l,b,x,l,b,x,l,b,x,l,b,x,l,b,x)))
gto['index'] = gto.index
gto.to_csv('gto.csv',index = False)
