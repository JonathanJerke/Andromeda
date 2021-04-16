from numpy import *
import pandas as pd
import os
from andromedaPy import galaxy
if os.getenv('LAUNCH')==None :
    print('warning,  bootup Jupyter after ignition.sh')
from scipy import linalg
import itertools

l = [0]
b = [0.2,1,5,25]
x = [0]
gto = pd.DataFrame(data = list(itertools.product(l,b,x,l,b,x,l,b,x)))
gto['index'] = gto.index
gto.to_csv('gto.csv',index = False)
