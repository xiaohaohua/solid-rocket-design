# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 10:07:40 2017

@author: yf182
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
#绘图设置
params={
        'axes.labelsize':'35',
        'xtick.labelsize':'27',
        'ytick.labelsize':'27',
        'lines.linewidth':2,
        'legend.fontsize':'27',
        'figure.figsize':'12,12'
        }
pylab.rcParams.update(params)

with open('pe.data','r') as f1:
     string=f1.readlines();
     lines=len(string)
     data=np.zeros((lines,2))
     for i in range(lines):
         middle=string[i]
         data[i,:]=middle.split()

t=np.zeros(data.shape[0])
p=np.zeros(data.shape[0])

for i in range(data.shape[0]):
    t[i]=data[i,1]
    p[i]=data[i,0]
plt.figure(1)
    
plt.plot(t,p,label='erod')
plt.legend(loc='best',labelspacing=0.25)

plt.xlabel('time')
plt.ylabel('pressue')
with open('pnoe.data','r') as f1:
     string=f1.readlines();
     lines=len(string)
     data=np.zeros((lines,2))
     for i in range(lines):
         middle=string[i]
         data[i,:]=middle.split()

t=np.zeros(data.shape[0])
p=np.zeros(data.shape[0])

for i in range(data.shape[0]):
    t[i]=data[i,1]
    p[i]=data[i,0]

plt.plot(t,p,label='no erod')
plt.legend(loc='best',labelspacing=0.25)

         
     
    
    


