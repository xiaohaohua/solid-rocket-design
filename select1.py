# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 15:54:02 2017

@author: yf182
"""

import numpy as np
from scipy.optimize import fsolve
import math as math
import sys as sys
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
#绘图设置
params={
        'axes.labelsize':'35',
        'xtick.labelsize':'27',
        'ytick.labelsize':'27',
        'lines.linewidth':2,
        'legend.fontsize':'27',
        'figure.figsize':'12,9'
        }
pylab.rcParams.update(params)
#define variables
with open('original.data','r') as file1:
    variables=file1.readlines()
    D=float(variables[0])
    Iz=float(variables[1])
    F=float(variables[2])
    pc=float(variables[3])
    k=float(variables[4])
    M=float(variables[5])
    Tfc=float(variables[6])
    density=float(variables[7])
    rspeed=float(variables[8])
    n=float(variables[9])
    r1=float(variables[10])
pe=101325#pa
R=8.314#j/mol/K
Tf=Tfc+273.15
a=rspeed/math.pow(6.9*10**6,n)
rs=a*math.pow(pc,n)
datasize=0
successsize=0
class  DynamicArray(object):
    def __init__(self,item_type):
        self._data=np.zeros(10,dtype=item_type)
        self._size=0
    def get_data(self,index):
        return self._data[index]
    def append(self,value):
        if len(self._data)==self._size:
            self._data=np.resize(self._data,int(len(self._data)*1.25))
        self._data[self._size]=value
        self._size+=1
item_type=np.dtype({"names":["n","s","verify1","angle0d2","sdl1","sdlmin","Vertify2"],"formats":["i4","f8","f8","f8","f8","f8","f8"]})       
da=DynamicArray(item_type)
item_type2=np.dtype({"names":["n","s","angle","L"],"formats":["i4","f8","f8","f8"]})
ssdata=DynamicArray(item_type2)
#热力参量计算
w=math.sqrt(k)*math.pow(2/(k+1),(k+1)/2/(k-1)) #
cf=w*math.sqrt(2*k/(k-1)*(1-math.pow(pe/pc,(k-1)/k)))
c=math.sqrt(R/M*Tf)/w
Isp=c*cf/9.8
At=F/cf/pc
kn=math.pow(pc,1-n)/(c*density*a)
S=kn*At
mpeff=1.05*Iz/(9.8*Isp)
e1=mpeff/(density*S)
r=0.02*D#m
l=D/2-e1-r
y0=r/l
y1=(e1+r)/l
if  e1/rs<4.0 or e1/rs>7.0:
    print('燃速不符合,请调整参数')
    print('当前e1为'+str(e1))
    sys.exit(0)
if  y1<0.8 or y1>1.2:
    print('y1不符合，请调整参数')
    print('当前e1为'+str(e1))
    print('当前y1为'+str(y1))
    sys.exit(0)
with open('process.data','a') as f2:
    f2.write(str(mpeff)+'------'+'有效装药质量'+'\n')
    f2.write(str(At)+'------'+'At'+'\n')
    f2.write(str(kn)+'------'+'kn'+'\n')
    f2.write(str(S)+'------'+'S'+'\n')
    f2.write(str(e1)+'------'+'e1'+'\n')
    f2.write(str(l)+'------'+'l'+'\n')
    f2.write(str(y0)+'------'+'y0'+'\n')
    f2.write(str(y1)+'------'+'y1'+'\n')
    f2.write(str(Isp)+'------'+'Isp'+'\n')
    f2.write('\n')
#第一次选择
print('第一次选择')
with open('process.data','a') as f:
    f.write('第一次选择'+'\n')
    f.write('\n')
for i in range(3,9):
    for j in range(3,9):
        s=j*0.1;
        vertify1=math.sin(s*math.pi/i)/y1
        if vertify1<=1:
           sdl1=2*i*((1-s)*math.pi/i+y1*(math.pi/i+math.asin(math.sin(s*math.pi/i)/y1)))
           def func1(x):
               return x/2+1.0/math.tan(x/2)-math.pi/2-math.pi/i
           angle0=fsolve(func1,0.001)
           sdlmin=2*i*(math.sin(s*math.pi/i)/math.sin(angle0/2)+(1-s)*math.pi/i)
           vertify2=sdl1/sdlmin
           if vertify2<=1.2:
               da.append((i,s,vertify1,angle0/2,sdl1,sdlmin,vertify2))
               datasize=datasize+1
               with open('process.data','a') as f:
                    f.write('right,第一次筛选成功'+'\n')
                    f.write('星角数n为'+'------'+str(i)+'\n')
                    f.write('角度系数为'+'-------'+str(s)+'\n')
                    f.write('反三角函数判据'+'------'+str(vertify1)+'\n')
                    f.write('后期最小燃面对应星边角/2'+'------'+str(angle0/2)+'\n')
                    f.write('后期最小燃面值'+'------'+str(sdlmin)+'\n')
                    f.write('后期最大燃面值'+'------'+str(sdl1)+'\n')
                    f.write('增面比'+'-----'+str(vertify2)+'\n')
                    f.write('\n')
           else:
               with open('process.data','a') as f:
                   f.write('false,增面比不符合要求'+'\n')
                   f.write('星角数n为'+'------'+str(i)+'\n')
                   f.write('角度系数为'+'-------'+str(s)+'\n')
                   f.write('反三角函数判据'+'------'+str(vertify1)+'\n')
                   f.write('后期最小燃面对应星边角/2'+'------'+str(angle0/2)+'\n')
                   f.write('后期最小燃面值'+'------'+str(sdlmin)+'\n')
                   f.write('后期最大燃面值'+'------'+str(sdl1)+'\n')
                   f.write('增面比'+'-----'+str(vertify2)+'\n')
                   f.write('\n')
        else:
            with open('process.data','a') as f:
                f.write('false,不满足反三角函数要求'+'\n')
                f.write('星角数n为'+'------'+str(i)+'\n')
                f.write('角度系数为'+'-------'+str(s)+'\n')
                f.write('反三角函数判据'+'------'+str(vertify1)+'\n')
                f.write('\n')
if datasize==0:
    print('第一次筛选全部失败,请调整参数')
    sys.exit(0)
#角的求解
angle=np.zeros(datasize)
for i in range(0,datasize):
    data=da.get_data(i)
    def func2(x):
            return 2*data[0]*(math.sin(data[1]*math.pi/data[0])/math.sin(x/2)+(1-data[1])*math.pi/data[0]+(r1+r)/l*(math.pi/2+math.pi/data[0]-x/2-1/math.tan(x/2)))-data[4]
    angle[i]=fsolve(func2,0.1)
#注意python中数组的索引取值从0开始
#通气面积通气参量，装填系数计算
print('装填与剩药选药开始')
with open('process.data','a') as f:
    f.write('装填与剩药选药开始'+'\n')
Ac=0.25*math.pi*math.pow(D,2)
for i in range(0,datasize):
    data=da.get_data(i)
    if data[1]*math.pi/data[0]>=angle[i]/2:
        with open('process.data','a') as f:
            f.write('第二次筛选失败，判据3不符合'+'\n')
            f.write('星角数为'+'-----'+str(data[0])+'\n')
            f.write('星角系数为'+'-----'+str(data[1])+'\n')
            f.write('前期星边夹角'+'-----'+str(angle[i])+'\n')
            f.write('\n')
    else:
        replace=data[1]*math.pi/data[0]
        ap1=data[0]*((1-data[1])*math.pi/data[0]+math.sin(replace)*(math.cos(replace)-math.sin(replace)/math.tan(angle[i]/2)))
        ap2=2*data[0]*y0*(math.sin(replace)/math.sin(angle[i]/2)+(1-data[1])*math.pi/data[0])
        ap3=data[0]*math.pow(y0,2)*(math.pi/2+math.pi/data[0]-angle[i]/2-1/math.tan(angle[i]/2))
        ap4=data[0]*math.pow((r1/l),2)*(angle[i]/2+1/math.tan(angle[i]/2)-math.pi/2)
        Ap0=math.pow(l,2)*(ap1+ap2+ap3+ap4)
        J=At/Ap0
        vertify3=1-Ap0/Ac
        if 0.75<=vertify3<=0.85:
            #计算剩药与剩药系数
            af1=data[1]*math.pi*math.pow(1+y1,2)
            af2=-data[0]*(math.sin(replace)*(math.sqrt(math.pow(y1,2)-math.pow(math.sin(replace),2))+math.cos(replace)))
            af3=-data[0]*math.pow(y1,2)*(replace+math.asin(math.sin(replace)/y1))
            Af=math.pow(l,2)*(af1+af2+af3)
            vertify4=Af/Ac
            if vertify4<=0.05:
                L=mpeff/density/(Ac-Ap0-Af)
                ai=S/Ap0
                if 1.3<=L<=1.8:
                    with open('process.data','a') as f:
                                    f.write('第二次筛选成功，各个标准全部合格'+'\n')
                                    f.write('星角数为'+'-----'+str(data[0])+'\n')
                                    f.write('星角系数为'+'-----'+str(data[1])+'\n')
                                    f.write('前期星边夹角'+'-----'+str(angle[i])+'\n')
                                    f.write('装填系数为'+'-----'+str(vertify3)+'\n')
                                    f.write('剩药系数为'+'-----'+str(vertify4)+'\n')
                                    f.write('长度为'+'------'+str(L)+'\n')
                                    f.write('\n') 
                                    ssdata.append((data[0],data[1],angle[i],L))
                                    successsize=successsize+1
                else:
                    with open('process.data','a') as f:
                                f.write('第二次筛选失败，装填系数合格，剩药合格，长度不合格'+'\n')
                                f.write('星角数为'+'-----'+str(data[0])+'\n')
                                f.write('星角系数为'+'-----'+str(data[1])+'\n')
                                f.write('前期星边夹角'+'-----'+str(angle[i])+'\n')
                                f.write('装填系数为'+'-----'+str(vertify3)+'\n')
                                f.write('剩药系数为'+'-----'+str(vertify4)+'\n')
                                f.write('长度为'+'------'+str(L)+'\n')
                                f.write('\n')   
            else:
                with open('process.data','a') as f:
                        f.write('第二次筛选失败，装填系数合格，剩药不合格'+'\n')
                        f.write('星角数为'+'-----'+str(data[0])+'\n')
                        f.write('星角系数为'+'-----'+str(data[1])+'\n')
                        f.write('前期星边夹角'+'-----'+str(angle[i])+'\n')
                        f.write('装填系数为'+'-----'+str(vertify3)+'\n')
                        f.write('剩药系数为'+'-----'+str(vertify4)+'\n')
                        f.write('\n')
        else:
            with open('process.data','a') as f:  
                f.write('第二次筛选失败，装填系数不合格'+'\n')
                f.write('星角数为'+'-----'+str(data[0])+'\n')
                f.write('星角系数为'+'-----'+str(data[1])+'\n')
                f.write('前期星边夹角'+'-----'+str(angle[i])+'\n')
                f.write('装填系数为'+'-----'+str(vertify3)+'\n')
                f.write('\n')  
#筛选判断  
if successsize==0:
    print('第二次筛选失败，没有合格组，请调整参数'+'\n')
    sys.exit(0)
#燃面曲线绘制
for i in range(0,successsize):
    tdata=ssdata.get_data(i)
    replace=tdata[1]*math.pi/tdata[0]
    ex=l*math.sin(replace)/math.cos(tdata[2]/2)-r
    def Sface(e):
        if 0<=e<=r1:
             s=l*2*tdata[0]*(math.sin(replace)/math.sin(tdata[2]/2)+(1-tdata[1])*math.pi/tdata[0]+(r1+r)/l*(math.pi/2+math.pi/tdata[0]-tdata[2]/2-1/math.tan(tdata[2]/2))-(r1-e)/l*math.pi/tdata[0])
        elif r1<e<=ex:
             s=l*2*tdata[0]*(math.sin(replace)/math.sin(tdata[2]/2)+(1-tdata[1])*math.pi/tdata[0]+(e+r)/l*(math.pi/2+math.pi/tdata[0]-tdata[2]/2-1/math.tan(tdata[2]/2)))
        elif ex<e<=e1:
             s=l*2*tdata[0]*((1-tdata[1])*math.pi/tdata[0]+(e+r)/l*(math.pi/tdata[0]+math.asin(math.sin(replace)*l/(e+r))))
        return tdata[3]*s
    xdata=np.linspace(0,e1,num=200)
    ydata=[Sface(element) for element in xdata]
    plt.figure(1)
    string=str(tdata[0])+"--"+str(tdata[1])+"--"+str(tdata[2])
    plt.plot(xdata,ydata,label=str(string))
    plt.legend(loc='upper center',labelspacing=0.25,bbox_to_anchor=(0.6,1))#必须添加图例，否则不会显示
    plt.xlabel('e')
    plt.ylabel('Ab')
    plt.savefig('graph1',dip=1000)
                 