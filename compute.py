# -*- coding: utf-8 -*-
"""
Created on Sat Oct  7 19:32:12 2017

@author: yf182
"""
import math as math
import numpy as np
from scipy.optimize import fsolve
#define functions
def get_flamda(lamda):#计算流量用到
    return (1+math.pow(lamda,2))*math.pow(1-(k-1)/(k+1)*math.pow(lamda,2),1/(k-1))

def get_rlamda(lamda):#计算各个截面上的静压用到
    return (1-(k-1)/(k+1)*math.pow(lamda,2))/(1+math.pow(lamda,2))


#获取速度系数
def func(x): #积分函数
    if x==0:
        return 0
    else:
        zx=0.5*(x+1/x)
        rx=get_rlamda(x)
        return (1/math.pow(x,2)-1)/(math.pow(zx,2)*math.pow(rx,n)*ratio)
def integrate(a,b):#使用牛顿-科茨求和公式计算积分
    dy=(b-a)/50
    sum=0
    sumall=0
    for i in range(50):
        y=0+i*dy
        sum==(b-a)/8*(func(y)+3*func(y+dy)+3*func(y+2*dy)+func(y+3*dy))
        sumall=sum+sumall
    return sumall
    
def get_lamda(lamda_L,location):#获取不同位置处的速度系数用于以后的计算
    count=0
    replace=0
    equation=1
    while equation>0.0001:
        dcount=lamda_L/100
        replace=replace+integrate(count,count+dcount)
        count=count+dcount
        equation=replace*L-location*integrate(0,lamda_L)
    return count


#获取几何参量
def get_s(e):
    replace=tdata[1]*math.pi/tdata[0]
    ex=l*math.sin(replace)/math.cos(tdata[2]/2)-r
    if 0<=e<=r1:
             s=l*2*tdata[0]*(math.sin(replace)/math.sin(tdata[2]/2)+(1-tdata[1])*math.pi/tdata[0]+(r1+r)/l*(math.pi/2+math.pi/tdata[0]-tdata[2]/2-1/math.tan(tdata[2]/2))-(r1-e)/l*math.pi/tdata[0])
    elif r1<e<=ex:
             s=l*2*tdata[0]*(math.sin(replace)/math.sin(tdata[2]/2)+(1-tdata[1])*math.pi/tdata[0]+(e+r)/l*(math.pi/2+math.pi/tdata[0]-tdata[2]/2-1/math.tan(tdata[2]/2)))
    elif ex<e<=e1:
             s=l*2*tdata[0]*((1-tdata[1])*math.pi/tdata[0]+(e+r)/l*(math.pi/tdata[0]+math.asin(math.sin(replace)*l/(e+r))))
    return s
def get_Ap(e):
    replace=tdata[1]*math.pi/tdata[0]
    ex=l*math.sin(replace)/math.cos(tdata[2]/2)-r
    if 0<=e<=r1:
        ap1=tdata[0]*((1-tdata[1])*math.pi/tdata[0]+math.sin(replace)*(math.cos(replace)-math.sin(replace)/math.tan(tdata[2]/2)))
        ap2=2*tdata[0]*(e+r)/l*(math.sin(replace)/math.sin(tdata[2]/2)+(1-tdata[1])*math.pi/tdata[0])
        ap3=tdata[0]*math.pow((e+r)/l,2)*(math.pi/2+math.pi/tdata[0]-tdata[2]/2-1/math.tan(tdata[2]/2))
        ap4=tdata[0]*math.pow((r1-e)/l,2)*(tdata[2]/2+1/math.tan(tdata[2]/2)-math.pi/2)
        Ap=math.pow(l,2)*(ap1+ap2+ap3+ap4)
    elif r1<e<=ex:
        ap1=tdata[0]*((1-tdata[1])*math.pi/tdata[0]+math.sin(replace)*(math.cos(replace)-math.sin(replace)/math.tan(tdata[2]/2)))
        ap2=2*tdata[0]*(e+r)/l*(math.sin(replace)/math.sin(tdata[2]/2)+(1-tdata[1])*math.pi/tdata[0])
        ap3=tdata[0]*math.pow((e+r)/l,2)*(math.pi/2+math.pi/tdata[0]-tdata[2]/2-1/math.tan(tdata[2]/2))
        Ap=math.pow(l,2)*(ap1+ap2+ap3)
    elif ex<e<=e1:
        ap1=tdata[0]*(math.pow(1+(e+r)/l,2)*(1-tdata[1])*math.pi/tdata[0])
        ap2=tdata[0]*math.sin(replace)*(math.sqrt(math.pow((e+r)/l,2)-math.pow(math.sin(replace),2))+math.cos(replace))
        ap3=tdata[0]*math.pow((e+r)/l,2)*(replace+math.asin(math.sin(replace)*l/(e+r)))
        Ap=math.pow(l,2)*(ap1+ap2+ap3)
    return Ap


def get_lamdaL(J): #计算得到入L
    def func(x):
        return math.pow((k+1)/2,1/(k-1))*x*math.pow(1-(k-1)/(k+1)*math.pow(x,2),1/(k-1))-J
    return fsolve(func,0.01)

def get_ratio(ai):
    if erod==1:
        if ai<=72.9:
            return 1
        else:
          return 1.3128-1.3249*math.pow(10,-2)*ai+1.5527*math.pow(10,-4)*math.pow(ai,2)-4.3868*math.pow(10,-7)*math.pow(ai,3)
    elif erod==0:
        return 1
#更新e
def get_eupdate(e,dt,rspeed,m):
    dx=L/m
    Ab=0
    s=np.zeros(m)
    for j in range(m):
        de=rspeed[j]*dt
        if e[j]+de<=e1:
            e[j]=e[j]+de
        else:
            e[j]=e1
        s[j]=get_s(e[j])
        dAb=s[j]*dx
        Ab=Ab+dAb
        if j==m-1:
            Ap=get_Ap(e[j])
    return e,Ab,Ap

def get_Vg(Vg,mb,dt):
    return Vg+mb/density*dt

def get_rspeed(lamda_L,p,ratio,m):
    dx=L/m
    speed=np.zeros(m)
    for j in range(m):
        location=j*dx
        lamda=get_lamda(location,lamda_L)
        rlamda=get_rlamda(lamda)
        speed[j]=ratio*a*math.pow(p*rlamda,n)
    return speed

def get_mb(e,rspeed,m):
    mb=0
    dx=L/m
    s=np.zeros(m)#内部使用，内部定义
    for j in range(m):
        if e[j]==e1:
            s[j]=0.0
        else:
            s[j]=get_s(e[j])
        dAb=s[j]*dx
        dmb=density*dAb*rspeed[j]
        mb=mb+dmb
    return mb

def get_md(p,lamda_L):
     flamda=get_flamda(lamda_L)
     return At*p/cx/flamda
 

def compute_k(Vg,mb,md):
    return R/M*Tf/Vg*(mb-md)*dt
    
def update(p,e,lamda_L,ratio,m,Vg,dt):#dt后数据?此时的压力与几何是否密切相关？
    rspeed=get_rspeed(lamda_L,p,ratio,m)
    tup=get_eupdate(e,dt,rspeed,m)
    e=tup[0]
    Ab=tup[1]
    Ap=tup[2]
    tup2=get_parametre(Ab,Ap)
    ratio=tup2[2]
    lamda_L=tup2[3]
    rspeed=get_rspeed(lamda_L,p,ratio,m)
    mb=get_mb(e,rspeed,m)
    Vg=get_Vg(Vg,mb,dt)
    md=get_md(p,lamda_L)
    return e,Ab,Ap,ratio,lamda_L,rspeed,Vg,mb,md
    

def get_parametre(Ab,Ap):
     ai=Ab/Ap
     J=At/Ap
     ratio=get_ratio(ai)
     lamda_L=get_lamdaL(J)
     return ai,J,ratio,lamda_L
 
     
#define variables
tdata=np.zeros(3)
with open('original2.data','r') as file1:
    variables=file1.readlines()
    D=float(variables[0])
    Iz=float(variables[1])
    F=float(variables[2])
    pc=float(variables[3])
    k=float(variables[4])
    M=float(variables[5])
    Tfc=float(variables[6])
    density=float(variables[7])
    rspeedfo=float(variables[8])
    n=float(variables[9])
    r1=float(variables[10])
    tdata[0]=float(variables[11])
    tdata[1]=float(variables[12])
    tdata[2]=float(variables[13])
    L=float(variables[14])
    erod=int(variables[15])
pe=101325#pa
R=8.314#j/mol/K
Tf=Tfc+273.15
a=rspeedfo/math.pow(6.9*10**6,n)
rs=a*math.pow(pc,n)
w=math.sqrt(k)*math.pow(2/(k+1),(k+1)/2/(k-1)) #
cf=w*math.sqrt(2*k/(k-1)*(1-math.pow(pe/pc,(k-1)/k)))
cx=math.sqrt(R/M*Tf)/w
Isp=cx*cf/9.8
At=F/cf/pc
kn=math.pow(pc,1-n)/(cx*density*a)
S=kn*At
mpeff=1.05*Iz/(9.8*Isp)
e1=mpeff/(density*S)
r=0.02*D#m
l=D/2-e1-r
y0=r/l
y1=(e1+r)/l


#给定初始条件
number=50#空间上的划分
e=np.zeros(number)
s=get_s(e[0])
Ab=s*L
Ap=get_Ap(e[-1])
ai=Ab/Ap
J=At/Ap
lamda_L=get_lamdaL(J)
flamdaL=get_flamda(lamda_L)
p0=0.1*math.pow(density*cx*a*Ab/At*flamdaL,1/(1-n))
rspeed=a*math.pow(p0,n)*np.ones(number)
mb=density*Ab*rspeed[0]
md=At*p0/cx
Vg=Ap*L
ratio=1
rspeed=a*math.pow(p0,n)*np.ones(number)
#计算压强随时间的变化
if erod==1:
    sting='pe.data'
elif erod==0:
        sting='pnoe.data'
    
with open(sting,'a') as f:
    f.write(str(p0)+'      '+str(0)+'\n')
t=0
p=p0
while(p>=p0):
    if t<=1:
        dt=0.001
    else:
        dt=0.01
    t=t+dt
    answer1=update(p,e,lamda_L,ratio,number,Vg,dt)
    p=p+R/M*Tf/Vg*(mb-md)*dt
    e=answer1[0]
    Ab=answer1[1]
    Ap=answer1[2]
    ratio=answer1[3]
    lamda_L=answer1[4]
    rspeed=answer1[5]
    Vg=answer1[6]
    mb=answer1[7]
    md=answer1[8]
    print('now the time is {0}'.format(t))
    with open(sting,'a') as f:
        f.write(str(p)+'     '+str(t)+'\n')


    
    


       
