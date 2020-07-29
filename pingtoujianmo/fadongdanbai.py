# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 16:19:45 2020
这是建立发动蛋白的粗粒度模型，在运行程序前
请大声喊出：平头傻逼，否则程序会启动自毁程序
@author: Wang
"""
# import numpy as np
import math
d = 40 #圆片直径
m = 1  #螺旋圈数
r = 100 #孔半径
h = 50 #螺旋高度（必须大于等于round(m)*d）
n = 133 #原片数
c = 4  #圆片珠子层数
s = 28 #第一层珠子数
fw = open('luoxvan.xyz','w')
for i in range(n):
    z = i * h / n
    x = r * math.cos(math.pi*2*m*i/n)
    y = r * math.sin(math.pi*2*m*i/n)
    R = math.sqrt(x**2 + y**2 + z**2) #原点极径
    theta = math.asin(z/R) #z轴角度
    fi = math.acos(x/math.sqrt(x**2+y**2)) #x轴角度
    print(x,y,z,end = ' ', file = fw)
    print(file = fw)
    for j in range(c,0,-1):
        x_1 = x * (R-d/2*j/c)/R #近原点x
        y_1 = y * (R-d/2*j/c)/R #近原点y
        z_1 = z * (R-d/2*j/c)/R #近原点z
        x_2 = x * (R+d/2*j/c)/R #近原点x
        y_2 = y * (R+d/2*j/c)/R #近原点y
        z_2 = z * (R+d/2*j/c)/R #近原点z
        print(x_1,y_1,z_1,end = ' ', file = fw)
        print(file = fw)
        print(x_2,y_2,z_2,end = ' ', file = fw)
        print(file = fw)
        p = int(round((s-2)/2/math.pow(1.472,(c-j)))+1)
        for k in range(1,p):
            R1 = math.sqrt(R**2+(d/2*j/c)**2-R*d*j/c*math.cos(k*math.pi/p)) #圆盘极径
            beta = math.acos((R**2+R1**2-(d/2*j/c)**2)/2/R/R1) #两极径夹角
            x1 = R1*math.cos(theta-beta)*math.cos(math.pi*2*m*i/n)
            y1 = R1*math.cos(theta-beta)*math.sin(math.pi*2*m*i/n)
            z1 = R1*math.sin(theta-beta) 
            x2 = R1*math.cos(theta+beta)*math.cos(math.pi*2*m*i/n)
            y2 = R1*math.cos(theta+beta)*math.sin(math.pi*2*m*i/n)
            z2 = R1*math.sin(theta+beta)
            print(x1,y1,z1,end = ' ', file = fw)
            print(file = fw)
            print(x2,y2,z2,end = ' ', file = fw)
            print(file = fw)
fw.close()
# 制作ipt文件  谈谈吧，啥感受
position = []
f_in = open('luoxvan.xyz','r')
f_out = open('SBZ.itp','w')
fout = open('SBZ.pdb','w')
count=0
for line in f_in.readlines():
    count=count+1
f_in.close()
f_in = open('luoxvan.xyz','r')
p = int(count/n)
print('[moleculetype]',file=f_out)
print('; molname      nrexcl',file=f_out)
print('  MOL          1',file=f_out)
print(file=f_out)
print('[atoms]',file=f_out)
print('; id 	type 	resnr 	residu 	atom 	cgnr 	charge', file=f_out)
for i in range(1,int(count+1)):
    print('  ', i, ' P4  ', 1, '  MOL ', 'SBZ  ', i, '  0', end=' ', file=f_out)
    print(file=f_out)
print('[bonds]',file=f_out)
print(';  i  j 	funct 	length 	force.c.',file=f_out)
for i in range(count):
    line = f_in.readline()
    atom = line.split()
    ATOM = [float(x) for x in atom]
    position.append(ATOM)
for i in range(n-1):
    start = int(i*p)
    end = int((i+1)*p)
    for j in range(start,end):
        X = position[int(j+p)][0] - position[j][0]
        Y = position[int(j+p)][1] - position[j][1]
        Z = position[int(j+p)][2] - position[j][2]
        R = math.sqrt(X**2+Y**2+Z**2)
        print('  ', j+1,'', j+1+p, '  1 ', round(R/10,4), ' 1000', end=' ', file=f_out)
        print(file=f_out)
        for k in range(int(j+1),end):
            x = position[k][0] - position[j][0]
            y = position[k][1] - position[j][1]
            z = position[k][2] - position[j][2]
            r = math.sqrt(x**2+y**2+z**2)
            if r < 8:
                print('  ', j+1,'', k+1, '  1 ', round(r/10,4), ' 125000', end=' ', file=f_out)
                print(file=f_out)
for i in range(count):
    print("ATOM{0:>7d}{1:>3s}{2:>6s}{3:>6d}{4:>12.3f}{5:>8.3f}{6:>8.3f}{7:>6.2f}{8:>6.2f}".format(i+1,'P4','MOL',1,position[i][0],position[i][1],position[i][2],1.00,0.00), file=fout)
print("TER{0:>8d}".format(count+1), file=fout)
f_in.close()
f_out.close()
fout.close()