# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
#import numpy as np
import glob
import os
import math
from itertools import groupby

os.getcwd()
os.chdir(r"C:\Users\Wang\Desktop\T380")

N = [] # 分子数
S = []   # 盒子边界
T = 0    # 体系数
A = 30  # 分子数
#读取文件及编辑数据
INCELL = sorted(glob.glob('*.gro'))
for f in  INCELL:
    [dirname, filename] = os.path.split(f) #文件名分割
    f_out = open(filename.split('.')[0]+'.xyz','w')
    f_in = open(filename.split('.')[0]+'.gro','r')
    line =  f_in.readline()
    line =  f_in.readline()
    length = int(line) #总原子数
    N.append(int(length/A)) # 总分子数
    for i in range(length):
        if i%A == 0:
            print(int(i/A), end = ' ', file = f_out)
        line =  f_in.readline()
        atom = [j.strip() for j in [line[:10], line[10:15], line[15:22], line[22:30], line[30:38], line[38:45]]]
        for k in range(3,6):
            print (atom[k],end = ' ', file = f_out)
        if (i+1)%A == 0:
            print (file = f_out)
    line = f_in.readline()
    box = [h.strip() for h in [line[:11], line[11:22], line[15:22], line[22:32], line[32:42], line[42:52], line[52:62], line[62:72], line[72:82], line[82:92]]]
    S.append(box[0])
f_out.close()
f_in.close()
#去周期性
HEIGOU = sorted(glob.glob('*.xyz'))
for r in HEIGOU:
    [dirname, filename] = os.path.split(r) #文件名分割
    fin = open(filename.split('.')[0]+'.xyz','r')
    fout = open(filename.split('.')[0]+'.txt','w')
    size = float(S[T])
    half = size/2.0
    n = N[T]
    T = T+1
    for l in range(n):
        line = fin.readline()   
        position = line.split()      
        for m in range(A-1):
            x = 3*m+1
            y = 3*m+2
            z = 3*m+3
            X = x+3
            Y = y+3
            Z = z+3
            lx = abs(float(position[x])-float(position[X]))
            Dx = float(position[x])-float(position[X])
            ly = abs(float(position[y])-float(position[Y]))
            Dy = float(position[y])-float(position[Y])
            lz = abs(float(position[z])-float(position[Z]))
            Dz = float(position[z])-float(position[Z])
            if lx > half and Dx > 0:
                position[X] = float(position[X]) + size
            if lx > half and Dx < 0:
                position[X] = float(position[X]) - size    
            if ly > half and Dy > 0:
                position[Y] = float(position[Y]) + size
            if ly > half and Dy < 0:
                position[Y] = float(position[Y]) - size 
            if lz > half and Dz > 0:
                position[Z] = float(position[Z]) + size
            if lz > half and Dz < 0:
                position[Z] = float(position[Z]) - size 
        for a in range(A*3+1):
            if a == 0:
                p = int (position[a])
                print(p, end = ' ', file = fout)
            else: 
                p = float (position[a])
                print ("%.3f" % p, end = ' ', file = fout)
        print (file = fout)
fout.close()
fin.close()
# 取三重心及末端距
T = 0  # 计数器归零 
MED = 0.000 # 末端距
Gx = 0.000 # 重心x
Gy = 0.000 # 重心y
Gz = 0.000 # 重心z
g1 = 0.000 # 重心x
g2 = 0.000 # 重心y
g3 = 0.000 # 重心z
gi = 0.000 # 重心x
gj = 0.000 # 重心y
gk = 0.000 # 重心z
PINGTOU = sorted(glob.glob('*.txt'))    
for z in PINGTOU:
    [dirname, filename] = os.path.split(z) #文件名分割
    fi = open(filename.split('.')[0]+'.txt','r')
    fo = open(filename.split('.')[0]+'.data','w')
    n = N[T]
    T = T+1
    for l in range(n):
        line = fi.readline()   
        position = line.split()
        print (position[0], end = ' ', file = fo)
        x = 3*A-2
        y = 3*A-1
        z = 3*A
        for m in range(0,int(A/3)):
            x = 3*m+1
            y = 3*m+2
            z = 3*m+3
            g1 = (g1 + float(position[x]))
            g2 = (g2 + float(position[y]))
            g3 = (g3 + float(position[z]))
            x1 = g1/A*3.0
            x2 = g2/A*3.0
            x3 = g3/A*3.0
        print ('%.3f' % x1, end = ' ', file = fo)
        print ('%.3f' % x2, end = ' ', file = fo)
        print ('%.3f' % x3, end = ' ', file = fo)
        g1 = 0.000 # 重心x归零
        g2 = 0.000 # 重心y归零
        g3 = 0.000 # 重心z归零
        for m in range(A):
            x = 3*m+1
            y = 3*m+2
            z = 3*m+3
            Gx = (Gx + float(position[x]))
            Gy = (Gy + float(position[y]))
            Gz = (Gz + float(position[z]))
            gx = Gx/A
            gy = Gy/A
            gz = Gz/A
        print ('%.3f' % gx, end = ' ', file = fo)
        print ('%.3f' % gy, end = ' ', file = fo)
        print ('%.3f' % gz, end = ' ', file = fo)
        Gx = 0.000 # 重心x归零
        Gy = 0.000 # 重心y归零
        Gz = 0.000 # 重心z归零
        for m in range(int(A/3*2),A):
            x = 3*m+1
            y = 3*m+2
            z = 3*m+3
            gi = (gi + float(position[x]))
            gj = (gj + float(position[y]))
            gk = (gk + float(position[z]))
            xi = gi/A*3.0
            xj = gj/A*3.0
            xk = gk/A*3.0
        print ('%.3f' % xi, end = ' ', file = fo)
        print ('%.3f' % xj, end = ' ', file = fo)
        print ('%.3f' % xk, end = ' ', file = fo)
        gi = 0.000 # 重心x归零
        gj = 0.000 # 重心y归零
        gk = 0.000 # 重心z归零
        print (file = fo)
fo.close()
fi.close()
T = 0
PANGDIDI = sorted(glob.glob('*.data')) 
for d in PANGDIDI:
    [dirname, filename] = os.path.split(d) #文件名分割
    f1 = open(filename.split('.')[0]+'.data','r')
    f2 = open(filename.split('.')[0]+'.txt','w')
    n = N[T]
    T = T+1
    for i in range(n):
        line = f1.readline()
        angle = line.split()
        x = float(angle[7])-float(angle[3])
        y = float(angle[8])-float(angle[4])
        z = float(angle[9])-float(angle[5])
        Angle = z/math.sqrt(x**2+y**2+z**2)
        jiaodu = math.degrees(math.acos(Angle))
        if jiaodu > 90:
            jiaodu = 180 - jiaodu
        print(jiaodu, end = '\n', file = f2)
f1.close()
f2.close()
sb = []
DAXIONG = sorted(glob.glob('*.txt')) 
for d in DAXIONG:
    [dirname, filename] = os.path.split(d) #文件名分割
    fb = open(filename.split('.')[0]+'.txt','r')
    line = fb.readlines()
    jiao = [float(x) for x in line]
    sb.append(jiao)
fb.close()    
fd = open('jun.pdb','w')
FRAME = len(sb)
Nmole = len(sb[0])
for i in range(Nmole):
    jj = 0
    for j in range(FRAME):
        jj = sb[j][i] + jj
    jj = jj/FRAME
    print(jj, end = '\n', file = fd)
fd.close()
f11 = open('jun.pdb','r')        
f22 = open('fenbu.pdb','w')
line = f11.readlines()
jiaofenbu = [float(x) for x in line] 
for i, j in groupby(sorted(jiaofenbu), key=lambda x: x//5):
    print('{:.2f} {:.0f}'.format(i*5, len(list(j))), file = f22)
f11.close()
f22.close()
exit()