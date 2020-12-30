# -*- coding: utf-8 -*-
"""
Created on Wed Dec 30 19:55:55 2020
分束器
@author: David Lyu
"""
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 26 19:41:57 2020
使用crank-nicolson方法求解BPM方程,波导的方向耦合器
@author: David Lyu
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 11:14:51 2020
2d WMF by PLC
@author: lenovo
"""
import numpy as np
from matplotlib import pyplot as plt


dx=0.1
dz=dx
wavelength=1.55
k0=2*np.pi/wavelength
nmax=1.45+1.45*0.0075/2
nmin=1.45-1.45*0.0075/2
#nmin=1.45-1.45*0.0075/2

beta=k0*nmax*0.9935#0.9925280199252802
betaa=0

lx=200
lz=800
d=3
dsourse=lx/13    #光束宽度
dwaveguide=lx/5        #波导宽度

phi=np.zeros((lx,lz),dtype=complex)
psi=np.zeros((lx,lz),dtype=complex)

a=np.zeros((lx),dtype=complex)
b=np.zeros((lx),dtype=complex)
c=np.zeros((lx),dtype=complex)

p=np.zeros((lx-2,lx-2),dtype=complex)
q=np.zeros((lx-2,lx-2),dtype=complex)
aa=np.zeros((lx-2,lx-2),dtype=complex)
aaa=np.zeros((lx-2,lx-2),dtype=complex)

n=np.zeros((lx,lz))
nn=np.zeros((lx,lz))
q0=np.zeros((lx,lz))
n[:,:]=nmax

h=complex(0,0)
j=complex(0,1)
"""
def sinn(x):
    if x>1.1:
        return 0.9
    else:
        return 0.9
"""
    
alpha=0.05

#设定折射率的区域
theta=0.1
d=200
for x in range(0,lx):
    for z in range(0,lz):
        if (lx/2-dwaveguide/2<=x<=lx/2+dwaveguide/2)and  (z<=d):
            n[x,z]=nmax
        elif (-theta*(z-d)-dwaveguide/2<=x-lx/2<=-theta*(z-d)+dwaveguide/2)and(z>=d):
            n[x,z]=nmax
        elif (theta*(z-d)-dwaveguide/2<=x-lx/2<=theta*(z-d)+dwaveguide/2)and(z>=d):
            n[x,z]=nmax
        else:
            n[x,z]=nmin

for t in range(0,1):
    #设定光源
    for x in range(0,lx-1):
        psi[x,0]=1*np.exp(-(x-lx/2)**2/(dsourse)**2)
        if psi[x,0]<0.001:
            psi[x,0]=0
        phi[x,0]=1*np.exp(-(x-lx*3/4)**2/(dsourse)**2)
        if phi[x,0]<0.001:
               phi[x,0]=0
    nn=np.fliplr(n)
    #计算P,Q两个矩阵
    for z in range(0,lz-1):
        b[:]=2*beta+j*dz/dx/dx-j*k0*k0*(n[:,z+1]**2-(beta/k0)**2)/2
        c[:]=2*beta-j*dz/dx/dx+j*k0*k0*(n[:,z]**2-(beta/k0)**2)/2
        a[:]=j*dz/dx/dx/2
        p[:,:]=np.diag(b[1:lx-1])
        aa[:,:]=np.diag(a[1:lx-1])
        aa[:,:]=np.vstack((aa[1:,:],np.zeros((lx-2),dtype=complex)))
        
        aaa[:,:]=np.diag(a[1:lx-1])
        aaa[:,:]=np.hstack((aaa[:,1:],np.zeros((lx-2,1),dtype=complex)))
        p=p-aa-aaa
        q[:,:]=np.diag(c[1:lx-1])
        q=q+aa+aaa
        psi[1:lx-1,z+1]=(np.linalg.inv(p)).dot(q.dot(psi[1:lx-1,z]))
    phi=np.fliplr(phi)
"""
    for z in range(0,lz):
        for x in range(0,lx):
            h=psi[x,z]*np.conjugate(phi[x,z])
            q[x,z]=h.imag
            n[x,z]=n[x,z]-alpha*q[x,z]
            if (n[x,z]>1.5):
                n[x,z]=1.5
            elif(n[x,z]<1):
                n[x,z]=1
"""
abspsi=abs(psi)
plt.figure(figsize = (16, 9))
"""
plt.subplot(2,2,1)
plt.imshow(np.real(psi),cmap='gray',vmin=-1, vmax=1)
plt.subplot(2,2,2)
plt.imshow(abs(phi),cmap='gray',vmin=0, vmax=1)
#plt.plot(range(0,lx),abs(psi[:,lz-1]))
"""
plt.subplot(2,1,1)
plt.imshow(abs(psi),cmap='gray',vmin=0, vmax=1)
plt.subplot(2,1,2)
plt.imshow(n,cmap='gray')
    

