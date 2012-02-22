#!/usr/bin/python

from scipy import fromfile
from pylab import *
from sys import argv
g=loadtxt('zgrid.out')
#a=fromfile('op.rho.0000005000.bin',float32)
#a=fromfile('op.vor.0000005000.bin',float32)
a=fromfile(argv[1],float32)
print a.shape
NI,NJ,NK = 48,32,16
NI,NJ,NK =48,48,24 
NI,NJ,NK =96,96,32
NI,NJ,NK =48,96,32
NI,NJ,NK =96,96,64
NI,NJ,NK = 96,192,64
NI,NJ,NK =96,96,24
NI,NJ,NK = 96,192,32
NI,NJ,NK =48,48,32
z1 = g[0:NK+2,1]
yz,zy = meshgrid(arange(NJ+2),z1)
xz,zx = meshgrid(arange(NI+2),z1)
xy,yx = meshgrid(arange(NI+2),arange(NJ+2))

b=a.reshape(NK+2,NJ+2)
figure(1)
contourf(b[1:-1,1:-1])
colorbar()
title('xy')
xlabel('x')
ylabel('y')
show()
