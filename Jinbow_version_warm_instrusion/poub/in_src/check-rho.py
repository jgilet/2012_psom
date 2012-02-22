#!/usr/bin/python
import scipy
#from scipy import *
from pylab import *
from sys import argv
g=loadtxt('zgrid.out')
#a=fromfile('op.rho.0000005000.bin',float32)
#a=fromfile('op.vor.0000005000.bin',float32)
y=fromfile('dyM.data','f8')
y=cumsum(y)/1e3 -y[0]/2e3


a=fromfile(argv[1],float32)
print a.shape
NI,NJ,NK = 48, 32,  16
NI,NJ,NK = 96, 96,  32
NI,NJ,NK = 48, 96,  32
NI,NJ,NK = 96, 96,  64
NI,NJ,NK = 96, 192, 64
NI,NJ,NK = 96, 96,  24
NI,NJ,NK = 48, 48,  32
NI,NJ,NK = 96, 192, 24
NI,NJ,NK = 96, 192, 32
NI,NJ,NK = 48, 48,  24 
NI,NJ,NK = 24, 48,  24 
NI,NJ,NK = 24, 24,  24 
y=arange(NJ+2)
z1 = g[0:NK+2,1]
yz,zy = meshgrid(y,z1)
xz,zx = meshgrid(arange(NI+2),z1)
xy,yx = meshgrid(arange(NI+2),y)

b=a.reshape(NK+2,NJ+2,NI+2)
figure(1)
#contourf(xy,yx,b[-5,:,:])
pcolor(b[-5,:,:])
colorbar()
title('xy')
xlabel('x')
ylabel('y')
figure(2)
contourf(yz,zy,b[:,:,NI/2])
colorbar()
xlabel('y')
ylabel('z')
title('yz')
figure(3)
subplot(2,2,1)
plot(b[:,NJ/2,NI/2].flatten(),zy[:,0])
subplot(2,2,2)
plot(y,b[-8,:,NI/2].flatten(),'-*')
subplot(2,2,3)
plot(-10./1025.*diff(b[:,NJ/2,NI/2].flatten())/diff(z1),z1[1::])
xlabel('y')
ylabel('z')
title('yz')
show()
#scipy.imsave("out.png",b)
