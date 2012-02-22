#!/usr/bin/python

from scipy import fromfile
from sys import argv
from pylab import *
#from enthought.mayavi.mlab import *

#fn = argv[1]
#a=fromfile(fn+'.bin')
a=fromfile('tmp-parti.bin')

print a.shape

NP=33
Nvar = 13

b=a.reshape(-1,NP,Nvar)
figure(1)
x,y=meshgrid(arange(49),arange(49))
#contour(sin(pi*x/48.)*sin(pi*y/48.))
for i in arange(NP):
    #aa=plot3d(b[:,i,0],b[:,i,1],b[:,i,2],colormap='hot',tube_radius=None)
    #aa=plot(b[0,i,0],b[0,i,1],'k.',markersize=0.5)
    aa=plot(b[:,i,1],b[:,i,2],'k.',markersize=0.5)
#axis('equal')
#xlabel('x')
#ylabel('y')
#title(r'$\psi = sin(\pi x/L_x) sin(\pi y/L_y)$')
#title('the chaotic mixer')
#savefig(fn+'.png')
show()
