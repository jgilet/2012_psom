#!/usr/bin/python

from scipy import *
from numpy import *
from pylab import *
import os


NJ = 96

NJ_bdy = 15

dy_bdy = 10e3
dy_int = 2e3

dyM = zeros(NJ+2)

y = linspace(0,1,NJ_bdy)
r = (tanh((-y+0.4)*5)+1.)/2

r = (4*(y-1))**4

r = (r-r.min())
r = r/r.max()

r = r*(dy_bdy-dy_int) + dy_int

dyM[0:NJ_bdy] = r
dyM[-NJ_bdy:] = r[::-1]
dyM[NJ_bdy:-NJ_bdy] = dy_int
dyM.astype('<f8').tofile('dyM.data')
print dyM

yc = cumsum(r_[-dyM[0]/2,(dyM[1:]+dyM[:-1])/2])
print yc

y = (yc-yc[0])/(yc[-1]-yc[0])
r_T = (tanh((-y+0.05)*20)+1.)/2.
r_T = r_T + r_T[::-1]

r_T = (r_T-r_T.min())/(r_T.max()-r_T.min())
print r_T
#r_T = r_T*0
r_T.astype('<f8').tofile('r.data')


print os.uname()
if 'Linux' in os.uname():
    subplot(2,1,1)
    plot(yc/1e3,r_T,'-o',yc/1e3,r_T[::-1],'r')
    ylabel('sponge r_T')
    xlabel('y km')
    subplot(2,1,2)
    plot(yc/1e3,dyM, '-o')
    ylabel('dyM')
    xlabel('y km')
    grid(True)
show()
