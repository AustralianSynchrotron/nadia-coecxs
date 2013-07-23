import os, sys, traceback
from pyNADIA.double2d import PyDouble2D
from pyNADIA.complex2d import PyComplex2D
from pyNADIA.fresnelcdi import PyFresnelCDI
from pyNADIA.phasediversecdi import PyPhaseDiverseCDI
from utils import Algorithms, MAG

# create the required objects
# read in datasets
# set up the fresnelcdi object
# iterate
# write out the result
os.environ["LD_RUN_PATH"]='../../lib'
try:
    nx=1024
    ny=1024
    max_iterations=15

    # make the support
    beam = PyDouble2D(nx,ny)
    beam_fraction=0.5
    i0=(nx-1)/2.0
    j0=(ny-1)/2.0

    # there's probably a clever pythonic way to do this but it matches the C example
    for i in range(0,nx):
        for j in range(0,ny):
            if sqrt((i-i0)*(i-i0)+(j-j0)*(j-j0)) < beam_fraction*sqrt(j0*j0+i0*i0):
                beam.set(i,j,100)
            else:
                beam.set(i,j,0)
    
    # make a new PhaseDiverseCDI object
    pd = PyPhaseDiverseCDI()
    
    numframes =7
    
    
except:
