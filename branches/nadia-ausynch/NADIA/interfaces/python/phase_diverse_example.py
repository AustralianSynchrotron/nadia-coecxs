import os, sys, traceback
from pyNADIA.double2d import PyDouble2D
from pyNADIA.complex2d import PyComplex2D
from pyNADIA.fresnelcdi import PyFresnelCDI
from pyNADIA.phasediversecdi import PyPhaseDiverseCDI, CROSSCORRELATION,MINIMUMERROR
from pyNADIA.transmissionconstraint import PyTransmissionConstraint
from utils import Algorithms, MAG
from math import sqrt



try:
    nx=1024
    ny=1024
    pd=PyPhaseDiverseCDI()
    numframes =7
    tc1=PyTransmissionConstraint()
    estimatelist=[]
    cdiobjs=[]
        #initialise some arrays to hold data.
        
    
      #set the experimental parameters (all are in meters)
    wavelength = 4.892e-10 #wavelength
    fz = 16.353e-3 #zone plate to focal distance
    ps = 13.5e-6 #pixel size
    print 'running phase diverse example'
    os.environ["LD_RUN_PATH"]='../../lib'
    max_iterations=15
    
    beam = PyDouble2D(nx,ny)
    beam_fraction=0.5
    i0=(nx-1)/2.0
    j0=(ny-1)/2.0
    
    for i in range(0,nx):
        for j in range(0,ny):
            if sqrt((i-i0)*(i-i0)+(j-j0)*(j-j0)) < beam_fraction*sqrt(j0*j0+i0*i0):
                beam.set(i,j,100)
            else:
                beam.set(i,j,0)
    
    #zone plate to detector distances
    fd = [0.909513, 0.909388, 0.909263, 0.909213, 0.909088, 0.909088, 0.909088]

  #sample to detector distances
    fs = [18.513e-3, 18.388e-3, 18.263e-3, 18.213e-3, 18.088e-3, 18.088e-3, 18.088e-3]
  
  #white-field normalisation
    norm = [0.984729833,0.97700431,0.986270638,0.967487825, 0.980945916, 0.97628279, 0.963066039 ]

  
  #reconstructed white-field file names
    wf_string = ["../../examples/image_files/phasediverse_data/wf_A_1024.cplx",
                   "../../examples/image_files/phasediverse_data/wf_B_1024.cplx",
                   "../../examples/image_files/phasediverse_data/wf_C_1024.cplx",
                   "../../examples/image_files/phasediverse_data/wf_D_1024.cplx",
                   "../../examples/image_files/phasediverse_data/wf_E_1024.cplx",
                   "../../examples/image_files/phasediverse_data/wf_F_1024.cplx",
                   "../../examples/image_files/phasediverse_data/wf_G_1024.cplx"]

  #data file names
    diff_string = ["../../examples/image_files/phasediverse_data/A.dbin",
                     "../../examples/image_files/phasediverse_data/B.dbin",
                     "../../examples/image_files/phasediverse_data/C.dbin",
                     "../../examples/image_files/phasediverse_data/D.dbin",
                     "../../examples/image_files/phasediverse_data/E.dbin",
                     "../../examples/image_files/phasediverse_data/F.dbin",
                     "../../examples/image_files/phasediverse_data/G.dbin"]

  #these are the correct coordinates
    #x_pos = [0,26,6,28,-14,34,38]
   # y_pos = [-150,-149,-131,-127,-134,-95,211]

  #these are the coordinates smeared a bit
    x_pos = [0,23,3,35,30,37,60]
    y_pos = [-150,-145,-134,-123,-150,-93,210]
    
    wf = PyComplex2D(nx,ny)
    diffraction = PyDouble2D(nx,ny)
    
    for j in range(0,numframes):
        object_estimate=PyComplex2D(nx,ny)
        estimatelist.append(object_estimate)
        #read the white-field from file
        #wf = PyComplex2D(nx,ny)
        wf.read_cplx(wf_string[j])

        #read the data from file
        #diffraction = PyDouble2D(nx,ny)
        diffraction.read_dbin(diff_string[j],nx,ny)

        #set-up the reconstruction for a single frame
        proj=PyFresnelCDI(object_estimate,wf,wavelength,fd[j]-fz,fs[j]-fz,ps,norm[j])
        
    
        #set the support and intensity and initialise
        proj.setIntensity(diffraction)
        proj.setSupport(beam)
        proj.initialiseEstimate()
    
        #add a complex constraint
        proj.setComplexConstraint(tc1)
       
        
        #New part.. Add the FresnelCDI to the PhaseDiverseCDI.
        pd.addNewPosition(proj, x_pos[j], y_pos[j])
        #return proj
        cdiobjs.append(proj)
  
    
    #initialise the phase diverse transmission function
    pd.initialiseEstimate()
    
    #now run the reconstruction
    for i in range(0,max_iterations):
        #do some position alignment
        if i==0:
            pd.adjustPositions(CROSSCORRELATION)
        if i==10:
            pd.adjustPositions(MINIMUMERROR)
    #iteration!
        pd.iterate()

    theobject = PyComplex2D(nx,ny)
    theobject=pd.getTransmission()
    x,y=pd.getFinalXPosition(1), pd.getFinalYPosition(1)
    print x,y
    
    result=theobject.get2dMAG()
    
    #write the magnitude and phase of it to an image file

    result.write_image("object_mag_recovered.tiff",False,0.4,1.0)

    result=theobject.get2dPHASE()
    result.write_image("object_phase_recovered.tiff",False, -1.2,0)
    
    #save the transmission function in case we want to use it later.
    theobject.write_cplx("trans.cplx")  

except:
    print 'there was an error'
    print  traceback.print_exc(file=sys.stdout)
    
