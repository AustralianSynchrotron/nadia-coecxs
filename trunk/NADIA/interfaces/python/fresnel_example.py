import os, sys, traceback
from pyNADIA.double2d import PyDouble2D
from pyNADIA.complex2d import PyComplex2D
from pyNADIA.fresnelcdi import PyFresnelCDI
from utils import Algorithms, MAG

# create the required objects
# read in datasets
# set up the fresnelcdi object
# iterate
# write out the result
os.environ["LD_RUN_PATH"]='../../lib'
try:
    # read all our files into data objects
    nx =1024
    ny =1024
    data = PyDouble2D()
    data.read_dbin('../../examples/image_files/FCDI_data.dbin',nx,ny)
    support = PyDouble2D(nx,ny)
    support.read_tiff("../../examples/image_files/FCDI_support.tiff")

    wf = PyComplex2D(nx,ny)
    # should include checking that this actually exists...
    wf.read_cplx('../../examples/wf_recovered.cplx')

    # set up the reconstruction
    result = PyDouble2D(nx,ny)
    object_estimate=PyComplex2D(nx,ny)
    wavelength=4.892e-10
    ftd=0.909513 - 16.353e-3 #focal-to-detector
    fsd=18.513e-3 - 16.353e-3 #focal-to-sample
    pixel_size=13.5e-6
    normalisation=0.984729833
    proj = PyFresnelCDI(object_estimate,wf,wavelength,ftd,fsd,pixel_size,normalisation)
    proj.setSupport(support)
    proj.setIntensity(data)
    proj.setAlgorithm(Algorithms.ER)
    proj.initialiseEstimate(0)

    # run the reconstruction
    for i in range(0,20):
        print "iteration %d" % i

        #apply the iterations  
        proj.iterate(); 
        print "Error: %f" % proj.getError()

        if i%5==0:
            #output the current estimate of the object
            result=object_estimate.get2dMAG()
            temp_str="fcdi_example_iter_%d.ppm" % i
            result.write_ppm(temp_str)
    
            #apply the shrinkwrap algorithm
            proj.applyShrinkwrap(2,0.1);

       
       

    # output results
    print "done iterating"
    trans=PyComplex2D(nx,ny)
    trans=proj.getTransmissionFunction()

   # result=trans.get2dMAG()
    result=proj.getTransmissionFunction().get2dMAG()

    result.write_ppm("fcdi_example_trans_mag.ppm",log_scale=True)

    result=trans.get2dPHASE()
    result.write_ppm("fcdi_example_trans_phase.ppm")

    print "done!"
except:
    print "something went wrong "
    print  traceback.print_exc(file=sys.stdout)
