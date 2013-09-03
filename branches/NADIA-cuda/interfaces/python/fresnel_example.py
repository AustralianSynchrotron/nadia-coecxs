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
    proj.setAlgorithm(Algorithms.HIO)
    proj.initialiseEstimate(0)

    # run the reconstruction
    for i in range(1,100):
        result=object_estimate.get2dMAG()
        result.write_ppm("erfcdi_example_iter_%i.ppm" % i)
        temp_str=""

        #output the current estimate of the object
        tmp = PyDouble2D(nx,ny)
        tmp=proj.getIlluminationAtSample().get2dMAG()
        result=object_estimate.get2dPHASE()
        data.write_ppm("fcdi_ill_example_iter_%i.ppm" % i)

        #apply the shrinkwrap algorithm
        proj.applyShrinkwrap(2,0.1)  
        #iterate
        proj.iterate()
        print i
        print proj.getError()

    # output results
    print "done iterating"
    trans=PyComplex2D(nx,ny)
    trans=proj.getTransmissionFunction()
    print trans

    result=trans.get2dMAG()
    result.write_ppm("fcdi_example_trans_mag.ppm");

    result=trans.get2dPHASE()
    result.write_ppm("fcdi_example_trans_phase.ppm");

    print "done!"
except:
    print "something went wrong "
    print  traceback.print_exc(file=sys.stdout)
