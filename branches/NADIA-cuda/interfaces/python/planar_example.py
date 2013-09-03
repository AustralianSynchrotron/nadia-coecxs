import os, sys, traceback
from pyNADIA.double2d import PyDouble2D
from pyNADIA.complex2d import PyComplex2D
from pyNADIA.planarcdi import PyPlanarCDI
from utils import Algorithms, MAG

# create the required objects
# read in datasets
# set up the fresnelcdi object
# iterate
# write out the result
os.environ["LD_RUN_PATH"]='../../lib'
try:
    data_file_name = "planar_data.tif"

  #the file which provides the support (pixels with the value 0
  #are considered as outside the object)
    support_file_name = "planar_support.tiff"

  #number of cycles of ER and HIO to repeat
    cycles=2

  #number of error reduction iterations to perform before the HIO.
    er_iterations1 = 50

  #number of hybrid input-out iterations to perform.
    hio_iterations = 100

  #number of error reduction iterations to perform after the HIO.
    er_iterations2 = 50

  #output the current image every "output_iterations"
    output_iterations = 10

  #apply the shrinkwrap algorithm every "shrinkwrap iterations"
    shrinkwrap_iterations = 50

  #the number of pixels in x and y
    nx = 1024
    ny = 1024

  #**** get the diffraction data from file and read into an array *****/
    data=PyDouble2D(nx,ny)
    data.read_tiff(data_file_name)  

  #****** get the support from a file and read it into an array *****/

  #Every pixel with a zero value is intepreted as being outside
  #the support

    support=PyDouble2D(nx,ny)
    support.read_tiff(support_file_name)

  #check that the support image has the same dimensions
  #as the data.
    if support.getSizeX()!=nx or support.getSizeY()!=ny:
        print "support file has the wrong dimensions"
        exit
  

  #/*******  set up the reconstuction *********************/

#Create a complex 2D field which will hold the result of
#the reconstruction.
    object_estimate=PyComplex2D(nx,ny)

  #create the planar CDI object which will be used to
  #perform the reconstuction.
    planar=PyPlanarCDI(object_estimate, 4)

  #//set the support and intensity
    planar.setSupport(support,False)
    planar.setIntensity(data)

  #//Initialise the current object ESW with a random numbers
  #//"0" is the seed to the random number generator
    planar.initialiseEstimate(0)

  #//planar.set_fftw_type(FFTW_ESTIMATE);
  #//read_cplx("Planar_trans.cplx", object_estimate);

  #//make a 2D object. This will be used to output the 
  #//image of the current estimate.
  #Double_2D result(nx,ny);
    result=PyDouble2D(nx,ny)
    result=object_estimate.get2dMAG()

  #/******* for fun, let's get the autocorrelation *****/
    autoc = PyDouble2D(nx,ny)
    autoc=planar.getIntensityAutocorrelation()
    autoc.write_ppm("test_autocorrelation.ppm", True)# //"true" means log scale

  #//set the algorithm to Error Reduction
    planar.setAlgorithm(Algorithms.ER);

  #//  ProfilerStart("profile");


  #/*** run the reconstruction ************/
    for a in range(0,cycles):
      for i in range(0,er_iterations1):
          print "iteration %d" % (i+a*(hio_iterations+er_iterations1+er_iterations2))
           #apply the set of planar CDI projections 
          planar.iterate()
          print "Current error is %f" % planar.getError()
     
          #every "output_iterations" 
          #output the current estimate of the object
          if i%output_iterations==0:
              result=object_estimate.get2dMAG()
             
              result.write_tiff("planar_example_iteration_%d.tiff" % (i+a*(hio_iterations+er_iterations1+er_iterations2)))
              

	#//uncomment to output the estimated diffraction pattern
	#//planar.propagate_to_detector(object_estimate);
	#//object_estimate.get_2d(MAG_SQ,result);
	#//temp_str << "diffraction.ppm";
	#//write_ppm(temp_str.str(), result, true);
	#//planar.propagate_from_detector(object_estimate);
	#//object_estimate.get_2d(MAG,result);

	#//apply the shrinkwrap algorithm
	#//1.5 is the gaussian width in pixels
	#//0.1 is the threshold (10% of the maximum pixel).

      
          if i%shrinkwrap_iterations==(shrinkwrap_iterations-1):
              planar.applyShrinkwrap(2.0,0.1)


    #now change to the error reduction algorithm 
      planar.setRelaxationParameter(0.5)
      planar.setAlgorithm(Algorithms.HIO)

      for i in range(er_iterations1,hio_iterations+er_iterations1):
          itnum=i+a*(hio_iterations+er_iterations1+er_iterations2)
          print "iteration %(itnum)d" % locals() 

          planar.iterate() 

          print "Current error is %f" % planar.getError()

          if i%output_iterations==0:
	#//output the current estimate of the object
	#ostringstream temp_str ( ostringstream::out ) ;
              result=object_estimate.get2dMAG()
              temp_str="planar_example_iteration_%d.tiff" % (i+a*(hio_iterations+er_iterations1+er_iterations2))
              result.write_tiff(temp_str)

	#//apply the shrinkwrap algorithm
	#//planar.apply_shrinkwrap(1.5,0.1);
      
          if i%shrinkwrap_iterations==(shrinkwrap_iterations-1):
              planar.applyShrinkwrap(2.0,0.1)

    

      planar.setAlgorithm(Algorithms.ER)

    for i in range(er_iterations1+hio_iterations, hio_iterations+er_iterations1+er_iterations2):

      print "iteration %d" %(i+a*(hio_iterations+er_iterations1+er_iterations2))

      planar.iterate();

      print "Current error is %f" % planar.getError()

      if i%output_iterations==0:
	#//output the current estimate of the object
	result = object_estimate.get2dMAG()
	temp_str ="planar_example_iteration_%d.tiff" %(i+a*(hio_iterations+er_iterations1+er_iterations2))
	result.write_tiff(temp_str)

	#//apply the shrinkwrap algorithm
	#//planar.apply_shrinkwrap(1.5,0.1);
      
      if i%shrinkwrap_iterations==(shrinkwrap_iterations-1):
	planar.applyShrinkwrap(2.0,0.1)

    
  


  #//And we are done. "object_estimate" contained the final estimate of
  #//the ESW.
    object_estimate.write_cplx("Planar_trans.cplx")
  

except:
    print "something bad happened"
    print  traceback.print_exc(file=sys.stdout)
