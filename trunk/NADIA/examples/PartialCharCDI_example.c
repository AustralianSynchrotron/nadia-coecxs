// Copyright 2013 Daniel Rodgers-Pryor for The ARC Centre of Excellence in 
// Coherent X-ray Science. This program is distributed under the GNU
// General Public License. We also ask that you cite this software in 
// publications where you made use of it for any part of the data 
// analysis.

/**
 * @file PartialCharCDI_example.c
 *
 * PartialCharCDI_example.c This example reconstructs some Partial
 * diffraction data (Diffractive imaging using partially coherent X rays).
 * This approach assumes a gaussian coherence function and can characterise,
 * the coherence of the beam as it reconstructs the image.
 * For more detail, see Applied Physics Letters 99 (2011) - 'Simultaneous sample and 
 * spatial coherence characterisation using diffractive imaging' [doi: 10.1063/1.3650265]
 * The shrinkwrap algorithm is used to improve the reconstruction. A 
 * combination of HIO and the error-reduction algorithm are used.
 *
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstring>
#include <stdlib.h>
#include <fftw3.h>
#include <cstdlib> 
#include <io.h>
#include <Complex_2D.h>
#include <PartialCharCDI.h>
#include <Double_2D.h>

//using namespace std;

int main(void){

 // ProfilerStart("partoutput.pprof");

  //Define some constants which will be used in the code.

  //the data file name
  //string data_file_name = "image_files/part_data.dbin";
  string data_file_name = "image_files/part_data.dbin";

  string support_file_name = "image_files/part_support.tiff";

  const int cycles = 2;
  //number of cycles of ER, HIO,ER to perform.

  const int er_iterations = 50;
  //number of hybrid input-out iterations to perform.
  const int hio_iterations = 50;
  //number of error reduction iterations to perform after the HIO.
  const int dm_iterations = 50;

  //output the current image ever "output_iterations"
  int output_iterations = 10;

  //apply the shrinkwrap algorithm every "shrinkwrap iterations"
  int shrinkwrap_iterations = 150;

  //the number of pixels in x and y
  int nx = 2048;
  int ny = 2048;
  
  //Pixel size detector in m
  double psize_x=13.5e-6;
  double psize_y=13.5e-6;
  
  //Energy of the beam in eV
  double e_beam=1400;

  //Distance between detector and sample in metres
  double z_sd=1.4;


  /**** get the diffraction data from file and read into an array *****/

  Double_2D data;
  read_image(data_file_name, data, nx, ny);

  /****** get the support from a file and read it into an array *****/

  //Every pixel with a zero value is intepreted as being outside
  //the support

  Double_2D support;
  read_image(support_file_name, support, nx, ny);

  //check that the support image has the same dimensions
  //as the data.
  if(support.get_size_x()!=nx || support.get_size_y()!=ny){
    cout << "support file has the wrong dimensions"  << endl;
    return(1);
  }

  /*******  set up the reconstuction *********************/

  //Create a complex 2D field which will hold the result of
  //the reconstruction.

  Complex_2D object_estimate(nx,ny);
  //read_image(initial_guess_name, object_estimate, nx, ny);

  PartialCharCDI partial(object_estimate, psize_x, psize_y, e_beam, z_sd, 4);

  
  
  //printf("%f %f\n", partial.get_x_coherence_length()*1e6, partial.get_y_coherence_length()*1e6);
  
  

  // 10000.0, 10000.0, 1, 4, 0);

  //set the support and intensity
  partial.set_support(support,false);

  partial.set_intensity(data);

  //set the algorithm to hybrid input-output
  partial.set_algorithm(ER);

  ostringstream temp_strsupp ( ostringstream::out ) ;
  temp_strsupp << "support_tmp.ppm";
  write_image(temp_strsupp.str(), support);

  //Initialise the current object ESW with a random numbers
  //"0" is the seed to the random number generator
  partial.initialise_estimate(0);

  //read_cplx("PCDI_trans.cplx", object_estimate);

  //make a 2D object. This will be used to output the 
  //image of the current estimate.
  Double_2D result(nx,ny);

  object_estimate.get_2d(MAG,result);

  //ProfilerStart("profile.txt");

  /*** run the reconstruction ************/
  for(int a=0; a<cycles; a++){

    for(int i=0; i<er_iterations; i++){

      cout << "iteration " << i+a*(hio_iterations+er_iterations+dm_iterations) << endl;

      //apply the set of partial CDI projections 
      partial.iterate(); 
      cout << "Current error is "<<partial.get_error()<<endl;

      //every "output_iterations" 
      //output the current estimate of the object
      if(i%output_iterations==0){

	ostringstream temp_str ( ostringstream::out ) ;
	object_estimate.get_2d(MAG,result);
	temp_str << "part_char_example_iteration_" << i+a*(hio_iterations+er_iterations+dm_iterations) << ".ppm";
	write_image(temp_str.str(), result);

	temp_str.clear();

	//uncomment to output the estimated diffraction pattern
	//partial.propagate_to_detector(object_estimate);
	//object_estimate.get_2d(MAG_SQ,result);
	//temp_str << "diffraction.ppm";
	//write_ppm(temp_str.str(), result, true);
	//partial.propagate_from_detector(object_estimate);
	//object_estimate.get_2d(MAG,result);

	//apply the shrinkwrap algorithm
	//1.5 is the gaussian width in pixels
	//0.1 is the threshold (10% of the maximum pixel).

      }
      if(i%shrinkwrap_iterations==(shrinkwrap_iterations-1))
	partial.apply_shrinkwrap(5.0,0.1);

    }

    //now change to the error reduction algorithm 
    partial.set_algorithm(HIO);

    for(int i=er_iterations; i<(hio_iterations+er_iterations); i++){

      cout << "iteration " << i+a*(hio_iterations+er_iterations+dm_iterations) << endl;

      partial.iterate(); 

      cout << "Current error is "<<partial.get_error()<<endl;

      if(i%output_iterations==0){
	//output the current estimate of the object
	ostringstream temp_str ( ostringstream::out ) ;
	object_estimate.get_2d(MAG,result);
	temp_str << "part_char_example_iteration_" << i+a*(hio_iterations+er_iterations+dm_iterations) << ".ppm";
	write_image(temp_str.str(), result);
	temp_str.clear();

	//apply the shrinkwrap algorithm
	//partial.apply_shrinkwrap(1.5,0.1);
      }
      if(i%shrinkwrap_iterations==(shrinkwrap_iterations-1))
	partial.apply_shrinkwrap(5.0,0.1);

    }

    partial.set_algorithm(ER);
    for(int i=er_iterations+hio_iterations; i<(hio_iterations+er_iterations+dm_iterations); i++){

      cout << "iteration " << i+a*(hio_iterations+er_iterations+dm_iterations) << endl;

      partial.iterate();

      cout << "Current error is "<<partial.get_error()<<endl;

      if(i%output_iterations==0){
	//output the current estimate of the object
	ostringstream temp_str ( ostringstream::out ) ;
	object_estimate.get_2d(MAG,result);
	temp_str << "part_char_example_iteration_" << i+a*(hio_iterations+er_iterations+dm_iterations) << ".ppm";
	write_image(temp_str.str(), result);
	temp_str.clear();

	//apply the shrinkwrap algorithm
	//partial.apply_shrinkwrap(1.5,0.1);
      }
      if(i%shrinkwrap_iterations==(shrinkwrap_iterations-1))
	partial.apply_shrinkwrap(5.0,0.1);

    }
  }

  //And we are done. "object_estimate" contained the final estimate of
  //the ESW.

//  ProfilerStop();
  write_cplx("PCCDI_trans.cplx", object_estimate);
  
  printf("\nFinal estimated coherence lengths:\n\tlx %fum\n\tly %fum\n", partial.get_x_coherence_length()*1e6, partial.get_y_coherence_length()*1e6);

  return 0;
}

