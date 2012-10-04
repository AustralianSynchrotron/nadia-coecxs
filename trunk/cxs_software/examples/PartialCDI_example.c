// Copyright 2011 Nadia Davidson 
// for The ARC Centre of Excellence in Coherent X-ray Science. 
//
// This program is distributed under the GNU General Public License. 
// We also ask that you cite this software in publications where you made 
// use of it for any part of the data analysis.

/**
 * @file PlanarCDI_example.c
 *
 * \a PlanarCDI_example.c This example reconstructs some planar
 * diffraction data (Lachie's data). The shrinkwrap algorithm is used
 * to improve the reconstruction. A combination of HIO and the
 * error-reduction algorithm are used.
 *
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <fftw3.h>
#include <cstdlib> 
#include "io.h"
#include "Complex_2D.h"
#include "PartialCDI.h"
#include "Double_2D.h"
//#include <gperftools/profiler.h>


//#include "google/profiler.h"

//using namespace std;

int main(void){

 // ProfilerStart("partoutput.pprof");

  //Define some constants which will be used in the code.

  //the data file name
  string data_file_name = "image_files/part_data.dbin";
			  //"part_sim_intensity.tiff";



  string support_file_name = "image_files/part_support.tiff";
			  //"image_files/planar_support.tiff";

  //the file with the initial guess
  //string initial_guess_name = "08280_800.tiff";
  const int cycles = 6;

  const int er_iterations1 = 50;
  //number of hybrid input-out iterations to perform.
  const int hio_iterations = 100;
  //number of error reduction iterations to perform after the HIO.
  const int er_iterations2 = 50;

  //output the current image ever "output_iterations"
  int output_iterations = 10;

  //apply the shrinkwrap algorithm every "shrinkwrap iterations"
  int shrinkwrap_iterations = 200;

  //the number of pixels in x and y
  int nx = 2048;
  int ny = 2048;

  //The number of legendre polynomials and the square root of the modes
  //The number of legendre polynomials must be greater than or equal to the 
  //number of modes
  int nleg = 32;
  int nmodes = 7;

  /**** get the diffraction data from file and read into an array *****/

  Double_2D data;
  read_image(data_file_name, data, nx, ny);  

  /*  ostringstream atemp_str ( ostringstream::out ) ;
      atemp_str << "08280_log.ppm";
      write_image(atemp_str.str(), data, true);
   */

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

  //create the planar CDI object which will be used to
  //perform the reconstuction.
  double beta = 0.9;

  //Coherence lengths x and y in m
  double lcx = 13.3e-6;
  double lcy = 40.0e-3;

  //Pixel size detector in m
  double psize_x=13.5e-6;
  double psize_y=13.5e-6;

  //Energy of the beam in eV
  double e_beam=1400.0;

  //Distance between detector and sample in metres
  double z_sd=1.4;

  PartialCDI partial(object_estimate, beta, lcx, lcy, psize_x, psize_y, e_beam, z_sd, 4, 0);

  // 10000.0, 10000.0, 1, 4, 0);

  //set the support and intensity
  partial.set_support(support,false);

  partial.set_intensity(data);

  partial.set_threshold(+0.4e-5);

  //set the algorithm to hybrid input-output
  partial.set_algorithm(ER);

  //Initialise the wave function
  partial.initialise_matrices(nleg, nmodes);

  Double_2D result(nx, ny);

  ostringstream temp_strsupp ( ostringstream::out ) ;
  temp_strsupp << "support_tmp.ppm";
  write_image(temp_strsupp.str(), support);
/*
  for(int i=0; i< nmodes*nmodes; i++){
    ostringstream temp_str0 ( ostringstream::out ) ;
    partial.get_mode(i).get_2d(MAG,result);
    temp_str0 << "mode"<<i<<".ppm";
    write_image(temp_str0.str(), result);
  }

  ostringstream temp_str0 ( ostringstream::out ) ;
  object_estimate.get_2d(MAG,result);
  temp_str0 << "modes.ppm";
  write_image(temp_str0.str(), result);
*/

  //Initialise the current object ESW with a random numbers
  //"0" is the se:ed to the random number generator
  partial.initialise_estimate(7);

  //  partial.set_fftw_type(FFTW_ESTIMATE);

  //make a 2D object. This will be used to output the 
  //image of the current estimate.

  //  Double_2D result(nx,ny);
  object_estimate.get_2d(MAG,result);

  /******* for fun, let's get the autocorrelation *****/

  //  Double_2D autoc(nx,ny);
  //  partial.get_intensity_autocorrelation(autoc);
  //  write_image("test_autocorrelation.ppm", autoc, true); //"true" means log scale


  //ProfilerStart("profile.txt");

  //partial.set_algorithm(HIO);

  /*** run the reconstruction ************/
  for(int a=0; a<cycles; a++){


    for(int i=0; i<er_iterations1; i++){

      cout << "iteration " << i << endl;

      //apply the set of partial CDI projections 
      partial.iterate(); 
      cout << "Current error is "<<partial.get_error()<<endl;


      //every "output_iterations" 
      //output the current estimate of the object
      if(i%output_iterations==0){


	ostringstream temp_str ( ostringstream::out ) ;
	object_estimate.get_2d(MAG,result);
	temp_str << "part_example_iteration_" << i+a*(hio_iterations+er_iterations1+er_iterations2) << ".ppm";
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
	partial.apply_shrinkwrap(1.5,0.1);

    }

    //now change to the error reduction algorithm 
    partial.set_algorithm(HIO);

    for(int i=er_iterations1; i<(hio_iterations+er_iterations1); i++){

      cout << "iteration " << i << endl;

      partial.iterate(); 

      cout << "Current error is "<<partial.get_error()<<endl;

      if(i%output_iterations==0){
	//output the current estimate of the object
	ostringstream temp_str ( ostringstream::out ) ;
	object_estimate.get_2d(MAG,result);
	temp_str << "part_example_iteration_" << i+a*(hio_iterations+er_iterations1+er_iterations2) << ".ppm";
	write_image(temp_str.str(), result);
	temp_str.clear();

	//apply the shrinkwrap algorithm
	//partial.apply_shrinkwrap(1.5,0.1);
      }
      if(i%shrinkwrap_iterations==(shrinkwrap_iterations-1))
	partial.apply_shrinkwrap(1.5,0.1);

    }

    partial.set_algorithm(ER);
    for(int i=er_iterations1+hio_iterations; i<(hio_iterations+er_iterations1+er_iterations2); i++){

      cout << "iteration " << i << endl;

      partial.iterate();

      cout << "Current error is "<<partial.get_error()<<endl;

      if(i%output_iterations==0){
	//output the current estimate of the object
	ostringstream temp_str ( ostringstream::out ) ;
	object_estimate.get_2d(MAG,result);
	temp_str << "part_example_iteration_" << i+a*(hio_iterations+er_iterations1+er_iterations2) << ".ppm";
	write_image(temp_str.str(), result);
	temp_str.clear();

	//apply the shrinkwrap algorithm
	//partial.apply_shrinkwrap(1.5,0.1);
      }
      if(i%shrinkwrap_iterations==(shrinkwrap_iterations-1))
	partial.apply_shrinkwrap(1.5,0.1);

    }
  }

  //And we are done. "object_estimate" contained the final estimate of
  //the ESW.

  /** ignore the stuff below  
    Double_2D result2(nx,ny);

    double error=0;
    partial.get_best_result(0,error)->get_2d(MAG,result2);
    write_ppm("best_error.ppm", result2);
    cout << "Best error 0 is "<< error <<endl;

    partial.get_best_result(1,error)->get_2d(MAG,result2);
    write_ppm("best_error_1.ppm", result2);
    cout << "Best error 1 is "<< error <<endl;

    partial.get_best_result(2,error)->get_2d(MAG,result2);
    write_ppm("best_error_2.ppm", result2);
    cout << "Best error 2 is "<< error <<endl;

    partial.get_best_result(3,error)->get_2d(MAG,result2);
    write_ppm("best_error_3.ppm", result2);
    cout << "Best error 3 is "<< error <<endl; **/

//  ProfilerStop();
  write_cplx("PCDI_trans.cplx", object_estimate);

  return 0;
}

