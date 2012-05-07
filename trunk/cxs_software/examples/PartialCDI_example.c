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

//#include "google/profiler.h"

using namespace std;

int main(void){

  //Define some constants which will be used in the code.

  //the data file name
  string data_file_name = "part_sim_intensity.tiff";//*/"image_files/planar_data.tif";//08280.hdf";
  //string data_file_name = "08280.hdf";

    //"part_sim_intensity.tiff";
    //"image_files/planar_data.tif";
  //part_sim_intensity.tiff";//08280.hdf";

  //the file which provides the support (pixels with the value 0
  //are considered as outside the object)
  string support_file_name = /*"support_02141.tiff";//*/"image_files/planar_support.tiff";

  //the file with the initial guess
  //string initial_guess_name = "08280_800.tiff";

  //number of hybrid input-out iterations to perform.
  const int hio_iterations = 51;
  
  //number of error reduction iterations to perform after the HIO.
  const int er_iterations = 50;

  //output the current image ever "output_iterations"
  int output_iterations = 10;

  //apply the shrinkwrap algorithm every "shrinkwrap iterations"
  int shrinkwrap_iterations = 50;

  //the number of pixels in x and y
  int nx =1024;
  int ny = 1024;

  //the number of legendre polynomials and the square root of the modes
  int nleg = 12;

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
  PartialCDI partial(object_estimate, 0.9, 1.5e+0, 1.5e+0, 4.0, 4.0, 4, 0);
  // 10000.0, 10000.0, 1, 4, 0);

  //set the support and intensity
  partial.set_support(support,false);

  partial.set_intensity(data);

//  partial.set_threshold(1.0e-5);

  //set the algorithm to hybrid input-output
  partial.set_algorithm(HIO);

  //Initialise the wave function
  partial.initialise_matrices(nleg, nmodes);

  Double_2D result(nx, ny);

  for(int i=0; i< nmodes*nmodes; i++){
    ostringstream temp_str0 ( ostringstream::out ) ;
    partial.get_mode(i).get_2d(REAL,result);
    temp_str0 << "mode"<<i<<".ppm";
    write_image(temp_str0.str(), result);
  }

  ostringstream temp_str0 ( ostringstream::out ) ;
  object_estimate.get_2d(MAG,result);
  temp_str0 << "modes.ppm";
  write_image(temp_str0.str(), result, false);


  //Initialise the current object ESW with a random numbers
  //"0" is the se:ed to the random number generator
  partial.initialise_estimate(0);

  //  partial.set_fftw_type(FFTW_ESTIMATE);

  //make a 2D object. This will be used to output the 
  //image of the current estimate.

  //  Double_2D result(nx,ny);
  object_estimate.get_2d(MAG,result);

  /******* for fun, let's get the autocorrelation *****/

  //Double_2D autoc(nx,ny);
  //partial.get_intensity_autocorrelation(autoc);
  //write_image("test_autocorrelation.ppm", autoc, true); //"true" means log scale


  //  ProfilerStart("profile");

  //partial.set_algorithm(HIO);

  /*** run the reconstruction ************/

  for(int i=0; i<hio_iterations; i++){

    cout << "iteration " << i << endl;

    //apply the set of partial CDI projections 
    partial.iterate(); 
    cout << "Current error is "<<partial.get_error()<<endl;


    //every "output_iterations" 
    //output the current estimate of the object
    if(i%output_iterations==0){


      ostringstream temp_str ( ostringstream::out ) ;
      object_estimate.get_2d(MAG,result);
      temp_str << "part_example_iteration_" << i << ".ppm";
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
  partial.set_algorithm(ER);

  for(int i=hio_iterations; i<(hio_iterations+er_iterations+1); i++){

    cout << "iteration " << i << endl;

    partial.iterate(); 

    cout << "Current error is "<<partial.get_error()<<endl;

    if(i%output_iterations==0){
      //output the current estimate of the object
      ostringstream temp_str ( ostringstream::out ) ;
      object_estimate.get_2d(MAG,result);
      temp_str << "part_example_iteration_" << i << ".ppm";
      write_image(temp_str.str(), result);
      temp_str.clear();

      //apply the shrinkwrap algorithm
      //partial.apply_shrinkwrap(1.5,0.1);
    }
    if(i%shrinkwrap_iterations==(shrinkwrap_iterations-1))
      partial.apply_shrinkwrap(1.5,0.1);

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
