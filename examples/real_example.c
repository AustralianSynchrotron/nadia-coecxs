/**
 * @file real_example.c
 * @author  Nadia Davidson <nadiamd@unimelb.edu.au>
 *
 * @section DESCRIPTION
 *
 * This file provides an example of running the planar CDI 
 * reconstruction on real data (from the file test_dat.tif)
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
#include "PlanarCDI.h"
#include "shrinkwrap.h"
#include "Double_2D.h"
#include <google/profiler.h>

using namespace std;

int main(void){

  //Config c("my.in");

  ProfilerStart("profiler.prof");

  //Define some constants which will be used in the code.

  //the data file name
  string data_file_name = "image_files/test_dat.tif";

  //the approx. level of background noise in the image
  double noise_level = 30;

  //the file which provides the support (pixels with the value 0
  //are considered as outside the object)
  string support_file_name = "image_files/support_2.tiff";

  //number of hybrid input-out iterations to perform.
  const int hio_iterations = 500;
  
  //number of error reduction iterations to perform after the HIO.
  const int er_iterations = 200;

  //output the current image ever "output_iterations"
  int output_iterations = 40;

  /**** get the diffraction data from file and read into an array *****/

  Double_2D data;
  int status = read_tiff(data_file_name, data);  

  //check that the file could be opened okay.
  //status = 0 means it failed, 1 means it opened okay.
  if(!status){
    cout << "failed to get data from "<< data_file_name 
	 <<".. exiting"  << endl;
    return(1);
  }

  /***** set up a 2d array which will hold the magnitude information ****/
  //The image is a bit big which would slow down the fourier transforms.
  //So lets crop the edges; the image will go from 2048x2048 to 1024x1024.

  int nx = data.get_size_x()/2;
  int ny = data.get_size_y()/2;

  //declare the array and fill it.
  Double_2D intensity(nx,ny);

  //loop over the pixels in the image
  for(int i=0; i < nx; i++){
    for(int j=0; j< ny; j++){
      //copy to the new array
      //apply a threshold to the data to remove background
      intensity.set(i,j,data.get(i+nx/2,j+ny/2)-noise_level);
      if(intensity.get(i,j)<0)
	intensity.set(i,j,0);
    }
  }

  //write out the image before and after the crop and threashold to see 
  //what they look. "true" is used to indicate we want it ouput on log scale.
  //write_ppm("data_before.ppm", *data, true);
  //write_ppm("data_after.ppm", intensity, true);


  /****** get the support from a file and read it into an array *****/
  
  //Every pixel with a zero value is intepreted as being outside
  //the support

  Double_2D support;
  status = read_tiff(support_file_name, support);
  if(!status){
    cout << "failed to get data from "<< support_file_name 
	 <<".. exiting"  << endl;
    return(1);
  }

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

  //create the planar CDI object which will be used to
  //perform the reconstuction.
  PlanarCDI planar(object_estimate,4);
 
  //set the support and intensity
  planar.set_support(support);
  planar.set_intensity(intensity);

  //set the algorithm to hybrid input-output
  planar.set_algorithm(HIO);

  //Initialise the current object ESW with a random numbers
  //"0" is the seed to the random number generator
  planar.initialise_estimate(0);
  
  //make a 2D object. This will be used to output the 
  //image of the current estimate.
  Double_2D result(nx,ny);

  /******* for fun, let's get the autocorrelation *****/

  Double_2D autoc(nx,ny);
  planar.get_intensity_autocorrelation(autoc);
  write_ppm("test_autocorrelation.ppm", autoc, true); //"true" means log scale


  /*** run the reconstruction ************/

  for(int i=0; i<hio_iterations; i++){

    cout << "iteration " << i << endl;

    //apply the set of planar CDI projections 
    planar.iterate(); 
    cout << "Current error is "<<planar.get_error()<<endl;
    
    //every "output_iterations" 
    //output the current estimate of the object
    if(i%output_iterations==0){

      ostringstream temp_str ( ostringstream::out ) ;
      object_estimate.get_2d(MAG,result);
      temp_str << "real_example_iteration_" << i << ".ppm";
      write_ppm(temp_str.str(), result);

      temp_str.clear();

      //uncomment to output the estimated diffraction pattern
      /**Complex_2D * temp = object_estimate.clone();
      fft.perform_forward_fft(temp);
      temp->get_2d(MAG_SQ,&result);
      temp_str << "diffraction.ppm";
      write_ppm(temp_str.str(), nx, ny, result, true); 
      delete temp;
      **/
      
      //apply the shrinkwrap algorithm
      //1.5 is the gaussian width in pixels
      //0.1 is the threshold (10% of the maximum pixel).
      planar.apply_shrinkwrap(1.5,0.1);

    }

  }
  
  //now change to the error reduction algorithm 
  planar.set_algorithm(ER);

  for(int i=hio_iterations; i<(hio_iterations+er_iterations+1); i++){

    cout << "iteration " << i << endl;
    
    planar.iterate(); 
    
    cout << "Current error is "<<planar.get_error()<<endl;

    if(i%output_iterations==0){
      //output the current estimate of the object
      ostringstream temp_str ( ostringstream::out ) ;
      object_estimate.get_2d(MAG,result);
      temp_str << "real_example_iteration_" << i << ".ppm";
      write_ppm(temp_str.str(), result);
      temp_str.clear();

      //apply the shrinkwrap algorithm
      planar.apply_shrinkwrap(1.5,0.1);
    }
    
  }

  //And we are done. "object_estimate" contained the final estimate of
  //the ESW.

  Double_2D result2(nx,ny);

  double error=0;
  planar.get_best_result(0,error)->get_2d(MAG,result2);
  write_ppm("best_error.ppm", result2);
  cout << "Best error 0 is "<< error <<endl;

  planar.get_best_result(1,error)->get_2d(MAG,result2);
  write_ppm("best_error_1.ppm", result2);
  cout << "Best error 1 is "<< error <<endl;

  planar.get_best_result(2,error)->get_2d(MAG,result2);
  write_ppm("best_error_2.ppm", result2);
  cout << "Best error 2 is "<< error <<endl;

  planar.get_best_result(3,error)->get_2d(MAG,result2);
  write_ppm("best_error_3.ppm", result2);
  cout << "Best error 3 is "<< error <<endl;

  ProfilerStop();

  return 0;
}

