/**
 * @file simulation_example.c
 * @author  Nadia Davidson <nadiamd@unimelb.edu.au>
 *
 * @section DESCRIPTION
 *
 * This file provides an example of running the planar CDI 
 * reconstruction on simulated data. Take a look at
 * real_example.c as well to see what the code does.
 *
 */

#include <iostream>
#include <math.h>
#include <string>
#include <cstdlib> 
#include "io.h"
#include "Complex_2D.h"
#include "Double_2D.h"
#include "PlanarCDI.h"
#include "FourierT.h"
#include <sstream>

using namespace std;


/**************************************/
int main(void){

  //define some constants which will be used in the code:

  //the data file name
  const static char * data_file_name = "image_files/object.tiff";

  //the approx. level of background noise in the image
  const double noise_level = 15;

  //the file which provides the support (pixels with the value 0
  //are considered as outside the object)
  const static char * support_file_name = "image_files/support_2.tiff";

  //number of hybrid input-out iterations to perform.
  const int hio_iterations = 1000;
   
  //output the current image ever "output_iterations"
  const int output_iterations = 40;


  /****** get the object from an image file ****************/

  //get the data from file
  
  Double_2D data;

  //read the data into an array
  int status = read_tiff(data_file_name, data);
  if(!status){
    cout << "failed.. exiting"  << endl;
    return(1);
  }

  int n_x = data.get_size_x();
  int n_y = data.get_size_y();
  
  //fill the complex no. with image data
  Complex_2D input(n_x,n_y);
  for(int i=0; i<n_x; i++){
    for(int j=0; j<n_y; j++){
      input.set_value(i,j,REAL, 1.0/sqrt(2.0)*data.get(i,j)*pow(-1,i + j));
      input.set_value(i,j,IMAG, 1.0/sqrt(2.0)*data.get(i,j)*pow(-1,i + j));
    }
  }

  /**** fourier transform the image to get the diffraction pattern *****/

  //fourier transform
  FourierT fft(n_x,n_y);
  fft.perform_forward_fft(input);

  //write the fourier transform to file.
  Double_2D intensity(n_x,n_y);
  input.get_2d(MAG_SQ,intensity);

  //apply a threashold to make the simulation a bit more realistic
  for(int i=0; i<n_x; i++){
    for(int j=0; j<n_y; j++){
      intensity.set(i,j,intensity.get(i,j)-noise_level);
      if(intensity.get(i,j)<0)
	intensity.set(i,j,0);
    }
  }

  //write the output to file (use log scale)
  write_ppm("sim_intensity.ppm",intensity,true);

  /******** get the support from file ****************************/

  Double_2D support(n_x,n_y);
  status = read_tiff(support_file_name, support);

  /*************** do the reconstruction *******************/

  //create a project object and set the options.
  Complex_2D first_guess(n_x,n_y);
  PlanarCDI my_planar(first_guess);
  my_planar.set_support(support);
  my_planar.set_intensity(intensity);
  my_planar.set_algorithm(HIO);

  //set the inital guess to be random inside the support
  //and zero outside. Note that this must be called
  //after "my_planar.set_support()"
  my_planar.initialise_estimate(0);
  
  //make a temporary arrary
  Double_2D result(n_x,n_y);

  //apply the projection operators
  for(int i=0; i<hio_iterations; i++){

    cout << "iteration " << i << endl;

    my_planar.iterate();

    if(i%output_iterations==0){

      ostringstream temp_str ( ostringstream::out ) ;
      first_guess.get_2d(MAG,result);
      temp_str << "sim_result_" << i << ".ppm";
      write_ppm(temp_str.str(), result);

    }
  }
  
  return 0;
}

