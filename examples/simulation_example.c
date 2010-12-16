/**
 * @file simulation_example.c
 * @author  Nadia Davidson <nadiamd@unimelb.edu.au>
 *
 * @section DESCRIPTION
 *
 * This file provides an example of running the planar CDI 
 * reconstruction on simulated data.
 *
 */

#include <iostream>
//#include <fstream>
#include <math.h>
#include <string>
//#include <stdlib.h>
//#include <fftw3.h>
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

  //write the output to file
  write_ppm("sim_intensity.ppm",intensity);


  /******** get the support from file ****************************/

  Double_2D support(n_x,n_y);
  status = read_tiff(support_file_name, support);

  /**** make the first guess: random with the support imposed ******/

  Complex_2D first_guess(n_x,n_y);
  for(int i=0; i<n_x; i++){
    for(int j=0; j<n_y; j++){
      if(!support.get(i,j)){
	first_guess.set_value(i,j,REAL,0); 
	first_guess.set_value(i,j,IMAG,0);
      }
      else{
	double r = (255.0*rand()/(double) RAND_MAX)* pow(-1,i + j);
	double im = (255.0*rand()/(double) RAND_MAX)* pow(-1,i + j);
	first_guess.set_value(i,j,REAL,r); 
	first_guess.set_value(i,j,IMAG,im);
      }
    }
  }

  /*************** do the reconstruction *******************/

  //create a project object and set the options.
  PlanarCDI proj(first_guess);
  proj.set_support(support);
  proj.set_intensity(intensity);
  proj.set_algorithm(HIO);
  
  //write the fourier transform to file.
  Double_2D result(n_x,n_y);

  //apply the iterations
  for(int i=0; i<hio_iterations; i++){

    cout << "iteration " << i << endl;

    proj.iterate();

    if(i%output_iterations==0){

      ostringstream temp_str ( ostringstream::out ) ;
      first_guess.get_2d(MAG,result);
      temp_str << "sim_result_" << i << ".ppm";
      write_ppm(temp_str.str(), result);
      //      temp_str.clear();

      /**Complex_2D * temp = first_guess.clone();
      fft.perform_forward_fft(temp);
      temp->get_2d(MAG,&result);
      temp_str << "diffraction.ppm";
      write_ppm(temp_str.str(), n_x, n_y, result, true);
      delete temp;**/
    }
  }

  

  for(int i=0; i<hio_iterations; i++){

    cout << "iteration " << i << endl;

    proj.iterate();

    if(i%output_iterations==0){

      ostringstream temp_str ( ostringstream::out ) ;
      first_guess.get_2d(MAG,result);
      temp_str << "sim_result_" << i << ".ppm";
      write_ppm(temp_str.str(), result);
      //      temp_str.clear();

      /**Complex_2D * temp = first_guess.clone();
      fft.perform_forward_fft(temp);
      temp->get_2d(MAG,&result);
      temp_str << "diffraction.ppm";
      write_ppm(temp_str.str(), n_x, n_y, result, true);
      delete temp;**/
    }
  }
  
  return 0;
}

