// Copyright 2011 Nadia Davidson 
// for The ARC Centre of Excellence in Coherent X-ray Science. 
//
// This program is distributed under the GNU General Public License. 
// We also ask that you cite this software in publications where you made 
// use of it for any part of the data analysis.
//
// date last modified: 17/01/2011

/**
 * @file PlanarCDI_simulation_example.c
 *
 * \a PlanarCDI_simulation_example.c The object from the
 * PlanarCDI_example (in fact I used the reconstructed image as the
 * object) is used to simulate a diffraction pattern for planar
 * CDI. The diffraction pattern is thresholded to make it more
 * realistic and then CDI reconstruction is performed.
 *
 */

#include <iostream>
#include <math.h>
#include <string>
#include <cstdlib> 
#include "io.h"
#include "Complex_2D.h"
#include "Double_2D.h"
#include "PolyCDI.h"
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
  const static char * support_file_name = "image_files/planar_support.tiff";

  const static char * data_spectrum_name = "image_files/spectrum.txt";
  //image_files/planar_support.tiff";

  //number of hybrid input-out iterations to perform.
  const int hio_iterations = 200;

  //the number of error-reduction iterations to perform.
  const int er_iterations = 100;

  //output the current image ever "output_iterations"
  const int output_iterations = 50;


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
      input.set_value(i,j,REAL, 1.0/sqrt(2.0)*data.get(i,j));
      input.set_value(i,j,IMAG, 1.0/sqrt(2.0)*data.get(i,j));
    }
  }

  /**** create the projection/reconstruction object *****/

  //Double_2D spectrum;
  //read_spec(data_spectrum_name, spectrum);

  double n=500.0;

  Double_2D spectrum(n, 2);

  double bw, sigma, mean, scale, del;

  mean=1.4;
  bw=0.1*mean;
  sigma=bw/2.35482; //FWHM in to sigma
  scale=1.0/(sqrt(2.0*3.14159));
  //    del=(1.6-1.2)/n;
  del=5.0*sigma/n;

  for(double i=0; i<n; i++){


    spectrum.set(int(i), 0, 1.0/(mean-(n/2.0-i)*del));
    spectrum.set(int(i), 1, scale*exp(-(del*(i-(n/2))*del*(i-(n/2))/2.0/sigma/sigma)));

    std::cout<<1/spectrum.get(i, 0)<<" "<<spectrum.get(i, 1)<<"\n";

  }



  Complex_2D first_guess(n_x,n_y);
  PolyCDI my_poly(first_guess, spectrum, 0.9, 4, 0);

  //propagate to the detector plane
  my_poly.propagate_to_detector(input);
  my_poly.expand_wl(input);


  //write the fourier transform to file.
  Double_2D intensity(n_x,n_y);
  intensity=my_poly.get_intensity(); 

  //apply a threshold to make the simulation a bit more realistic
  for(int i=0; i<n_x; i++){
    for(int j=0; j<n_y; j++){
      intensity.set(i,j,intensity.get(i,j)*intensity.get(i, j)-noise_level);
      if(intensity.get(i,j)<0)
	intensity.set(i,j,0);
    }
  }

  //write the output to file (use log scale)
  write_tiff("poly_sim_intensity.tiff",intensity,false);
  write_ppm("poly_sim_intensity.ppm", intensity, true);

  /******** get the support from file ****************************/
  /*
     Double_2D support(n_x,n_y);
     status = read_tiff(support_file_name, support);

  /*************** do the reconstruction *******************/
  /*
  //create a project object and set the options.
  my_poly.set_support(support);
  my_poly.set_intensity(intensity);
  my_poly.set_algorithm(HIO);

  //set the inital guess to be random inside the support
  //and zero outside. Note that this must be called
  //after "my_poly.set_support()"
  my_poly.initialise_estimate(0);

  //make a temporary arrary
  Double_2D result(n_x,n_y);

  //apply the projection operators
  for(int i=0; i<hio_iterations; i++){

  cout << "iteration " << i << endl;

  my_poly.iterate();

  if(i%output_iterations==0){

  ostringstream temp_str ( ostringstream::out ) ;
  first_guess.get_2d(MAG,result);
  temp_str << "sim_result_" << i << ".ppm";
  write_ppm(temp_str.str(), result);

  my_poly.apply_shrinkwrap();

  }
  }

  my_poly.set_algorithm(ER);

  for(int i=hio_iterations; i<(er_iterations+hio_iterations+1); i++){

  cout << "iteration " << i << endl;

  my_poly.iterate();

  if(i%output_iterations==0){

  ostringstream temp_str ( ostringstream::out ) ;
  first_guess.get_2d(MAG,result);
  temp_str << "sim_result_" << i << ".ppm";
  write_ppm(temp_str.str(), result);

  my_poly.apply_shrinkwrap();

  }
  }


   */

  return 0;
}

