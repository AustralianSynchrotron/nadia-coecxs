// Copyright 2012 T'Mir Julius 
// for The ARC Centre of Excellence in Coherent X-ray Science. 
//
// This program is distributed under the GNU General Public License. 
// We also ask that you cite this software in publications where you made 
// use of it for any part of the data analysis.
//
// date last modified: 07/12/2012

/**
 * @file PolyCDI_simulation_example.c
 *
 * \a PolycCDI_simulation_example.c The object from the
 * PlanarCDI_example (in fact I used the reconstructed image as the
 * object) is used to simulate a diffraction pattern for polychromatic
 * CDI. A beam is simulated according to the properties of the beam of
 * thee Australain Synchrotron. The spectrum for this was generated using 
 * SPECTRA (http://radiant.harima.riken.go.jp/spectra/) and saved as a 
 * text file. The diffraction pattern is thresholded to make it more 
 * realistic.
 *
 */

#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdlib> 
#include <io.h>
#include <Complex_2D.h>
#include <Double_2D.h>
#include <PolyCDI.h>
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

  //read in the spectrum from a Spectra text file
  const static char * data_spectrum_name = "image_files/spectrum.txt";


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

  Double_2D spectrum;
  read_spec(data_spectrum_name, spectrum);

  Complex_2D first_guess(n_x,n_y);
  PolyCDI my_poly(first_guess, 0.9, 4, 0);

  my_poly.set_spectrum(data_spectrum_name);
  first_guess=input;

  //propagate to the detector plane
  my_poly.propagate_to_detector(input);
  my_poly.expand_wl(input);


  //write the fourier transform to file.
  Double_2D intensity(n_x,n_y);
  intensity=my_poly.get_intensity(); 

  //apply a threshold to make the simulation a bit more realistic
/*  for(int i=0; i<n_x; i++){
    for(int j=0; j<n_y; j++){
      intensity.set(i,j,intensity.get(i,j)*intensity.get(i, j)-noise_level);
      if(intensity.get(i,j)<0)
	intensity.set(i,j,0);
    }
  }
*/
  //write the output to file (use log scale)
  write_tiff("poly_sim_intensity.tiff",intensity,false);
  write_ppm("poly_sim_intensity.ppm", intensity, true);

  return 0;
}

