#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <fftw3.h>
#include <cstdlib> 
#include "Complex_2D.h"
#include "FourierT.h"
#include "FresnelCDI_WF.h"
#include "io.h" //
#include <sstream>

using namespace std;

FresnelCDI_WF::FresnelCDI_WF(Complex_2D & initial_guess,
				 double beam_wavelength,
				 double zone_focal_length,
				 double focal_detector_length,
				 double pixel_size)
                                 :PlanarCDI(initial_guess){


  wavelength = beam_wavelength;
  zone_to_focal_length = zone_focal_length;
  focal_to_detector_length = focal_detector_length;
  pixel_length = pixel_size;
  
  //set-up the coefficients
  //it's easier to do it once and reuse the matrix.

  int nx = initial_guess.get_size_x();
  int ny = initial_guess.get_size_y();

  forward_coefficients_const = new Complex_2D(nx,ny);
  backward_coefficients_const = new Complex_2D(nx,ny);

  double x_mid = (nx-1)/2;
  double y_mid = (ny-1)/2;

  double z12 = zone_to_focal_length;
  double z23 = focal_to_detector_length;

  double scaling_x = beam_wavelength*z23/(pixel_length*nx);
  double scaling_y = beam_wavelength*z23/(pixel_length*ny);


  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){

      double rho_2 = pow(scaling_x*(x_mid-i),2) + pow(scaling_y*(y_mid-j),2);

      double phi = 2*z12 + 2*z23 + rho_2/z12 + rho_2/z23 ;
      phi*= M_PI/beam_wavelength;

      forward_coefficients_const->set_real(i,j,sin(phi));
      forward_coefficients_const->set_imag(i,j,-cos(phi));

      backward_coefficients_const->set_real(i,j,sin(phi));
      backward_coefficients_const->set_imag(i,j,cos(phi));

    }
  }

}


FresnelCDI_WF::~FresnelCDI_WF(){
  delete forward_coefficients_const;
  delete backward_coefficients_const;

  cout << "In the FresnelCDI_WF destructor" << endl;
}

void FresnelCDI_WF::initialise_estimate(int seed){
  //initialise the random number generator
  srand(seed);

  int nx = complex.get_size_x();
  int ny = complex.get_size_y();

  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      if(!support[i][j]){ //enforce the support condition on the inital guess
	complex.set_value(i,j,REAL,0); 
	complex.set_value(i,j,IMAG,0);
      }
      else{
	double r = intensity_sqrt[i][j]/sqrt(2.0);//(65000.0*rand()/(double) RAND_MAX) ;//* pow(-1,i + j);
	double im = intensity_sqrt[i][j]/sqrt(2.0);//(65000.0*rand()/(double) RAND_MAX) ;//* pow(-1,i + j);
	complex.set_value(i,j,REAL,r); 
	complex.set_value(i,j,IMAG,im);
      }
    }
  }
}

int FresnelCDI_WF::iterate(){
  
  int nx = complex.get_size_x();
  int ny = complex.get_size_y();

  double ** result = new double*[nx];
  for(int i=0; i < nx; i++)
    result[i]= new double[ny];

  scale_intensity(complex);
  cout << "Current Error is: " << current_error<<endl;

  //start with wavefield in the detector plane
  //go to the focal plane
  complex.invert();
  fft.perform_forward_fft(&complex);
  complex.invert();

  complex.multiply(*backward_coefficients_const);

  complex.invert();
  //go back to zone plate plane. 
  fft.perform_forward_fft(&complex);
  complex.invert();

  //complex.get_2d(MAG,&result);
  //write_ppm("zone_reco.ppm", nx, ny, result);
 
  //apply support constraint
  apply_support(complex);
  
  //go to the focal plane again.
  complex.invert();
  fft.perform_forward_fft(&complex);
  complex.invert();
  complex.multiply(*forward_coefficients_const);

  //complex.get_2d(MAG,&result);
  //write_ppm("focal_reco.ppm", nx, ny, result,true);

  complex.invert();

  //and back to the detector plane
  fft.perform_forward_fft(&complex);

  complex.invert();

  for(int i=0; i < nx; i++)
    delete [] result[i];
  delete [] result;

  return SUCCESS;
}
