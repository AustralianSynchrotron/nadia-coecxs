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
			     double pixel_size,
			     int n_best)
  :PlanarCDI(initial_guess,n_best),
   wavelength(beam_wavelength),
   zone_to_focal_length(zone_focal_length),
   focal_to_detector_length(focal_detector_length),
   pixel_length(pixel_size),
   forward_coefficient(nx,ny),
   backward_coefficient(nx,ny)
{


  //set-up the coefficients
  //it's easier to do it once and reuse the matrix.

  double x_mid = (nx-1)/2;
  double y_mid = (ny-1)/2;

  double z12 = zone_to_focal_length;
  double z23 = focal_to_detector_length;

  double scaling_x = beam_wavelength*z23/(pixel_length*nx);
  double scaling_y = beam_wavelength*z23/(pixel_length*ny);

  double norm = 1/(sqrt(nx*ny));

  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){

      double rho_2 = pow(scaling_x*(x_mid-i),2) + pow(scaling_y*(y_mid-j),2);
      //double phi = 2*z12 + 2*z23 + rho_2/z12 + rho_2/z23 ;

      double phi = rho_2/z12 + rho_2/z23 ;

      phi*= M_PI/beam_wavelength;

      forward_coefficient.set_real(i,j,sin(phi)*norm);
      forward_coefficient.set_imag(i,j,-cos(phi)*norm);

      backward_coefficient.set_real(i,j,-sin(-phi)*norm);
      backward_coefficient.set_imag(i,j,cos(-phi)*norm);

    }
  }
  
  /**  Double_2D result(nx,ny);
  forward_coefficients_const->get_2d(PHASE,result);
  write_ppm("forward_coefficients_phase.ppm",result);**/

}


void FresnelCDI_WF::initialise_estimate(int seed){
  //initialise the random number generator
  srand(seed);

  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      if(!support.get(i,j)){ //enforce the support condition on the inital guess
	complex.set_value(i,j,REAL,0); 
	complex.set_value(i,j,IMAG,0);
      }
      else{
	//double r = intensity_sqrt.get(i,j)///sqrt(2.0);//(65000.0*rand()/(double) RAND_MAX) ;//* pow(-1,i + j);
	//double im = 0;
	//double im = intensity_sqrt.get(i,j)/sqrt(2.0);//(65000.0*rand()/(double) RAND_MAX) ;//* pow(-1,i + j);
	double random = (65000.0*rand()/(double) RAND_MAX);

	//double r = (65000.0*rand()/(double) RAND_MAX) ;//* pow(-1,i + j);
	//double im = (65000.0*rand()/(double) RAND_MAX) ;//* pow(-1,i + j);
	//	complex.set_value(i,j,REAL,r); 
	//	complex.set_value(i,j,IMAG,im);
	complex.set_value(i,j,REAL,random); 
	complex.set_value(i,j,IMAG,0);
      }
    }
  }
}

int FresnelCDI_WF::iterate(){
  
  propagate_from_detector(complex);

  apply_support(complex);
  
  propagate_to_detector(complex);

  scale_intensity(complex);
  
  return SUCCESS;
}

void FresnelCDI_WF::propagate_from_detector(Complex_2D & c){
  //go to the focal plane
  c.perform_backward_fft();
  c.invert(true);  
  c.multiply(backward_coefficient);

  //go back to zone plate plane. 
  c.perform_backward_fft();

}

void FresnelCDI_WF::propagate_to_detector(Complex_2D & c){

  //go to the focal plane again.
  c.perform_forward_fft();
  c.multiply(forward_coefficient);

  //and back to the detector plane
  c.invert(true);
  c.perform_forward_fft();

}
