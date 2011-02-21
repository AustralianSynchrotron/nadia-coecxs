#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <fftw3.h>
#include <cstdlib> 
#include "Complex_2D.h"
#include "Double_2D.h"
#include "FresnelCDI.h"
#include "io.h" //
#include <sstream>

using namespace std;

FresnelCDI::FresnelCDI(Complex_2D & initial_guess,
		       Complex_2D & white_field,
		       double beam_wavelength,
		       double focal_detector_length,
		       double focal_sample_length,
		       double pixel_size,
		       double normalisation,
		       int n_best)
  :PlanarCDI(initial_guess,n_best),
   wavelength(beam_wavelength),
   pixel_length(pixel_size),
   norm(normalisation),
   illumination(nx,ny),
   B_s(nx,ny),
   B_d(ny,ny){

  illumination.copy(white_field);

  double x_mid = (nx-1)/2;
  double y_mid = (ny-1)/2;

  double zfd = focal_detector_length;
  double zfs = focal_sample_length;
  double zsd = focal_detector_length - focal_sample_length;

  double zsd_ = 1/(1/zsd - 1/focal_detector_length);

  double scaling_x = beam_wavelength*zsd/(pixel_length*nx);
  double scaling_y = beam_wavelength*zsd/(pixel_length*ny);

  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){

      double rho_2_d = pow(pixel_length*(x_mid-i),2) + pow(pixel_length*(y_mid-j),2);

      double phi_B_d = -(M_PI*rho_2_d)/(beam_wavelength)*((1/focal_detector_length)-(1/zsd));
      double phi_B_s = -phi_B_d;

      B_s.set_real(i,j,cos(phi_B_s));
      B_s.set_imag(i,j,sin(phi_B_s));

      B_d.set_real(i,j,cos(phi_B_d));
      B_d.set_imag(i,j,sin(phi_B_d));

    }
  }

  set_algorithm(ER);

}

void FresnelCDI::initialise_estimate(int seed){
  //initialise the random number generator

  srand(seed);

  //start in the detector plane and use eq. 137 from
  //Harry's review paper.

  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){

      if(illumination.get_value(i,j,REAL)==0 || norm==0){
      	complex.set_real(i,j,0); 
	complex.set_imag(i,j,0); 
      }
      else{
	double real_value = intensity_sqrt.get(i,j)*intensity_sqrt.get(i,j) 
	  - norm*norm*illumination.get_value(i,j,MAG_SQ);

	real_value = real_value / (2*norm*illumination.get_value(i,j,REAL));

	complex.set_real(i,j,real_value); 
	complex.set_imag(i,j,0);
      }

    }
  }

  propagate_from_detector(complex);
  apply_support(complex);

}


void FresnelCDI::scale_intensity(Complex_2D & c){
  c.add(illumination,norm); //add the white field
  PlanarCDI::scale_intensity(c);
  c.add(illumination,-norm);//subtract the white field
}

void FresnelCDI::propagate_from_detector(Complex_2D & c){
  c.multiply(B_d);
  c.perform_forward_fft(); 
  c.invert(true); 
}

void FresnelCDI::propagate_to_detector(Complex_2D & c){
  c.invert(true);  
  c.perform_backward_fft();
  c.multiply(B_s);
}


void FresnelCDI::get_transmission_function(Complex_2D & result){

  //get the illuminating wavefield in the sample plane
  propagate_from_detector(illumination);

  //divide the estimate by the illuminating wavefield
  //and add unity.
  for(int i=0; i<nx; i++){
    for(int j=0; j<nx; j++){
      double new_mag = complex.get_mag(i,j)/illumination.get_mag(i,j);
      double new_phi = complex.get_value(i,j,PHASE) 
	- illumination.get_value(i,j,PHASE);
      result.set_real(i,j,new_mag*cos(new_phi)+1);
      result.set_imag(i,j,new_mag*sin(new_phi));
    }
  }

  //go back to the detector plane
  propagate_to_detector(illumination);

}
