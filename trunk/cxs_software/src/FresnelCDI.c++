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
#include "TransmissionConstraint.h"
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
   illumination_at_sample(nx,ny),
   transmission(nx,ny),
   B_s(nx,ny),
   B_d(ny,ny){

  illumination.copy(white_field);
  set_experimental_parameters(wavelength,focal_detector_length,
			      focal_sample_length,pixel_length);

  set_algorithm(ER);
}

void FresnelCDI::set_experimental_parameters(double beam_wavelength,
					     double focal_detector_length,
					     double focal_sample_length,
					     double pixel_size){
  
  wavelength = beam_wavelength;
  pixel_length = pixel_size;

  double x_mid = (nx-1)/2;
  double y_mid = (ny-1)/2;

  double zfd = focal_detector_length;
  double zfs = focal_sample_length;
  double zsd = focal_detector_length - focal_sample_length;

  double zsd_ = 1/(1/zsd - 1/focal_detector_length);

  double scaling_x = beam_wavelength*zsd/(pixel_length*nx);
  double scaling_y = beam_wavelength*zsd/(pixel_length*ny);

  double rho_2_d;
  double phi_B_d;
  double phi_B_s;

  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){

      rho_2_d = pow(pixel_length*(x_mid-i),2) + 
	pow(pixel_length*(y_mid-j),2);

      phi_B_d = (M_PI*rho_2_d)/(beam_wavelength)*((1/focal_detector_length)-(1/zsd));
      phi_B_s = -phi_B_d;

      B_s.set_real(i,j,cos(phi_B_s));
      B_s.set_imag(i,j,sin(phi_B_s));

      B_d.set_real(i,j,cos(phi_B_d));
      B_d.set_imag(i,j,sin(phi_B_d));

    }
  }

  illumination_at_sample.copy(illumination);
  propagate_from_detector(illumination_at_sample);

}

void FresnelCDI::initialise_estimate(int seed){
  //initialise the random number generator
  srand(seed);
  
  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      
      //do a simple initialisation
      
      //make the magnitude the difference of the intensity and
      //white-field
      double amp = intensity_sqrt.get(i,j)-norm*illumination.get_mag(i,j);
      
      //perterb the phase a bit about zero to allow random starts.
      double phase = 0.2*2*M_PI*(0.5-rand()/((double) RAND_MAX) );
      complex.set_mag(i,j,amp);
      complex.set_phase(i,j,phase);
      
    }
  }

  //take the result to the sample plane and apply the support
  propagate_from_detector(complex);
  apply_support(complex);
  
}

void FresnelCDI::apply_support(Complex_2D & c){
  
  support_constraint(c);
  
  if(transmission_constraint){
    get_transmission_function(transmission,&c);
    //    (*custom_complex_constraint)(transmission);
    set_transmission_function(transmission,&c);
  }

}

void FresnelCDI::scale_intensity(Complex_2D & c){

  c.add(illumination,norm); //add the white field

  PlanarCDI::scale_intensity(c);
  
  c.add(illumination,-norm);//subtract the white field

}

void FresnelCDI::propagate_from_detector(Complex_2D & c){
  c.multiply(B_d);
  c.perform_backward_fft();
  c.invert(true);
}

void FresnelCDI::propagate_to_detector(Complex_2D & c){
  c.invert(true); 
  c.perform_forward_fft();
  c.multiply(B_s);
}

void FresnelCDI::set_transmission_function(Complex_2D & transmission,
					   Complex_2D * esw){

  if(!esw)
    esw=&complex;

  double ill_r;
  double trans_r;
  double ill_i;
  double trans_i;
  
  for(int i=0; i<nx; i++){
    for(int j=0; j<nx; j++){
      ill_r = norm*illumination_at_sample.get_real(i,j);
      trans_r = transmission.get_real(i,j);
      ill_i = norm*illumination_at_sample.get_imag(i,j);
      trans_i = transmission.get_imag(i,j);
   
      // ESW = TL - L
      // T - transmission function, L - illumination
      esw->set_real(i,j,ill_r*trans_r - ill_i*trans_i - ill_r);
      esw->set_imag(i,j,ill_r*trans_i + ill_i*trans_r - ill_i);      
    }
  }
}

void FresnelCDI::get_transmission_function(Complex_2D & result, 
					   Complex_2D * esw,
					   bool inforce_unity_mag){

  //divide the estimate by the illuminating wavefield and add unity.
 
  //the code below could written more eligantly with get_mag and get_phase etc.
  //but the code below is faster because it doesn't use the maths tan function.
  //This is important when complex constraints are applied

  if(!esw)
    esw=&complex;

  double ill_r;
  double esw_r;
  double ill_i;
  double esw_i;
  
  double real_numerator;
  double imag_numerator;
  double denom;

  for(int i=0; i<nx; i++){
    for(int j=0; j<nx; j++){
      if(norm!=0&&illumination_at_sample.get_mag(i,j)!=0){

	ill_r = norm*illumination_at_sample.get_real(i,j);
	esw_r = esw->get_real(i,j);
	ill_i = norm*illumination_at_sample.get_imag(i,j);
	esw_i = esw->get_imag(i,j);

	denom = ill_r*ill_r + ill_i*ill_i;
	real_numerator = ill_r*(esw_r+ill_r) + ill_i*(esw_i+ill_i);
	imag_numerator = ill_r*(esw_i+ill_i) - ill_i*(esw_r+ill_r);

	result.set_real(i,j,real_numerator/denom);
	result.set_imag(i,j,imag_numerator/denom);
	
	if(inforce_unity_mag && result.get_mag(i,j) > 1)
	  result.set_mag(i,j,1);      
      }
      else{
	result.set_real(i,j,0);
	result.set_imag(i,j,0);

      }
    }
  }

  if(transmission_constraint){
    transmission_constraint->apply_constraint(result);
  }
  
}
