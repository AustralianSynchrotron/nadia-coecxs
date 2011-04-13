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

#define FORWARD  -1
#define BACKWARD +1

FresnelCDI::FresnelCDI(Complex_2D & initial_guess,
		       Complex_2D & white_field,
		       double beam_wavelength,
		       double focal_detector_length,
		       double focal_sample_length,
		       double pixel_size,
		       double normalisation,
		       int n_best)
  :PlanarCDI(initial_guess,n_best),
   illumination(nx,ny),
   norm(normalisation),
   coefficient(nx/2,ny/2)
   // B_d(ny,ny)
{


   //   illumination_at_sample(nx,ny),
   //    transmission(nx,ny),


  illumination.copy(white_field);
  set_experimental_parameters(beam_wavelength,focal_detector_length,
			      focal_sample_length,pixel_size);

  illumination_at_sample=0;
  transmission = 0;
  
  set_algorithm(ER);
}

FresnelCDI::~FresnelCDI(){

  if(illumination_at_sample)
    delete illumination_at_sample;
  
  if(transmission)
    delete transmission;

}


void FresnelCDI::set_experimental_parameters(double beam_wavelength,
					     double focal_detector_length,
					     double focal_sample_length,
					     double pixel_size){
  
  wavelength = beam_wavelength;
  pixel_length = pixel_size;
  this->focal_detector_length = focal_detector_length;
  this->focal_sample_length = focal_sample_length;
    
  double x_mid = (nx-1)/2;//.0;
  double y_mid = (ny-1)/2;//.0;

  double zfd = focal_detector_length;
  double zfs = focal_sample_length;
  double zsd = focal_detector_length - focal_sample_length;

  double factor = pixel_length*pixel_length*M_PI/(beam_wavelength)*((1/zfd) - (1/zsd));

  double phi;

  for(int i=0; i<nx/2.0; i++){
    for(int j=0; j<ny/2.0; j++){

      phi= factor*((x_mid-i)*(x_mid-i) + (y_mid-j)*(y_mid-j)); 

      coefficient.set_real(i,j,cos(phi));
      coefficient.set_imag(i,j,sin(phi));

    }
  }

}


void FresnelCDI::multiply_factors(Complex_2D & c, int direction){
  
  double x_mid = (nx-1)/2;//.0;
  double y_mid = (ny-1)/2;//.0;

  double old_real, old_imag;
  double coef_real, coef_imag;
  int i_, j_;

  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){

      old_real = c.get_real(i,j);
      old_imag = c.get_imag(i,j);

      int i_ = x_mid - fabs(x_mid - i);
      int j_ = y_mid - fabs(y_mid - j);

      coef_real = coefficient.get_real(i_,j_);
      coef_imag = direction*coefficient.get_imag(i_,j_);

      c.set_real(i,j,old_real*coef_real - old_imag*coef_imag);
      c.set_imag(i,j,old_imag*coef_real + old_real*coef_imag);

    }
  }

  
}


void FresnelCDI::auto_set_norm(){
  Double_2D temp(nx,ny);
  illumination.get_2d(MAG,temp);
  double wf_norm = temp.get_sum();
  double int_norm = intensity_sqrt.get_sum();
  
  norm = int_norm/wf_norm;

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
    if(!transmission)
      transmission = new Complex_2D(nx,ny);

    get_transmission_function(*transmission,&c);
    set_transmission_function(*transmission,&c);
  }

}

void FresnelCDI::scale_intensity(Complex_2D & c){

  c.add(illumination,norm); //add the white field

  PlanarCDI::scale_intensity(c);
  
  c.add(illumination,-norm);//subtract the white field

}

void FresnelCDI::propagate_from_detector(Complex_2D & c){
  //  c.multiply(B_d);
  multiply_factors(c,BACKWARD);
  c.perform_backward_fft();
  c.invert(true);
}

void FresnelCDI::propagate_to_detector(Complex_2D & c){
  c.invert(true); 
  c.perform_forward_fft();
  //c.multiply(B_s);
  multiply_factors(c,FORWARD);
}

void FresnelCDI::set_transmission_function(Complex_2D & transmission,
					   Complex_2D * esw){
  if(!esw)
    esw=&complex;

  if(!illumination_at_sample){
    illumination_at_sample = new Complex_2D(nx,ny);
    illumination_at_sample->copy(illumination);
    propagate_from_detector(*illumination_at_sample);
  }

  double ill_r;
  double trans_r;
  double ill_i;
  double trans_i;
  
  for(int i=0; i<nx; i++){
    for(int j=0; j<nx; j++){
      ill_r = norm*illumination_at_sample->get_real(i,j);
      trans_r = transmission.get_real(i,j);
      ill_i = norm*illumination_at_sample->get_imag(i,j);
      trans_i = transmission.get_imag(i,j);
   
      // ESW = TL - L
      // T - transmission function, L - illumination
      esw->set_real(i,j,ill_r*trans_r - ill_i*trans_i - ill_r);
      esw->set_imag(i,j,ill_r*trans_i + ill_i*trans_r - ill_i);      
    }
  }
}

const Complex_2D & FresnelCDI::get_illumination_at_sample(){ 
  if(!illumination_at_sample){
    illumination_at_sample = new Complex_2D(nx,ny);
    illumination_at_sample->copy(illumination);
    propagate_from_detector(*illumination_at_sample);
  }
  
  return *illumination_at_sample;
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

  if(!illumination_at_sample){
    illumination_at_sample = new Complex_2D(nx,ny);
    illumination_at_sample->copy(illumination);
    propagate_from_detector(*illumination_at_sample);
  }

  double ill_r;
  double esw_r;
  double ill_i;
  double esw_i;
  
  double real_numerator;
  double imag_numerator;
  double denom;

  for(int i=0; i<nx; i++){
    for(int j=0; j<nx; j++){
      if(norm!=0&&illumination_at_sample->get_mag(i,j)!=0){

	ill_r = norm*illumination_at_sample->get_real(i,j);
	esw_r = esw->get_real(i,j);
	ill_i = norm*illumination_at_sample->get_imag(i,j);
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
