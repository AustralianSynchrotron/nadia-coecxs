#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <fftw3.h>
#include <cstdlib> 
#include "Complex_2D.h"
#include "Double_2D.h"
#include "FourierT.h"
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
		       double normalisation)
  :PlanarCDI(initial_guess),
   illumination(white_field),
   wavelength(beam_wavelength),
   pixel_length(pixel_size),
   norm(normalisation),
   B_s(nx,ny),
   B_d(ny,ny)
{

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

}

FresnelCDI::~FresnelCDI(){
  //  delete B_s;
  //  delete B_d;
}

void FresnelCDI::initialise_estimate(int seed){
  //initialise the random number generator
  srand(seed);

  //start in the detector plane and use eq. 137 from
  //Harry's review paper.

  double sum_int = 0;

  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
 
      double real_value = (intensity_sqrt[i][j]*intensity_sqrt[i][j] - norm*norm*illumination.get_value(i,j,MAG_SQ));
      sum_int +=intensity_sqrt[i][j];

      real_value = real_value / (2*norm*illumination.get_value(i,j,REAL));

      if(illumination.get_value(i,j,REAL)==0){
      	complex.set_real(i,j,0); 
	complex.set_imag(i,j,0); 
      }
      else{
	complex.set_real(i,j,real_value); 
	complex.set_imag(i,j,0);
      }
    }
  }

  Double_2D result(nx,ny);
  complex.get_2d(MAG,result);
  write_ppm("first_guess_bf.ppm",result);

  complex.multiply(B_d);
  fft.perform_forward_fft(&complex);
  complex.invert();  

  apply_support(complex);

  complex.get_2d(MAG,result);
  write_ppm("first_guess_after.ppm",result); 
}



//int FresnelCDI::iterate(){

//  PhaseRetrievalBase::iterate();
  
  //project_intensity(complex);

//  return SUCCESS;
// }

void FresnelCDI::project_intensity(Complex_2D & c){

  Double_2D result(complex.get_size_x(),complex.get_size_y());

  //Forward propgate
  c.invert();
  fft.perform_backward_fft(&c);
  c.multiply(B_s);

  c.get_2d(MAG,result);
  write_ppm("1-forward.ppm",result);
  c.get_2d(PHASE,result);
  write_ppm("1-forward_p.ppm",result);

  c.add(illumination,norm);

  c.get_2d(MAG,result);
  write_ppm("2-with_illum.ppm",result);
  c.get_2d(PHASE,result);
  write_ppm("2-with_illum_p.ppm",result);

  scale_intensity(c);

  c.get_2d(MAG,result);
  write_ppm("3-scaled.ppm",result);
  c.get_2d(PHASE,result);
  write_ppm("3-scaled_p.ppm",result);

  c.add(illumination,-norm);

  c.get_2d(MAG,result);
  write_ppm("4-subtracted.ppm",result);
  c.get_2d(PHASE,result);
  write_ppm("4-subtracted_p.ppm",result);

  //backward propogate
  c.multiply(B_d);
  fft.perform_forward_fft(&c);
  c.invert();  

  c.get_2d(MAG,result);
  write_ppm("5-backward.ppm",result);
  c.get_2d(PHASE,result);
  write_ppm("5-backward_p.ppm",result);

}

