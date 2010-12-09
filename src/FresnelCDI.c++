#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <fftw3.h>
#include <cstdlib> 
#include "Complex_2D.h"
#include "FourierT.h"
#include "FresnelCDI.h"
#include "io.h" //
#include <sstream>

using namespace std;

FresnelCDI::FresnelCDI(Complex_2D * initial_guess,
		       Complex_2D * white_field,
		       double beam_wavelength,
		       double focal_detector_length,
		       double focal_sample_length,
		       double pixel_size,
		       double normalisation)
                       :PhaseRetrievalBase(initial_guess){


  norm = normalisation;

  wavelength = beam_wavelength;
  pixel_length = pixel_size;
  
  illumination = white_field;

  //cout << sample_to_detector_length <<endl;
  //set-up the coefficients
  //it's easier to do it once and reuse the matrix.

  int nx = initial_guess->get_size_x();
  int ny = initial_guess->get_size_y();

  A_s = new Complex_2D(nx,ny);
  B_s = new Complex_2D(nx,ny);
  A_d = new Complex_2D(nx,ny);
  B_d = new Complex_2D(nx,ny);

  wf_curvature = new Complex_2D(nx,ny);

  double x_mid = (nx-1)/2;
  double y_mid = (ny-1)/2;

  double zfd = focal_detector_length;
  double zfs = focal_sample_length;
  double zsd = focal_detector_length - focal_sample_length;

  double zsd_ = 1/(1/zsd - 1/focal_detector_length);
  //cout << "zsd " <<zsd<<endl;

  double scaling_x = beam_wavelength*zsd/(pixel_length*nx);
  double scaling_y = beam_wavelength*zsd/(pixel_length*ny);

  cout << "scaling = "<<scaling_x<<endl;

  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){

      double rho_2_s = pow(scaling_x*(x_mid-i),2) + pow(scaling_y*(y_mid-j),2);
      double rho_2_d = pow(pixel_length*(x_mid-i),2) + pow(pixel_length*(y_mid-j),2);

      double phi_B_s = (M_PI*rho_2_s)/(beam_wavelength*zsd);
      double phi_B_d = (-M_PI*rho_2_d)/(beam_wavelength*zsd);

      double phi_A_s = (M_PI/beam_wavelength)*(-2*zsd-rho_2_s/zsd);
      double phi_A_d = (M_PI/beam_wavelength)*(2*zsd+rho_2_d/zsd);

      double phi_wf_curv = (M_PI*rho_2_d/beam_wavelength)*(1/zsd_);

      B_s->set_real(i,j,cos(phi_B_s));
      B_s->set_imag(i,j,sin(phi_B_s));

      B_d->set_real(i,j,cos(phi_B_d));
      B_d->set_imag(i,j,sin(phi_B_d));

      A_s->set_real(i,j,-sin(phi_A_s));
      A_s->set_imag(i,j,cos(phi_A_s));

      A_d->set_real(i,j,sin(phi_A_d));
      A_d->set_imag(i,j,-cos(phi_A_d));

      wf_curvature->set_real(i,j,cos(phi_wf_curv));
      wf_curvature->set_imag(i,j,sin(phi_wf_curv));

    }
  }
  //illumination->multiply(wf_curvature);
}

FresnelCDI::~FresnelCDI(){
  delete B_s;
  delete B_d;
  delete A_s;
  delete A_d;
  delete wf_curvature;

  cout << "In the FresnelCDI destructor" << endl;
}

void FresnelCDI::initialise_estimate(int seed){
  //initialise the random number generator
  srand(seed);

  int nx = complex->get_size_x();
  int ny = complex->get_size_y();

  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      if(!support[i][j]){ //enforce the support condition on the inital guess
	complex->set_value(i,j,REAL,0); 
	complex->set_value(i,j,IMAG,0);
      }
      else{
	//double r = intensity_sqrt[i][j]/sqrt(2.0);//(65000.0*rand()/(double) RAND_MAX) ;//* pow(-1,i + j);
	//double im = intensity_sqrt[i][j]/sqrt(2.0);//(65000.0*rand()/(double) RAND_MAX) ;//* pow(-1,i + j);
	
	double r = (65000.0*rand()/(double) RAND_MAX);//* pow(-1,i + j) ;
	double im = (65000.0*rand()/(double) RAND_MAX);//* pow(-1,i + j) ;

	complex->set_value(i,j,REAL,r); 
	complex->set_value(i,j,IMAG,im);
      }
    }
  }
}


void FresnelCDI::project_intensity(Complex_2D * c){
  int nx = complex->get_size_x();
  int ny = complex->get_size_y();

  double ** result = new double *[nx];
  for(int i=0; i<nx; i++)
    result[i]=new double[ny];

  /**illumination->get_2d(MAG,&result);
  write_ppm("wf_mag_d.ppm",nx,ny,result); 
  illumination->get_2d(PHASE,&result);
  write_ppm("wf_phase_d.ppm",nx,ny,result); 
  
  backward_propogate(illumination);

  illumination->get_2d(MAG,&result);
  write_ppm("wf_mag_s.ppm",nx,ny,result,true); 
  illumination->get_2d(PHASE,&result);
  write_ppm("wf_phase_s.ppm",nx,ny,result); 
  
  forward_propogate(illumination);**/

  /**B_d->get_2d(PHASE,&result);
  write_ppm("B_d.ppm",nx,ny,result);
  B_s->get_2d(PHASE,&result);
  write_ppm("B_s.ppm",nx,ny,result);
  A_d->get_2d(PHASE,&result);
  write_ppm("A_d.ppm",nx,ny,result);
  A_s->get_2d(PHASE,&result);
  write_ppm("A_s.ppm",nx,ny,result);**/

  c->invert();

  c->multiply(wf_curvature);

  forward_propogate(c);
  
  c->get_2d(MAG,&result);
  write_ppm("1-forward.ppm",nx,ny,result);

  c->add(illumination,norm);

  c->get_2d(MAG,&result);
  write_ppm("2-with_illum.ppm",nx,ny,result);

  scale_intensity(c);

  c->get_2d(MAG,&result);
  write_ppm("3-scaled.ppm",nx,ny,result);

  c->add(illumination,-norm);

  c->get_2d(MAG,&result);
  write_ppm("4-subtracted.ppm",nx,ny,result);

  backward_propogate(c);

  c->get_2d(MAG,&result);
  write_ppm("5-backward.ppm",nx,ny,result);


}

//FIXME! sample to detector
void FresnelCDI::forward_propogate(Complex_2D * c){
  //c->multiply(B_s);
  //c->invert();
  fft->perform_forward_fft(c);
  //c->invert();
  //c->multiply(A_d);
  cout << "we are in this function"<<endl;
}

//FIXME! detector to sample
void FresnelCDI::backward_propogate(Complex_2D * c){

  //c->multiply(B_d);

  int nx = complex->get_size_x();
  int ny = complex->get_size_y();

  double ** result = new double *[nx];
  for(int i=0; i<nx; i++)
    result[i]=new double[ny];

  /**illumination->get_2d(MAG,&result);
  write_ppm("wf_mag_d_1.ppm",nx,ny,result); 
  illumination->get_2d(PHASE,&result);
  write_ppm("wf_phase_d_1.ppm",nx,ny,result); **/

  //c->invert();

  /**illumination->get_2d(MAG,&result);
  write_ppm("wf_mag_d_2.ppm",nx,ny,result); 
  illumination->get_2d(PHASE,&result);
  write_ppm("wf_phase_d_2.ppm",nx,ny,result); **/

  fft->perform_backward_fft(c);

  /**illumination->get_2d(MAG,&result);
  write_ppm("wf_mag_d_3.ppm",nx,ny,result); 
  illumination->get_2d(PHASE,&result);
  write_ppm("wf_phase_d_3.ppm",nx,ny,result); **/

  c->invert();

  /**illumination->get_2d(MAG,&result);
  write_ppm("wf_mag_d_4.ppm",nx,ny,result); 
  illumination->get_2d(PHASE,&result);
  write_ppm("wf_phase_d_4.ppm",nx,ny,result); **/

  //c->multiply(A_s);

  /**illumination->get_2d(MAG,&result);
  write_ppm("wf_mag_d_5.ppm",nx,ny,result); 
  illumination->get_2d(PHASE,&result);
  write_ppm("wf_phase_d_5.ppm",nx,ny,result); **/
}

