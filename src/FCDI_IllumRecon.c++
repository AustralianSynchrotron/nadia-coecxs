#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <fftw3.h>
#include <cstdlib> 
#include "Complex_2D.h"
#include "FourierT.h"
#include "FCDI_IllumRecon.h"
#include "io.h" //
#include <sstream>

using namespace std;

FCDI_IllumRecon::FCDI_IllumRecon(Complex_2D * initial_guess,
				 double beam_wavelength,
				 double zone_focal_length,
				 double focal_detector_length,
				 double pixel_size){
  complex = initial_guess;
  wavelength = beam_wavelength;
  zone_to_focal_length = zone_focal_length;
  focal_to_detector_length = focal_detector_length;
  pixel_length = pixel_size;
  
  int nx = complex->get_size_x();
  int ny = complex->get_size_y();

  support = new double*[nx];
  intensity_sqrt = new double*[nx];
  for(int i=0; i< nx; i++){
    support[i]=new double[ny];
    intensity_sqrt[i]=new double[ny];
  }

  //make the fourier transform object
  fft = new FourierT(complex->get_size_x(), complex->get_size_y());


  //forward_coefficients = new Complex_2D(nx,ny);
  //backward_coefficients = new Complex_2D(nx,ny);
  forward_coefficients_const = new Complex_2D(nx,ny);
  backward_coefficients_const = new Complex_2D(nx,ny);
  A_forward = new Complex_2D(nx,ny);
  A_backward = new Complex_2D(nx,ny);
  B_forward = new Complex_2D(nx,ny);
  B_backward = new Complex_2D(nx,ny);

  //set-up the coefficients
  //it's easier to do it once and reuse the matrix.
  double x_mid = (nx-1)/2;
  double y_mid = (ny-1)/2;

  double z12 = zone_to_focal_length;
  double z21 = -1*zone_to_focal_length;
  double z23 = focal_to_detector_length;
  double z32 = -1*focal_to_detector_length;

  double scaling_x = beam_wavelength*z23/(pixel_length*nx);
  double scaling_y = beam_wavelength*z23/(pixel_length*ny);

  //double scaling = pixel_length;

  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){

      double rho_2 = pow(scaling_x*(x_mid-i),2) + pow(scaling_y*(y_mid-j),2);

      double phi_A_f = (M_PI/beam_wavelength)*(2*z12 + rho_2/z12);
      double phi_B_f = (M_PI/beam_wavelength)*(rho_2/z23);

      double phi_A_b = (M_PI/beam_wavelength)*(2*z32 + rho_2/z32);
      double phi_B_b = (M_PI/beam_wavelength)*(rho_2/z21);


      //      cout << "phi_B_f="<< phi_B_f;
      // cout << "sanity check: " << cos(phi_B_f) << endl;
      //cout << "sanity check: " << cos(phi_B_f) << endl;

      
      A_forward->set_real(i,j,sin(phi_A_f));
      A_forward->set_imag(i,j,-1*cos(phi_A_f));

      A_backward->set_real(i,j,sin(phi_A_b));
      A_backward->set_imag(i,j,-1*cos(phi_A_b));

      B_forward->set_real(i,j,rho_2);//cos(phi_B_f));
      B_forward->set_imag(i,j,rho_2);//sin(phi_B_f));

      B_backward->set_real(i,j,cos(phi_B_b));
      B_backward->set_imag(i,j,sin(phi_B_b));
      
      forward_coefficients_const->set_real(i,j,sin(2*M_PI*z23/beam_wavelength));
      forward_coefficients_const->set_imag(i,j,-1*cos(2*M_PI*z23/beam_wavelength));

      backward_coefficients_const->set_real(i,j,sin(2*M_PI*z21/beam_wavelength));
      backward_coefficients_const->set_imag(i,j,-1*cos(2*M_PI*z21/beam_wavelength));

      /**  double phi_a_12 = (M_PI/beam_wavelength)*(rho_2/fabs(z12));
      double phi_a_23 = (M_PI/beam_wavelength)*(rho_2/fabs(z23));

      A_forward->set_real(i,j,cos(phi_a_12));
      A_forward->set_imag(i,j,sin(phi_a_12));
      A_backward->set_real(i,j,cos(phi_a_12));
      A_backward->set_imag(i,j,-sin(phi_a_12));

      B_forward->set_real(i,j,cos(phi_a_23));
      B_forward->set_imag(i,j,sin(phi_a_23));
      B_backward->set_real(i,j,cos(phi_a_23));
      B_backward->set_imag(i,j,-sin(phi_a_23));**/
    }
  }
  
}


FCDI_IllumRecon::~FCDI_IllumRecon(){
  delete fft;
  //delete forward_coefficients;
  //delete backward_coefficients;
  delete forward_coefficients_const;
  delete backward_coefficients_const;
  delete A_forward;
  delete B_forward;
  delete A_backward;
  delete B_backward;

  for(int i=0; i<complex->get_size_x() ; i++){
    delete[] support[i];
    delete[] intensity_sqrt[i];
  }
  delete[] support;
  delete[] intensity_sqrt;


  //  for(int i=0; i<NTERMS-1; i++)
  //  delete x[i];

}

double ** FCDI_IllumRecon::get_intensity_autocorrelation(){

  //make a temporary object
  int nx = complex->get_size_x();
  int ny = complex->get_size_y();

  Complex_2D temp_intensity(nx,ny);
  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      //set the real and imaginary components using the magnitude.
      double component = (1.0/sqrt(2.0))*(intensity_sqrt[i][j]*intensity_sqrt[i][j]);
      //scale by "pow(-1,i + j)" to make the fourier transform centered.
      component*=pow(-1,i + j);
      temp_intensity.set_value(i,j,REAL, component);
      temp_intensity.set_value(i,j,IMAG, component);
    }
  }
  
  // fourier transform the intensity 
  fft->perform_backward_fft(&temp_intensity);  

  //allocate some memory for the output.
  double ** autoc = new double*[nx];
  for(int i=0; i < nx; i++)
    autoc[i]= new double[ny];

  //get the magnitude of the fourier transformed data.
  temp_intensity.get_2d(MAG, &autoc);

  return autoc;
}


void FCDI_IllumRecon::initialise_estimate(int seed){
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
	double r = (65000.0*rand()/(double) RAND_MAX) ;//* pow(-1,i + j);
	double im = (65000.0*rand()/(double) RAND_MAX) ;//* pow(-1,i + j);
	complex->set_value(i,j,REAL,r); 
	complex->set_value(i,j,IMAG,im);
      }
    }
  }
}

void FCDI_IllumRecon::set_support(double ** object_support){
  for(int i=0; i< complex->get_size_x(); i++){
    for(int j=0; j< complex->get_size_y(); j++){
      support[i][j] = object_support[i][j];
    }
  }
}


void FCDI_IllumRecon::set_intensity(double ** detector_intensity){
  for(int i=0; i< complex->get_size_x(); i++){
    for(int j=0; j< complex->get_size_y(); j++){
      intensity_sqrt[i][j] = detector_intensity[i][j]; //sqrt(detector_intensity[i][j]);
    }
  }
}

/**double FCDI_IllumRecon::calculate_error(bool error){
  do_error = error;
  }**/

double FCDI_IllumRecon::get_error(){
    
    return current_error;
}


int FCDI_IllumRecon::iterate(){
  
  //start with the guess in the zone plate plane.

  //assume the wavefield at the zone plate is initalised to Pexp(...)

  int nx = complex->get_size_x();
  int ny = complex->get_size_y();
  
  double ** result = new double*[nx];
  for(int i=0; i < nx; i++)
    result[i]= new double[ny];

  complex->get_2d(MAG,&result);
  write_ppm("before.ppm", nx, ny, result);
  complex->get_2d(PHASE,&result);
  write_ppm("before_phase.ppm", nx, ny, result);
  
  //cout << "0: "<< complex->get_norm() << endl;

  //propogate to the focal plane.
  fft->perform_forward_fft(complex);
  //complex->invert();

  /**for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      double temp = complex->get_real(i,j);
      complex->set_real(i,j,complex->get_imag(i,j));
      complex->set_imag(i,j,-1*temp);
    }
    }**/

  complex->get_2d(MAG,&result);
  write_ppm("A_forward_mag.ppm", nx, ny, result);
  complex->get_2d(PHASE,&result);
  write_ppm("A_forward_phase.ppm", nx, ny, result);

  complex->multiply(A_forward);
  
  complex->get_2d(MAG,&result);
  write_ppm("1.ppm", nx, ny, result);
  cout << "1: "<< complex->get_norm() << endl;

  complex->multiply(B_forward);
  complex->get_2d(MAG,&result);
  write_ppm("B_forward_mag.ppm", nx, ny, result);
  complex->get_2d(PHASE,&result);
  write_ppm("B_forward_phase.ppm", nx, ny, result);

  //multiply by A and B and constants (from nature paper 2006)
  //complex->multiply(forward_coefficients,(1.0/zone_to_focal_length));
  //complex->multiply(forward_coefficients,(1.0/zone_to_focal_length));

  //forward_coefficients->get_2d(MAG,&result);
  //write_ppm("1.forward.ppm", nx, ny, result);
  //cout << "forward_norm: "<< forward_coefficients->get_norm() << endl;

  //complex->get_2d(MAG,&result);
  //write_ppm("1.5.ppm", nx, ny, result);

  //propogate to the detector plane.
  fft->perform_forward_fft(complex);

  /** for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      double temp = complex->get_real(i,j);
      complex->set_real(i,j,complex->get_imag(i,j));
      complex->set_imag(i,j,-1*temp);
    }
    }**/
  

  complex->multiply(forward_coefficients_const);

  complex->get_2d(MAG,&result);
  write_ppm("2.ppm", nx, ny, result);
  cout << "2: "<< complex->get_norm() << endl;

  complex->get_2d(PHASE,&result);
  write_ppm("phase_bf.ppm", nx, ny, result);

  //impose the amplitude constraint
  scale_intensity(complex);
  cout << "Current Error is: " << current_error<<endl;

  complex->get_2d(MAG,&result);
  write_ppm("3.ppm", nx, ny, result);
  cout << "3: "<< complex->get_norm() << endl;

  //complex->get_2d(PHASE,&result);
  //write_ppm("phase_after.ppm", nx, ny, result);

  //go back to the focal plane

  fft->perform_forward_fft(complex);

  //complex->multiply(B_backward);

  /**
  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      double temp = complex->get_real(i,j);
      complex->set_real(i,j,complex->get_imag(i,j));
      complex->set_imag(i,j,-1*temp);
    }
    }**/

  complex->multiply(A_backward);

  //complex->multiply(forward_coefficients_const);


  complex->get_2d(MAG,&result);
  write_ppm("4.ppm", nx, ny, result,log);
  cout << "4: "<< complex->get_norm() << endl;
  complex->get_2d(PHASE,&result);
  write_ppm("4_phase.ppm", nx, ny, result);

  //multiply by constants
  //complex->multiply(backward_coefficients,1.0/focal_to_detector_length);

  //complex->get_2d(MAG,&result);
  //write_ppm("4.5.ppm", nx, ny, result);

  complex->multiply(B_backward);

  //go back to zone plate plane. 

  //complex->multiply(A_backward);

  fft->perform_forward_fft(complex);
 
  /**for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      double temp = complex->get_real(i,j);
      complex->set_real(i,j,complex->get_imag(i,j));
      complex->set_imag(i,j,-1*temp);
    }
    }**/

  
  complex->multiply(backward_coefficients_const);

  complex->get_2d(MAG,&result);
  write_ppm("5.ppm", nx, ny, result);
  cout << "5: "<< complex->get_norm() << endl;

  //apply support constraint
  project_support(complex);

  complex->get_2d(MAG,&result);
  write_ppm("6.ppm", nx, ny, result);
  cout << "6: "<< complex->get_norm() << endl;

  complex->get_2d(PHASE,&result);
  write_ppm("6_phase.ppm", nx, ny, result);

  for(int i=0; i < nx; i++)
    delete [] result[i];
  delete [] result;

  return SUCCESS;
}

void FCDI_IllumRecon::project_support(Complex_2D * c){
  int nx = c->get_size_x();
  int ny = c->get_size_y();

  for(int i=0; i< nx; ++i){
    for(int j=0; j< ny; ++j){
      double mag = c->get_mag(i,j);
      //c->set_real(i,j,mag);
      //c->set_imag(i,j,0);
      if(support[i][j]==0){
	c->set_real(i,j,0);
	c->set_imag(i,j,0);
      }
    }
  }
}


void FCDI_IllumRecon::project_intensity(Complex_2D * c){
  fft->perform_forward_fft(c);    
  scale_intensity(c);
  fft->perform_backward_fft(c);  
}


void FCDI_IllumRecon::scale_intensity(Complex_2D * c){

  double norm2_mag=0;
  double norm2_diff=0;

  int nx = c->get_size_x();
  int ny = c->get_size_y();

  for(int i=0; i< nx; ++i){
    for(int j=0; j< ny; ++j){
      //scale
      double current_mag = c->get_mag(i,j);
      if(current_mag > 0.0)
	c->set_mag(i,j,intensity_sqrt[i][j]/current_mag);

      //calculate the error
      norm2_mag += intensity_sqrt[i][j]*intensity_sqrt[i][j];
      norm2_diff += pow(current_mag-intensity_sqrt[i][j],2);
    }
  }
  current_error = (norm2_diff/norm2_mag);
}

