// Copyright 2012 T'Mir Danger Julius for The ARC Centre of Excellence in 
// Coherent X-ray Science. This program is distributed under the GNU  
// General Public License. We also ask that you cite this software in 
// publications where you made use of it for any part of the data     
// analysis. 

#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <cstdlib> 
#include "Complex_2D.h"
#include "Double_2D.h"
#include "PartialCDI.h"
#include "io.h" 
#include <sstream>
#include <typeinfo>
#include "utils.h"

using namespace std;

//constructor for the class which handles partially coherent diffractive imaging.

PartialCDI::PartialCDI(Complex_2D & initial_guess,
    double beta, 
    double lcx,
    double lcy,
    double lx,
    double ly,
    int n_best,
    bool parallel
    )
:BaseCDI(initial_guess,n_best),
  beta(beta), 
  lcx(lcx),
  lcy(lcy),
  parallel(parallel),
  transmission(nx,ny),
  threshold(-1.0),
  iterations_per_cycle(1){

    x_position.clear();
    y_position.clear();

    //Create arrays for the x and y positions of pixels within the detector.

    pxsize=lx/nx;
    pysize=ly/ny;
    
    double x_0=(nx/2*(-pxsize))+(pxsize/2);
    double y_0=(ny/2*(-pysize))+(pysize/2);

    x_position.push_back(x_0);
    y_position.push_back(y_0);

    for(int i=1; i<nx; i++){
      x_position.push_back(x_position.at(i-1)+pxsize);
    }
    for(int i=1; i<ny; i++){
      y_position.push_back(y_position.at(i-1)+pysize);
    }

  }


//destructor for cleaning up
PartialCDI::~PartialCDI(){

}

//return the transmision function.
Complex_2D  PartialCDI::get_transmission(){
  return transmission;
}

//set the transmission function 
//this function can be used to set the initial estimate.
void PartialCDI::set_transmission(Complex_2D & new_transmission){

  int new_nx = new_transmission.get_size_x();
  int new_ny = new_transmission.get_size_y();

  if(new_nx!=nx||new_ny!=ny){
    cout << "Can not set the transmission function in "
      << "PartialCDI because the complex array "
      << "given does not have the same dimensions as "
      << "the current transmission function"<<endl;
    return;
  }

  transmission.copy(new_transmission);

}

//set the initial guess.
//fills the transmission function with a random 
//value between 0 and 1 within the support, and
//0 everywhere else. Complex is then set to be 
//the value of the transmission function
void PartialCDI::initialise_estimate(int seed){

  //initialise the random number generator
  srand(seed);

  int max_value = intensity_sqrt.get_max();

  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){

      if(!support.get(i,j)){ 

	transmission.set_value(i,j,REAL,0);
	transmission.set_value(i,j,IMAG,0);
      }
      else{
	double re = support.get(i,j)*(max_value*rand()/(double) RAND_MAX);
	double im = support.get(i,j)*(max_value*rand()/(double) RAND_MAX);

	transmission.set_value(i,j,REAL,re);
	transmission.set_value(i,j,IMAG,im);
      }
    }
  }
  complex=transmission;
}

//this iterate function overrides that of BaseCDI, 
//and handles the multiple modes.
int PartialCDI::iterate(){

  //below is the code for the special case of ER
  //this is faster than using the generic algorithm code
  //further down in this function.
  
  vector<Complex_2D> temp_complex_PFS;
  vector<Complex_2D> temp_complex_PF;
  vector<Complex_2D> temp_complex_PS;
  vector<Complex_2D> temp_complex_PSF;

  singleCDI.clear();
  for(unsigned int mode=0; mode<singlemode.size(); mode++){
    singleCDI.push_back(singlemode.at(mode));
    apply_transmission(singleCDI.at(mode));
  }

  if(algorithm==ER){
    for(int mode=0; mode<singleCDI.size(); mode++){
      propagate_to_detector(singleCDI.at(mode));
    }
    scale_intensity(singleCDI);
    propagate_from_detector(singleCDI.back());
    apply_support(singleCDI.back());
    update_transmission();

    return SUCCESS;
  }

  //start of the generic algorithm code

  //PFS
  if(algorithm_structure[PFS]!=0){
    for(int mode=0; mode<singleCDI.size(); mode++){
      temp_complex_PFS.push_back(singleCDI.at(mode));
      apply_support(temp_complex_PFS.at(mode));
      propagate_to_detector(temp_complex_PFS.at(mode));
    } 
    scale_intensity(temp_complex_PFS);
    propagate_from_detector(temp_complex_PFS.back());
  }

  //F
  if(algorithm_structure[PF]!=0){

    for(int mode=0; mode<singleCDI.size(); mode++){
      temp_complex_PF.push_back(singleCDI.at(mode));
      propagate_to_detector(temp_complex_PF.at(mode));
    }

    scale_intensity(temp_complex_PF);
    propagate_from_detector(temp_complex_PF.back());
  }

  //S
  if(algorithm_structure[PS]!=0){
    temp_complex_PS.push_back(singleCDI.back());
    apply_support(temp_complex_PS.back());
  }

  //SF
  if(algorithm_structure[PSF]!=0){
    if(algorithm_structure[PF]!=0){
      temp_complex_PSF.push_back(temp_complex_PF.back());
      apply_support(temp_complex_PSF.back());
    }else{
      for(int mode=0; mode<singleCDI.size(); mode++){
	temp_complex_PSF.push_back(singleCDI.at(mode));
	propagate_to_detector(temp_complex_PSF.at(mode));
      }           
      scale_intensity(temp_complex_PSF);
      propagate_from_detector(temp_complex_PSF.back());
      apply_support(temp_complex_PSF.back());
    }
  }

  //combine the result of the separate operators
  double value_real, value_imag;
  for(int i=0; i < nx; ++i){
    for(int j=0; j < ny; ++j){

      //Add the identity
      value_real = (1+algorithm_structure[PI])*singleCDI.back().get_real(i,j);
      value_imag = (1+algorithm_structure[PI])*singleCDI.back().get_imag(i,j);

      //Add the component from the PfPs operator
      if(algorithm_structure[PFS]!=0){
	value_real+=algorithm_structure[PFS]*temp_complex_PFS.back().get_real(i,j);
	value_imag+=algorithm_structure[PFS]*temp_complex_PFS.back().get_imag(i,j);
      }

      //Add the component from the Pf operator
      if(algorithm_structure[PF]!=0){
	value_real+=algorithm_structure[PF]*temp_complex_PF.back().get_real(i,j);
	value_imag+=algorithm_structure[PF]*temp_complex_PF.back().get_imag(i,j);

      }

      //Add the component from the Ps operator
      if(algorithm_structure[PS]!=0){
	value_real+=algorithm_structure[PS]*temp_complex_PS.back().get_real(i,j);
	value_imag+=algorithm_structure[PS]*temp_complex_PS.back().get_imag(i,j);
      }

      //Add the component from the PsPf operator
      if(algorithm_structure[PSF]!=0){
	value_real+=algorithm_structure[PSF]*temp_complex_PSF.back().get_real(i,j);
	value_imag+=algorithm_structure[PSF]*temp_complex_PSF.back().get_imag(i,j);
      }

      singleCDI.back().set_real(i,j,value_real);
      singleCDI.back().set_imag(i,j,value_imag);

    }
  }

  //Update the transmission using the dominant mode
  update_transmission();

  return SUCCESS;
}

//uses the complex_2d multiply function to apply 
//the transmission function
void PartialCDI::apply_transmission(Complex_2D & c){
    c.multiply(transmission, 1);
}

//scale the highest occupancy mode 
//this overwrites the function of the same
//name in BaseCDI
void PartialCDI::scale_intensity(vector<Complex_2D> & c){
  double norm2_mag=0;
  double norm2_diff=0;
  double measured_int=0;
  double calc_int=0;
  double scaled_mag=0;

  Double_2D magnitude=sum_intensity(c);

  //set the intensities to 0
/*  for(int i=0; i< nx; i++){
    for(int j=0; j< ny; j++){
      magnitude.set(i,j,0.0);
    }
  }

  for(unsigned int mode=0; mode<c.size(); mode++){

    for(int i=0; i< nx; i++){
      for(int j=0; j< ny; j++){

	double mag_total = magnitude.get(i,j)+eigen.at(mode)*c.at(mode).get_mag(i,j)*c.at(mode).get_mag(i,j);
	magnitude.set(i,j,mag_total);
      }
    }
  }
*/
  for(int i=0; i< nx; i++){
    for(int j=0; j< ny; j++){
      //reset the magnitude
      if(beam_stop==0 || beam_stop->get(i,j)>0){

	measured_int=intensity_sqrt.get(i,j);

	calc_int=sqrt(magnitude.get(i,j)); 

	scaled_mag=c.back().get_mag(i,j)*measured_int/calc_int;
	c.back().set_mag(i,j,scaled_mag);

	//calculate the error
	norm2_mag += measured_int*measured_int;

	norm2_diff += (measured_int-calc_int)*(measured_int-calc_int);
      }
    }
  }
  current_error = norm2_diff/norm2_mag;

}


//add the intensities across all modes scaled by the eigenvalues
Double_2D PartialCDI::sum_intensity(vector<Complex_2D> & c){

  Double_2D magnitude(nx, ny);

  //set the intensities to 0
/*  for(int i=0; i< nx; i++){
    for(int j=0; j< ny; j++){
      magnitude.set(i,j,0.0);
    }
  }*/


  for(int i=0; i< nx; i++){
    for(int j=0; j< ny; j++){

      double mag_total=0;
      for(int mode=0; mode<c.size(); mode++){

	mag_total += eigen.at(mode)*c.at(mode).get_mag(i,j)*c.at(mode).get_mag(i,j);

      }
      magnitude.set(i,j,mag_total);
    }
  }

  return(magnitude);
};


//calculate the transmission function by dividing 
//the highest occupacy mode at the source by the 
//highest occupancy mode at the detector
void PartialCDI::update_transmission(){

  double transreal, transimag, rd, id, rs, is;

  for(int i=0; i< nx; i++){
    for(int j=0; j< ny; j++){

      rd = singlemode.back().get_real(i,j);
      id = singlemode.back().get_imag(i,j);
      rs = singleCDI.back().get_real(i,j);
      is = singleCDI.back().get_imag(i,j);

      if(abs(rd*rd+id*id) > 0.0000000001){ 

	transreal = (rd*rs+id*is)/(rd*rd+id*id);
	transimag = (rd*is-rs*id)/(rd*rd+id*id);

	transmission.set_real(i, j, transreal);
	transmission.set_imag(i, j, transimag);

      }else{
	transmission.set_real(i, j, 1);
	transmission.set_imag(i, j, 1);
      }
    }
  }

  complex=transmission;
}

//generate the S and J matrices for the decomposition 
//of the partially coherent wave where JC=nSC where
//H = integral(P*l(r1)J(r1, r2)Pm(r2)) dr1 dr2 and 
//S=integral(P*l(r)pm(r))dr where Pl is an orhtonormal
//basis set, in this case, the Legendre polynomials
void PartialCDI::initialise_matrices(int leg, int modes){

  nleg = leg;
  nmode = modes;

  Double_2D roots = legroots(nleg);

  vector<double> rootval;


  for(int i=0; i< roots.get_size_x(); i++){
    rootval.push_back(roots.get(i, 0));
  }

  Double_2D legmatrix = fill_legmatrix(rootval, nmode); 

  fill_jmatrix(legmatrix, roots);
  fill_smatrix(legmatrix, roots);
  solve_gep(*jmatrix, *smatrix, eigen);
  fill_modes(*jmatrix);

  Complex_2D mags(nx, ny);

  //set the intensities to 0
  for(int i=0; i< nx; i++){
    for(int j=0; j< ny; j++){
      mags.set_real(i,j,0.0);
      mags.set_imag(i, j, 0.0);
    }
  }


  for(int i=0; i< nx; i++){
    for(int j=0; j< ny; j++){

      double mag_total=0;
      for(unsigned int mode=0; mode<singlemode.size(); mode++){
	mag_total += eigen.at(mode)*singlemode.at(mode).get_mag(i,j)*singlemode.at(mode).get_mag(i,j);
      }
      mags.set_real(i,j,mag_total);
    }
  }

  complex=mags;
  return;
}

//the S matrix = integral(P*l(r)pm(r))dr = 2/(2n+1)
//from the orthogonality requirments of Legendre 
//Polynomials. We then turn it in to a 2D matrix
//for the x and y dimensions
void PartialCDI::fill_smatrix(Double_2D legmatrix, Double_2D roots){

  Complex_2D s1d(nleg,nleg);

  smatrix = new Complex_2D(nmode*nmode, nmode*nmode);

  for(int i = 0; i < nmode; i++){
    for(int j = 0; j < nmode; j++){
      //      if(i==j){
      //	s1d.set_real(i, j, 1.0);//2.0/(2.0*i+1.0));
      double s1_real = 0;
      for(int k=0; k<nleg; k++){
	s1_real+= roots.get(k,1)*legmatrix.get(k,i)*legmatrix.get(k,j)/sqrt(2.0/(2*j+1))/sqrt(2.0/(2*i+1));
      }
      if(fabs(s1_real)>1.0e-5){
	s1d.set_real(i, j, s1_real);
      }else{
	s1d.set_real(i, j, 0);
      }
      s1d.set_imag(i, j, 0);
    }
    }

    for(int i=0; i<nmode; i++){
      for(int j=0; j<nmode; j++){
	for(int k=0; k<nmode; k++){
	  for(int l=0; l<nmode; l++){

	    double val=s1d.get_real(i,k)*s1d.get_real(j,l);
	    smatrix->set_real(i*nmode+j,k*nmode+l, val);
	    smatrix->set_imag(i*nmode+j,k*nmode+l, 0.0);

	  }
	}
      }
    }

/*    for(int i=0; i<nmode*nmode; i++){
      for(int j=0; j<nmode*nmode; j++){
	std::cout<<smatrix->get_real(i, j)<<" ";
      }
      std::cout<<"\n";
    }*/
    return;
  }

  //the J matrix where J = integral(P*l(r1)J(r1, r2)Pm(r2))dr1dr2 
  //the x and y are computed seperately, then multiplied together
  //to take the 1D matrix to the 2D matrix.
  void PartialCDI::fill_jmatrix(Double_2D legmatrix, Double_2D roots){

    Complex_2D xjmatrix(nmode, nmode);
    Complex_2D yjmatrix(nmode, nmode);

    double xj_real;
    double yj_real;
    double xj_imag;
    double yj_imag;

    double scalex;
    double scaley;


    for(int i=0; i<nmode; i++){
      for(int j=0; j<nmode; j++){

	xj_real=0.0;
	yj_real=0.0;
	xj_imag=0.0;
	yj_imag=0.0;

	for(int k=0; k<nleg; k++){
	  for(int l=0; l<nleg; l++){

	    scalex=exp(-(1/(2*lcx*lcx))*(roots.get(k,0)-roots.get(l,0))*(roots.get(k,0)-roots.get(l,0)));
	    scaley=exp(-(1/(2*lcy*lcy))*(roots.get(k,0)-roots.get(l,0))*(roots.get(k,0)-roots.get(l,0)));

	    xj_real+= scalex*roots.get(k,1)*roots.get(l,1)*legmatrix.get(k,i)*legmatrix.get(l,j)/sqrt(2.0/(2*j+1))/sqrt(2.0/(2*i+1));
	    yj_real+= scaley*roots.get(k,1)*roots.get(l,1)*legmatrix.get(k,i)*legmatrix.get(l,j)/sqrt(2.0/(2*j+1))/sqrt(2.0/(2*i+1));

	    xj_imag+=0;
	    yj_imag+=0;

	  }
	}

	xjmatrix.set_real(i,j, xj_real);
	xjmatrix.set_imag(i,j, xj_imag);
	yjmatrix.set_real(i,j, yj_real);
	yjmatrix.set_imag(i,j, yj_imag);

      }
    }

    jmatrix = new Complex_2D(nmode*nmode, nmode*nmode);

    for(int i=0; i<nmode; i++){
      for(int j=0; j<nmode; j++){
	for(int k=0; k<nmode; k++){
	  for(int l=0; l<nmode; l++){

	    double val_real=xjmatrix.get_real(i,k)*yjmatrix.get_real(j,l);
	    double val_imag=xjmatrix.get_imag(i,k)*yjmatrix.get_imag(j,l);

	    jmatrix->set_real(i*nmode+j,k*nmode+l, val_real);
	    jmatrix->set_imag(i*nmode+j,k*nmode+l, val_imag);

	  }
	}
      }
    }
    return;
  }

  //fill a vector of Complex_2D for single modes. These 
  //modes do not evolve over time, and so are not BaseCDI's
  void PartialCDI::fill_modes(Complex_2D & c){

    Double_2D x_legmatrix = fill_legmatrix(x_position, nmode);
    Double_2D y_legmatrix = fill_legmatrix(y_position, nmode);

    double val_imag; 
    double val_real;

    singlemode.clear();

    int e_mode=0;

    for(int mode=0; mode<nmode*nmode; mode++){

      if(eigen.at(e_mode)/eigen.back() > threshold){

	Complex_2D tmp(nx, ny);

	for(int i=0; i<nx; i++){
	  for(int j=0; j<ny; j++){

	    val_imag=0;
	    val_real=0;
	    for(int k=0; k<nmode; k++){
	      for(int l=0; l<nmode; l++){

		val_real +=c.get_real(l+nmode*k, mode)*x_legmatrix.get(i, k)*y_legmatrix.get(j, l);
		val_imag +=c.get_imag(l+nmode*k, mode)*x_legmatrix.get(i, k)*y_legmatrix.get(j, l);

	      }
	    }
	    tmp.set_real(i, j, val_real);
	    tmp.set_imag(i, j, val_imag);
	  }
	}
	singlemode.push_back(tmp);

	e_mode++;
      }else{
	eigen.erase(eigen.begin()+e_mode);
      }
    }
  }

  //Propagate from the object plane to the detector
  void PartialCDI::propagate_to_detector(Complex_2D & c){

    c.perform_forward_fft();
    c.invert(true);

  }

  //Propagate form the detector plane to the object
  void PartialCDI::propagate_from_detector(Complex_2D & c){

    c.invert(true);
    c.perform_backward_fft();

  }

  //find and return the current intensity
  Double_2D PartialCDI::get_intensity(){

    Double_2D magnitude(nx, ny);

    singleCDI.clear();
    for(unsigned int i=0; i<singlemode.size(); i++){
      singleCDI.push_back(singlemode.at(i));
      apply_transmission(singleCDI.at(i));
    }

    //set the intensities to 0
    for(int i=0; i< nx; i++){
      for(int j=0; j< ny; j++){
	magnitude.set(i,j,0.0);
      }
    }

    for(unsigned int mode=0; mode<singleCDI.size(); mode++){

      for(int i=0; i< nx; i++){
	for(int j=0; j< ny; j++){

	  double mag_total = magnitude.get(i,j)+eigen.at(mode)*singleCDI.at(mode).get_mag(i,j)*singleCDI.at(mode).get_mag(i,j);
	  magnitude.set(i,j,mag_total);
	}
      }
    }

    return(magnitude);
  }


  //propagate each mode to the detector. This is specifically
  //for use in simulations.
  Double_2D PartialCDI::propagate_modes_to_detector(){

    singleCDI.clear();
    for(unsigned int i=0; i<singlemode.size(); i++){
      singleCDI.push_back(singlemode.at(i));
    }

    for(int mode=0; mode<singleCDI.size(); mode++){
      apply_transmission(singleCDI.at(mode));
      propagate_to_detector(singleCDI.at(mode));
    }

    Double_2D magnitude(nx, ny);
    for(unsigned int mode=0; mode<singleCDI.size(); mode++){

      for(int i=0; i< nx; i++){
	for(int j=0; j< ny; j++){

	  double mag_total = magnitude.get(i,j)+eigen.at(mode)*singleCDI.at(mode).get_mag(i,j)*singleCDI.at(mode).get_mag(i,j);
	  magnitude.set(i,j,mag_total);
	}
      }
    }
    return(magnitude);
  }

  //a function that returns a single mode.
  Complex_2D PartialCDI::get_mode(int mode){
    if(mode>=singlemode.size()){
      std::cout<<"Maximum mode is "<<singlemode.size()<<". Returning this mode instead\n";
      return(singlemode.back());
    }
    return(singlemode.at(mode));
  }


