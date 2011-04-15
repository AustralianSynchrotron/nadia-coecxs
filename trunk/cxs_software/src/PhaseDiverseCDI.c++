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
#include "PlanarCDI.h"
#include "PhaseDiverseCDI.h"
#include "io.h" //
#include <sstream>
#include <typeinfo>

using namespace std;


void get_center(Double_2D & mag, int * x, int * y){
  
  double offset=0;
  double mean_x = 0;
  double mean_y = 0;
  double total = 0;

  for(int i=0; i<mag.get_size_x(); i++){
    for(int j=0; j<mag.get_size_y(); j++){

      mean_x += fabs(mag.get(i,j)-offset)*i;
      mean_y += fabs(mag.get(i,j)-offset)*j;
      
      total += fabs(mag.get(i,j)-offset);

    }
  }

  *x = mean_x/total;
  *y = mean_y/total;


}




PhaseDiverseCDI::PhaseDiverseCDI(Complex_2D & object,
				 double beta, 
				 double gamma, 
				 bool parallel,
				 int granularity):object(object), 
						  scale(granularity),
						  beta(beta), 
						  gamma(gamma), 
						  parallel(parallel),
						  nx(object.get_size_x()),
						  ny(object.get_size_y()){
  
  iterations_per_cycle = 1;
  weight_norm=0;

};


PhaseDiverseCDI::~PhaseDiverseCDI(){

  while(single_result.size()!=0){
    Complex_2D * temp = single_result.back();
    single_result.pop_back();
    delete temp;
  }

  if(!weight_norm)
    delete weight_norm;

}

void PhaseDiverseCDI::add_new_position(PlanarCDI * local,
				       double x, 
				       double y,
				       double alpha){

  int nx = local->get_size_x() ;
  int ny = local->get_size_y() ;
    
  singleCDI.push_back(local);
  x_position.push_back(x);
  y_position.push_back(y);
  this->alpha.push_back(alpha);
  single_result.push_back(new Complex_2D(nx,ny));

  cout << "Added position "<<singleCDI.size()-1<<endl;
  
}


void PhaseDiverseCDI::set_up_weight_norm(){

  if(!weight_norm)
    weight_norm=new Double_2D(nx,ny);

  for(int n=0; n<singleCDI.size(); n++){
    
    double x = x_position.at(n);
    double y = y_position.at(n);

    int this_size_x = single_result.at(n)->get_size_x();
    int this_size_y = single_result.at(n)->get_size_y();

    Double_2D this_single_weights(this_size_x,this_size_y);
    get_weights(n,this_single_weights);
    
    for(int i_=0; i_< this_size_x; i_++){
      for(int j_=0; j_< this_size_y; j_++){
	
	double i = (i_-x)*scale;
	double j = (j_-y)*scale;
      
	if(i>0&&j>0&&i<nx&&j<ny){

	  double new_weight = this_single_weights.get(i_,j_);
	  double old_weight = weight_norm->get(i,j);

	  if(n==0)
	    weight_norm->set(i,j,new_weight);
	  else
	    weight_norm->set(i,j,new_weight+old_weight);
	}
      }
    }
  }
};

void PhaseDiverseCDI::initialise_estimate(){
  
  for(int i=0; i<object.get_size_x(); i++){
    for(int j=0; j<object.get_size_y(); j++){
      object.set_real(i,j,1);
    }
  }

  for(int i=0; i<singleCDI.size(); i++){
    add_to_object(i,1,0);
  }
  
}


void PhaseDiverseCDI::iterate(){

  static int total_iterations = 0;

  cout << "Iteration "<<total_iterations<<endl;

  bool parallel = true;  
    
  Double_2D result(1024,1024);

  for(int i=0; i<singleCDI.size(); i++){

    int x, y;

    //if(i!=0 && (total_iterations==2 || total_iterations==4))
    //  check_position(i);
    update_from_object(i);

    for(int j=0; j<iterations_per_cycle; j++){
      singleCDI.at(i)->iterate();
      cout << "This Error="<<singleCDI.at(i)->get_error()<<endl;
    }

    if(!parallel)
      add_to_object(i,beta,0);

  }
  
  if(parallel){
    for(int i=0; i<singleCDI.size(); i++){
      add_to_object(i,beta,0);
    }
  }

  total_iterations++;

};

void PhaseDiverseCDI::get_result(PlanarCDI * local, Complex_2D & result){
  if(typeid(*local)==typeid(FresnelCDI)){
    ((FresnelCDI*)local)->get_transmission_function(result);
  }
  else
    result.copy(local->get_exit_surface_wave());
};

void PhaseDiverseCDI::set_result(PlanarCDI * local, Complex_2D & result){
  
  if(typeid(*local)==typeid(FresnelCDI)){
    ((FresnelCDI*)local)->set_transmission_function(result);
  }
  else
    local->set_exit_surface_wave(result);
  
}



int PhaseDiverseCDI::check_position(int n_probe, double step_size, int tries){
  
  double x = x_position.at(n_probe);
  double y = y_position.at(n_probe);

  cout << "checking probe "<< n_probe << " position, " 
       << x << "," << y << " with step size " 
       << step_size << ". Try no.: "<<tries<< endl;

  //done, we found the best position.
  if(step_size<1.0/scale)
    return SUCCESS;

  //failed, we moved around a bit, but couldn't find a local minima
  //in the error metric.
  if(tries>10){
    cout << "giving up on probe "<< n_probe << ". Coule not find " 
	 << "it's position. Returning to the original. " 
	 << endl;
      return FAILURE;
  }
  
  PlanarCDI * single = singleCDI.at(n_probe);
  
  double best_x=x;
  double best_y=y;
  double best_error=100;
   
  double new_x;
  double new_y;
  double size_x = single_result.at(n_probe)->get_size_x();
  double size_y = single_result.at(n_probe)->get_size_y();
  
  //try the 9 positions around the current one
  //record the one with the lowest error metric.
  for(int i=-1; i<2; i++){
    for(int j=-1; j<2; j++){

      new_x=x+i*step_size;
      new_y=y+j*step_size;

      //checking whether the corrds are still within the 
      //boundary of the image
      //      if(new_x>0 && new_x<size_x && new_y>0 && new_y<size_y){

	x_position.at(n_probe)=new_x;
	y_position.at(n_probe)=new_y;
	
	update_from_object(n_probe);
	set_result(single,*(single_result.at(n_probe)));
	
	single->iterate();
	
	if(single->get_error()<best_error){
	  best_error = single->get_error();
	  best_x = x+i*step_size;
	  best_y = y+j*step_size;
	}
	  // }
    }
  }
  

  //set x and y to the best one :
  x_position.at(n_probe)=best_x;
  y_position.at(n_probe)=best_y;

  //recursively call this function with a smaller step size.
  if(best_x==x && best_y==y)
    step_size=step_size/2.0;
  
  int status = check_position(n_probe, step_size, ++tries);

  if(status == FAILURE){
    //return to orginal coordinates
    x_position.at(n_probe) = x ;
    y_position.at(n_probe) = y;
    return FAILURE;
  }

  
  cout << "moving probe "<< n_probe << " by " 
       << x_position.at(n_probe)-x <<" in x "
       << "and "<< y_position.at(n_probe)-y<<" in y." 
       << endl;
  
  return SUCCESS;
  
}

void PhaseDiverseCDI::get_weights(int n_probe, Double_2D & weights){

  int nx = single_result.at(n_probe)->get_size_x();
  int ny = single_result.at(n_probe)->get_size_y();
  PlanarCDI * this_CDI = singleCDI.at(n_probe);
  double value;
  Double_2D illum_mag(nx,ny);
  double max;

  if(typeid(*this_CDI)==typeid(FresnelCDI)){
    const Complex_2D & illum = ((FresnelCDI*)(this_CDI))->get_illumination_at_sample();
    illum.get_2d(MAG,illum_mag);
    max = illum_mag.get_max();
  }

  for(int i=0; i< nx; i++){
    for(int j=0; j< ny; j++){
      
      if(typeid(*this_CDI)==typeid(FresnelCDI))
	value = illum_mag.get(i,j)/max;
      else
	value = 1;

      value=alpha.at(n_probe)*pow(value,gamma);
      weights.set(i,j,value);
    }
  }
  
}


void PhaseDiverseCDI::add_to_object(int n_probe, double new_fraction, 
				    double old_fraction){
  if(!weight_norm)
    set_up_weight_norm();

  get_result(singleCDI.at(n_probe),*(single_result.at(n_probe)));
  
  double x_offset = x_position.at(n_probe);
  double y_offset = y_position.at(n_probe);
  Complex_2D & small = *(single_result.at(n_probe));
  //Double_2D & beam = *beam_shape.at(n_probe);
  
  //round off to a first approx. Will be fixed later.
  int i, j; //local coors (i,j) are global (object) coors

  Double_2D temp_weights(small.get_size_x(),
			 small.get_size_y());

  get_weights(n_probe, temp_weights);

  //using the local coordinate system
  for(int i_=0; i_< small.get_size_x(); i_++){
    for(int j_=0; j_< small.get_size_y(); j_++){
      
      i = (i_-x_offset)*scale;
      j = (j_-y_offset)*scale;
      
      if(i>0&&j>0&&i<nx&&j<ny){
	double weight = new_fraction*temp_weights.get(i_,j_)/weight_norm->get(i,j);
      
	if(weight_norm->get(i,j)<=0){
	  weight=0;
	}

	if(weight>1){
	  cout << weight << endl;
	}

	double new_real = weight*small.get_real(i_,j_)
	  +(1-weight)*object.get_real(i,j);
	
	double new_imag = weight*small.get_imag(i_,j_)
	  +(1-weight)*object.get_imag(i,j);
	
	update_to_object_sub_grid(i,j,new_real,new_imag); 
      }
      
    }
  } 
  
}

void PhaseDiverseCDI::update_to_object_sub_grid(int i, int j, 
			       double real_value, 
			       double imag_value){
  
  for(int di=0; di < scale; di++){
    for(int dj=0; dj < scale; dj++){
      object.set_real(i+di,j+dj,real_value);
      object.set_imag(i+di,j+dj,imag_value);
    }
  }
}

void PhaseDiverseCDI::update_from_object_sub_grid(int i, int j, 
				 double & real_value, 
				 double & imag_value){
  
  real_value=0;
  imag_value=0;

  for(int di=0; di < scale; di++){
    for(int dj=0; dj < scale; dj++){
      real_value+=object.get_real(i+di,j+dj);
      imag_value+=object.get_imag(i+di,j+dj);
    }
  }

  real_value=real_value/(scale*scale);
  imag_value=imag_value/(scale*scale);

}


void PhaseDiverseCDI::update_from_object(int n_probe){

  double x_offset = x_position.at(n_probe);
  double y_offset = y_position.at(n_probe);
  Complex_2D * small = single_result.at(n_probe);
  //  Double_2D * beam = beam_shape.at(n_probe);
  
  //round off to a first approx. Will be fixed later.
  int i, j; //local coors (i,j) are global (object) coors

  //using the local coordinate system
  for(int i_=0; i_< small->get_size_x(); i_++){
    for(int j_=0; j_< small->get_size_y(); j_++){

      i = (i_-x_offset)*scale;
      j = (j_-y_offset)*scale;
      
      //      if(beam->get(i_,j_)>0 && i>=0 && j>=0 && 
      //	 i<object.get_size_x() && j<object.get_size_y()){

      if(i>=0 && j>=0 && i<nx && j<ny){

	double new_real; //= object.get_real(i,j);
	double new_imag; //= object.get_imag(i,j);
	
	update_from_object_sub_grid(i, j, 
				    new_real, 
				    new_imag);

	small->set_real(i_,j_,new_real);
	small->set_imag(i_,j_,new_imag);
      }
    }
  }

  set_result(singleCDI.at(n_probe),*(single_result.at(n_probe)));

};

