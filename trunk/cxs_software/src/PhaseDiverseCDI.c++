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
#include "utils.h"

using namespace std;


//constructor for the class which handles phase diversity.
PhaseDiverseCDI::PhaseDiverseCDI(
				 double beta, 
				 double gamma, 
				 bool parallel,
				 int granularity):scale(granularity),
						  beta(beta), 
						  gamma(gamma), 
						  parallel(parallel),
						  nx(0),
						  ny(0),
						  x_min(0),
						  y_min(0),
						  weights_set(false),
                                                  iterations_per_cycle(1),
                                                  object(0){
  
};


//destructor for cleaning up
PhaseDiverseCDI::~PhaseDiverseCDI(){

  while(single_result.size()!=0){
    Complex_2D * temp = single_result.back();
    single_result.pop_back();
    delete temp;

    Double_2D * temp_d = weights.back();
    weights.pop_back();
    delete temp_d;
  }

  if(object)
    delete object;
  
}

//change the size of the global transmission function 
//this function is used for dynamically determining the size
//of the glocal transmission function array.
void PhaseDiverseCDI::reallocate_object_memory(int new_nx,int new_ny){

  if(object)
    delete object;

  object = new Complex_2D(new_nx,new_ny);
  
  nx = new_nx;
  ny = new_ny;

}

//return the global transmision function.
Complex_2D * PhaseDiverseCDI::get_transmission(){
  return object;
}


//set the global transmission function (note that
//local transmission function will not be automatically updated).
//this function can be used to set the initial estimate.
void PhaseDiverseCDI::set_transmission(Complex_2D & new_transmission){

  int new_nx = new_transmission.get_size_x();
  int new_ny = new_transmission.get_size_y();

  if(!object)
    reallocate_object_memory(new_nx,new_ny);

  if(new_nx!=nx||new_ny!=ny){
    cout << "Can not set the transmission function in "
	 << "PhaseDiverseCDI because the complex array "
	 << "given does not have the same dimensions as "
	 << "the current transmission function"<<endl;
    return;
  }

  object->copy(new_transmission);
  
}

void PhaseDiverseCDI::add_new_position(BaseCDI * local,
				       double x, 
				       double y,
				       double alpha){
  
  weights_set = false;


  int lnx = local->get_size_x() ;
  int lny = local->get_size_y() ;
  
  singleCDI.push_back(local);
  x_position.push_back(x);
  y_position.push_back(y);
  this->alpha.push_back(alpha);
  single_result.push_back(new Complex_2D(lnx,lny));
  weights.push_back(new Double_2D(lnx,lny));

  cout << "Added position "<<singleCDI.size()-1<<endl;
  
  if(!object){
    x_min = -x;
    y_min = -y;
    reallocate_object_memory(lnx*scale,lny*scale);
  }

  //dynamically increase the object size to fit this frame in

  //what is the position (globally) of the first 
  //pixel if we don't increase the frame size  
  int extra_x=0;
  int extra_y=0;

  int global_x_min = get_global_x_pos(0,x)*scale;  //-x-x_min;
  int global_y_min = get_global_y_pos(0,y)*scale;  //-y-y_min;
  int global_x_max = get_global_x_pos(lnx,x)*scale; //lnx-x-x_min;
  int global_y_max = get_global_y_pos(lny,y)*scale; //lny-y-y_min;

  if(global_x_min<0){
    x_min = -x;
    extra_x += -global_x_min;
  }
  if(global_y_min<0){
    y_min = -y;
    extra_y += -global_y_min;
  }
  if(global_x_max>nx)
    extra_x += global_x_max-nx;
  if(global_y_max>ny)
    extra_y += global_y_max-ny;

  if(extra_x || extra_y)
    reallocate_object_memory(nx+extra_x, ny+extra_y);

}




void PhaseDiverseCDI::initialise_estimate(){

  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      object->set_real(i,j,1);
      object->set_imag(i,j,0);
    }
  }

  for(int i=0; i<singleCDI.size(); i++){
    add_to_object(i);
  }


}


void PhaseDiverseCDI::iterate(){

  static int total_iterations = 0;

  cout << "Iteration "<<total_iterations<<endl;

  for(int i=0; i<singleCDI.size(); i++){

    int x, y;

    //if(i!=0 && (total_iterations==2 || total_iterations==4))
    //  check_position(i);
    update_from_object(i);

    for(int j=0; j<iterations_per_cycle; j++){

      //temp
      /**      static bool flag = false;
      if(!flag){
	char buf[50];
	Double_2D result(nx,ny);
	Complex_2D trans(nx,ny);
	sprintf(buf,"itr_%i_%i_1.tiff",i,j);
	//      ((FresnelCDI*)singleCDI.at(i))->get_transmission_function(trans);
	trans.copy(*get_transmission());
	trans.get_2d(MAG,result);
	write_image(buf,result,false,0,1.0);
	flag=true;
	}**/
      
      singleCDI.at(i)->iterate();
      cout << "This Error="<<singleCDI.at(i)->get_error()<<endl;

    }

    if(!parallel)
      add_to_object(i);

  }
  
  if(parallel){

    scale_object(1-beta); 

    for(int i=0; i<singleCDI.size(); i++){
      add_to_object(i);
    }
  }

  total_iterations++;

};


void PhaseDiverseCDI::scale_object(double factor){

  set_up_weights();

  int frames = singleCDI.size();

  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){

      bool not_in_image = true;
      int n_probe = 0;

      while(not_in_image&&n_probe<frames){

	int i_ = (i/scale) + x_position.at(n_probe) + x_min; //get_local_x_pos(i,x_position.at(n_probe));  
	int j_ = (j/scale) + y_position.at(n_probe) + y_min; //get_local_y_pos(j,y_position.at(n_probe));
	
	if(i_>=0&&j_>=0&&i_<weights.at(n_probe)->get_size_x()
	   &&j_<weights.at(n_probe)->get_size_y())
	  not_in_image=(weights.at(n_probe)->get(i_,j_)==0); 
	
	n_probe++;

      }

      if(not_in_image){
	object->set_real(i,j,1.0);
	object->set_imag(i,j,0.0);
      }
      else{
	object->set_real(i,j,factor*object->get_real(i,j));
	object->set_imag(i,j,factor*object->get_imag(i,j));
      }	
      
    }
  }
  
}


void PhaseDiverseCDI::get_result(BaseCDI * local, Complex_2D & result){
  if(typeid(*local)==typeid(FresnelCDI)){
    ((FresnelCDI*)local)->get_transmission_function(result);
  }
  else
    result.copy(local->get_exit_surface_wave());
};

void PhaseDiverseCDI::set_result(BaseCDI * local, Complex_2D & result){
  
  if(typeid(*local)==typeid(FresnelCDI)){
    ((FresnelCDI*)local)->set_transmission_function(result);
  }
  else
    local->set_exit_surface_wave(result);
  
}

//frames need to be in order for this to work.
void PhaseDiverseCDI::adjust_positions(double step_size, bool forward){

  //make a copy of the tranmission function
  Complex_2D * pointer_to_object = object;
  int nx_c = nx;
  int ny_c = ny;
  double x_min_c = x_min; 
  double y_min_c = y_min; 
  bool parallel_c = parallel;  
  double beta_c = beta;

  //reset the parameters to make it work in parrallel with weights of one.
  parallel = false;

  //leave the first frame, as everything will be aligned to it.
  int n=1;
  int limit=singleCDI.size();
  if(!forward){
    n=singleCDI.size()-2;
    limit=-1;
  }


  Complex_2D temp_object(nx,ny);
  object = &temp_object;
  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      temp_object.set_real(i,j,1);
      temp_object.set_imag(i,j,0);
    }
  }

  beta = 1.0;
  weights_set=false;

  if(forward)
    add_to_object(0);
  else
    add_to_object(singleCDI.size()-1);
  
  //  beta = 0.5;
  //weights_set=false;


  while(n!=limit){
    //for(int n=1; n<singleCDI.size(); n++){

    //store some before information for later use
    double before_x = x_position.at(n);
    double before_y = y_position.at(n);

    //replace the transmission function with
    //just the 1st frame information

    int lnx = single_result.at(0)->get_size_x();
    int lny = single_result.at(0)->get_size_y();

    Complex_2D temp(lnx,lny);
    get_object_sub_grid(temp,before_x,before_y);

    object = &temp;
    nx = lnx;
    ny = lny;
    x_min = -x_position.at(n);
    y_min = -y_position.at(n);
    
    /**    for(int i=0; i<lnx; i++){
      for(int j=0; j<lny; j++){
	temp.set_real(i,j,1);
	temp.set_imag(i,j,0);
      }
      }**/

    /**    char name[80];
    Double_2D pic(1024,1024);
    object->get_2d(MAG,pic);
    write_image("pic_0.tiff",pic,false,0,1);**/

    //if frame m overlaps with frame n
    //add it to the current array

    /**    if(forward){
      add_to_object(0);
      beta=0.5;
      weights_set=false;
      for(int m=1; m < n; m++){
	add_to_object(m);
      }
    }
    else{
      add_to_object(singleCDI.size()-1);
      beta=0.5;
      weights_set=false;
      for(int m=singleCDI.size()-2; m > n; m--){
	add_to_object(m);
      }
      }**/

    //    char name[80];
    //Double_2D pic(1024,1024);
    /**    object->get_2d(MAG,pic);
    sprintf(name,"pic_%i.tiff",n);
    write_image(name,pic,false, 0,1);**/

    //do an iteration of the frame on it's own.
    singleCDI.at(n)->iterate();
    get_result(singleCDI.at(n),*(single_result.at(n)));

    Double_2D temp_single(lnx,lny);
    Double_2D temp_others(lnx,lny);
    for(int i=0; i<lnx; i++){
      for(int j=0; j<lny; j++){
	temp_single.set(i,j,
			1-single_result.at(n)->get_mag(i,j));
	temp_others.set(i,j,
			1-object->get_mag(i,j));
      }
    }

    /**    static int counter = 0;
    char buff[90];
    sprintf(buff,"temp_single_%i.tiff",counter);
    write_image(buff,temp_single);
    sprintf(buff,"temp_object_%i.tiff",counter);
    write_image(buff,temp_others);
    counter++;**/

    int new_x, new_y;

    /** align(temp_single, temp_others, 
	  new_x, new_y,
	  8,
	  -100,+100,
	  -100,+100); **/
    
    align_even_better(temp_single, temp_others, 
		      new_x, new_y,
		      -20, +20,
		      -20, +20,
		      &temp_single,
		      //weights.at(n),
		      &temp_others,
		      0.2); 

    x_position.at(n) = before_x+new_x;
    y_position.at(n) = before_y+new_y;

    cout << "cross-correlation: moving prob " << n << " from ("
	 << before_x << "," << before_y << ") "
	 << "-> (" << x_position.at(n) 
	 << "," <<  y_position.at(n) << ")."
	 << endl;

    /**    check_position(n,8,0);
    
    cout << "error: moving prob " << n << " from ("
	 << before_x << "," << before_y << ") "
	 << "-> (" << x_position.at(n) 
	 << "," <<  y_position.at(n) << ")."
	 << endl;**/

    //if we moved the positions, we need to recalculate the 
    //weights for that frame.

    /**    double x_difference = x_position.at(n) - before_x;
	   double y_difference = y_position.at(n) - before_y;
	   
	   if( x_difference!=0 || y_difference!=0 ){
      //      weights_set=false;
      
      //change all the consecutive positions accordingly:
      if(forward){
	for(int m=n+1; m < singleCDI.size(); m++){
	  x_position.at(m) +=  x_difference;
	  y_position.at(m) +=  y_difference;
	}
      }
      else{
	for(int m=n-1; m >= 0; m--){
	  x_position.at(m) +=  x_difference;
	  y_position.at(m) +=  y_difference;
	}
      }
      
      }**/

    object = &temp_object;
    nx = nx_c;
    ny = ny_c;
    x_min = x_min_c; 
    y_min = y_min_c;
    add_to_object(n);

    if(forward)
      n++;
    else
      n--;

  }
  
  //set all the parameters back to normal
  object = pointer_to_object;
  /**  nx = nx_c;
  ny = ny_c;
  x_min = x_min_c; 
  y_min = y_min_c; **/
  parallel = parallel_c;
  beta=beta_c;

  //force the weights to be recalculated
  weights_set=false;


  //reset the object global transmission function
  scale_object(1-beta); 
  
  for(int i=0; i<singleCDI.size(); i++){
    add_to_object(i);
  }


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
    cout << "Giving up on probe "<< n_probe << ". Could not find " 
	 << "it's position. Returning to the original. " 
	 << endl;
      return FAILURE;
  }
  
  BaseCDI * single = singleCDI.at(n_probe);
  
  double best_x=x;
  double best_y=y;
  double best_error=100;
   
  double new_x;
  double new_y;
  double size_x = single_result.at(n_probe)->get_size_x();
  double size_y = single_result.at(n_probe)->get_size_y();
  
  //try the 9 positions around the current one
  //record the one with the lowest error metric.

  //copy some local stuff.
  Complex_2D * object_copy = new Complex_2D(object->get_size_x(),object->get_size_y());
  object_copy->copy(*object);
  Complex_2D * single_copy = new Complex_2D(size_x,size_y);
  single_copy->copy(*single_result.at(n_probe));
  double beta_c = beta;

  beta = 0.5;
  weights_set=false;

  for(int i=-1; i<2; i++){
    for(int j=-1; j<2; j++){

      new_x=x+i*step_size;
      new_y=y+j*step_size;

      //checking whether the corrds are still within the 
      //boundary of the image
      //      if(new_x>0 && new_x<size_x && new_y>0 && new_y<size_y){

      x_position.at(n_probe)=new_x;
      y_position.at(n_probe)=new_y;

      add_to_object(n_probe);
      update_from_object(n_probe);

      /**      if(n_probe==4&&new_x==3294&&new_y==2917){
	static int h = 0;
	char blah[80];
	Double_2D temp(1024,1024);
	single_result.at(n_probe)->get_2d(PHASE,temp);
	sprintf(blah,"single_prob%i_try%i.tiff",n_probe,h);
	write_image(blah,temp,true,-0.1,0);

	object->get_2d(PHASE,temp);
	sprintf(blah,"object_prob%i_try%i.tiff",n_probe,h);
	write_image(blah,temp,true,-0.1,0);

	h++;

	cout << "error at ("<<new_x<<","<<new_y
	     << ")"<<" is "<<single->get_error()<<endl;

	     }**/

      single->iterate();
	
      if(single->get_error()<best_error){
	best_error = single->get_error();
	best_x = new_x;
	best_y = new_y;
      }

      //reset the object estimate:
      object->copy(*object_copy);
      //and
      set_result(single,*single_copy);
      single_result.at(n_probe)->copy(*single_copy);

    }
  }

  delete object_copy;
  delete single_copy;

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
  

  beta = 1.0;
  weights_set=false;

  return SUCCESS;
  
}


void PhaseDiverseCDI::set_up_weights(){
  
  if(weights_set)
    return;
  else
    weights_set=true;

  for(int n_probe=0; n_probe<singleCDI.size(); n_probe++){
    
    int lnx = single_result.at(n_probe)->get_size_x();
    int lny = single_result.at(n_probe)->get_size_y();
    BaseCDI * this_CDI = singleCDI.at(n_probe);
    double value;
    Double_2D illum_mag(lnx,lny);
    double max;

    if(typeid(*this_CDI)==typeid(FresnelCDI)){
      const Complex_2D & illum = ((FresnelCDI*)(this_CDI))->get_illumination_at_sample();
      illum.get_2d(MAG,illum_mag);
      max = illum_mag.get_max();
    }
    
    //    Double_2D support(lnx,lny);
    //this_CDI->get_support(support);
    
    for(int i_=0; i_< lnx; i_++){
      for(int j_=0; j_< lny; j_++){
	
	if(typeid(*this_CDI)==typeid(FresnelCDI))
	  value = illum_mag.get(i_,j_)/max;
	else
	  value = 1;

	value=beta*alpha.at(n_probe)*pow(value,gamma);

	if(this_CDI->get_support().get(i_,j_)<=0.0)
	  value = 0;
	//else
	//  cout << support.get(i_,j_) <<endl;

	weights.at(n_probe)->set(i_,j_,value);
      }
    } 
  }
  
  if(parallel){
    
    Double_2D weight_norm(object->get_size_x(),object->get_size_y());
    
    for(int n=0; n<singleCDI.size(); n++){
      
      double x = x_position.at(n);
      double y = y_position.at(n);
      
      int this_size_x = single_result.at(n)->get_size_x();
      int this_size_y = single_result.at(n)->get_size_y();
      
      
      for(int i_=0; i_< this_size_x; i_++){
	for(int j_=0; j_< this_size_y; j_++){
	  
	  int i = get_global_x_pos(i_,x); // (i_-x-x_min)*scale;
	  int j = get_global_y_pos(j_,y); // (j_-y-y_min)*scale;
	  
	  if(i>=0&&j>=0&&i<nx&&j<ny){
	    
	    double new_weight = weights.at(n)->get(i_,j_);
	    double old_weight = weight_norm.get(i,j);
	    
	    if(n==0)
	      weight_norm.set(i,j,new_weight);
	    else
	      weight_norm.set(i,j,new_weight+old_weight);
	  }
	}
      }
    }
    
    for(int n=0; n<singleCDI.size(); n++){
      
      double x = x_position.at(n);
      double y = y_position.at(n);
      
      int this_size_x = single_result.at(n)->get_size_x();
      int this_size_y = single_result.at(n)->get_size_y();
      
      for(int i_=0; i_< this_size_x; i_++){
	for(int j_=0; j_< this_size_y; j_++){
	  
	  int i = get_global_x_pos(i_,x); 
	  int j = get_global_y_pos(j_,y); 

	  //	  int i = (i_-x-x_min)*scale;
	  // int j = (j_-y-y_min)*scale;
	  
	  if(i>=0&&j>=0&&i<nx&&j<ny){
	    
	    double old_weight = weights.at(n)->get(i_,j_);
	    double norm = weight_norm.get(i,j);

	    if(norm<=0)
	      weights.at(n)->set(i_,j_,0);
	    else
	      weights.at(n)->set(i_,j_,old_weight/norm);

	  }
	}
      }
    } 
    //    write_image("weight_norm.tiff",weight_norm);
 
  }

  //  write_image("weight_0.tiff",*(weights.at(0)));
 };


  


void PhaseDiverseCDI::add_to_object(int n_probe){
  
  get_result(singleCDI.at(n_probe),*(single_result.at(n_probe)));  
  set_up_weights();
  
  double x_offset = x_position.at(n_probe);
  double y_offset = y_position.at(n_probe);
  
  Complex_2D * small = 0;
  if(scale!=1){
    small = new Complex_2D(single_result.at(n_probe)->get_size_x()*scale,
			   single_result.at(n_probe)->get_size_y()*scale);
    interpolate( *(single_result.at(n_probe)), *small);  
  }
  else
    small = single_result.at(n_probe);
	      
  Double_2D & this_weight = *(weights.at(n_probe));
  
  int lnx = small->get_size_x();
  int lny = small->get_size_y();
  
  //round off to a first approx. Will be fixed later.
  int i, j; //local coors (i,j) are global (object) coors

  //using the local coordinate system
  for(int i_=0; i_< lnx; i_++){
    for(int j_=0; j_< lny; j_++){
      
      double weight = this_weight.get(i_/scale,j_/scale);
      
      if(weight!=0){

	i = i_ - x_offset*scale -x_min*scale; //get_global_x_pos(i_,x_offset); //(i_-x_offset-x_min)*scale;
	j = j_ - y_offset*scale -y_min*scale; //get_global_y_pos(j_,y_offset); //(j_-y_offset-y_min)*scale;
	
	//	i = (i_-x_offset-x_min)*scale;
	//	j = (j_-y_offset-y_min)*scale;
      
	if(i>=0&&j>=0&&i<nx&&j<ny){

	  double new_real = weight*small->get_real(i_,j_)
	    +object->get_real(i,j);
	  
	  double new_imag = weight*small->get_imag(i_,j_)
	    +object->get_imag(i,j);
	  
	  if(!parallel){
	    new_real -= weight*object->get_real(i,j);
	    new_imag -= weight*object->get_imag(i,j);
	  }
	  
	  object->set_real(i,j,new_real);
	  object->set_imag(i,j,new_imag);
	  
	  //update_to_object_sub_grid(i,j,new_real,new_imag);
	}
      }
      
    }
  } 

  if(scale!=1)
    delete small;
  
}



/** void PhaseDiverseCDI::update_to_object_sub_grid(int i, int j, 
						double real_value, 
						double imag_value,
						

){
  
  for(int di=0; di < scale; di++){
    for(int dj=0; dj < scale; dj++){
      object->set_real(i+di,j+dj,real_value);
      object->set_imag(i+di,j+dj,imag_value);
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
      real_value+=object->get_real(i+di,j+dj);
      imag_value+=object->get_imag(i+di,j+dj);
    }
  }

  real_value=real_value/((double)scale*scale);
  imag_value=imag_value/((double)scale*scale);

  }**/


void PhaseDiverseCDI::get_object_sub_grid(Complex_2D & result,
					  double x_offset,
					  double y_offset){
  
  

  //using the local coordinate system
  for(int i_=0; i_< result.get_size_x(); i_++){
    for(int j_=0; j_< result.get_size_y(); j_++){

      //local coors (i,j) are global (object) coors  
      int i = get_global_x_pos(i_,x_offset);
      int j = get_global_y_pos(j_,y_offset);
      
      if(i>=0 && j>=0 && i<nx && j<ny){

	double new_real = object->get_real(i,j); 
	double new_imag = object->get_imag(i,j); 
	
	/**	update_from_object_sub_grid(i, j, 
				    new_real, 
				    new_imag);**/

	result.set_real(i_,j_,new_real);
	result.set_imag(i_,j_,new_imag);

      }
    }
  }
  
};





void PhaseDiverseCDI::update_from_object(int n_probe){

  set_up_weights();

  double x_offset = x_position.at(n_probe);
  double y_offset = y_position.at(n_probe);

  Complex_2D * small = single_result.at(n_probe);
  Complex_2D * scaled_object = 0;

  if(scale!=1){
    scaled_object = new Complex_2D(object->get_size_x()/scale,
				   object->get_size_y()/scale);
    shrink(*object, *scaled_object);  
  }
  else
    scaled_object = object;
  
  Double_2D & this_weight = *weights.at(n_probe);

  int lnx = small->get_size_x();
  int lny = small->get_size_y();
  
  //round off to a first approx. Will be fixed later.
  int i, j; //local coors (i,j) are global (object) coors

  //using the local coordinate system
  for(int i_=0; i_< lnx; i_++){
      for(int j_=0; j_< lny; j_++){
	
	if(this_weight.get(i_,j_)!=0){

	  i = get_global_x_pos(i_,x_offset); //(i_-x_offset-x_min)*scale;
	  j = get_global_y_pos(j_,y_offset); //(j_-y_offset-y_min)*scale;

	  //	  i = (i_-x_offset-x_min)*scale;
	  //  j = (j_-y_offset-y_min)*scale;
      
	  if(i>=0 && j>=0 && i<nx && j<ny){
	    
	    double new_real = scaled_object->get_real(i,j);
	    double new_imag = scaled_object->get_imag(i,j);
	    
	    //update_from_object_sub_grid(i,j,new_real,new_imag);
	  
	    small->set_real(i_,j_,new_real);
	    small->set_imag(i_,j_,new_imag);
	  }
	}
      }
  }	    


  set_result(singleCDI.at(n_probe),*(single_result.at(n_probe)));
  if(scale!=1)
    delete scaled_object;
  
};
 
