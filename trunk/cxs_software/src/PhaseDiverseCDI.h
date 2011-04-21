/**
 * @file FresnelCDI.h
 * @class FresnelCDI
 * @author  Nadia Davidson <nadiamd@unimelb.edu.au> 
 *
 * @brief  The class which performs Fresnel CDI reconstruction.
 *
 * The class used for performing Fresnel CDI reconstruction (for
 * white-field reconstruction see FresnelCDI_WF). It inherits most
 * methods from PlanarCDI, so please look at the documentation of this
 * class also. Although there are some differences in the underlying
 * code between this class and the planar case, the interface is
 * generally unchanged. Therefore users should refer to the
 * instructions for PlanarCDI to understand how to use a FresnelCDI
 * object in their own code. Only the differences relevant to users
 * will be documented here.
 *
 */

#ifndef PHASED_H
#define PHASED_H

#include "PlanarCDI.h"
#include <vector>

//forward declarations
class Complex_2D;
class PhaseDiverseCDI{

 protected:

  std::vector<PlanarCDI*> singleCDI;
  std::vector<Complex_2D*> single_result;
  //std::vector<Double_2D *> beam_shape;
  std::vector<Double_2D *> weights;
  std::vector<double> x_position;
  std::vector<double> y_position;
  
  //parameters controlling the feedback
  double beta;
  double gamma;
  std::vector<double> alpha;

  Complex_2D * object;

  int iterations_per_cycle;
  int scale;
  int nx,ny;
  
  int x_min;
  int y_min;
  bool parallel; 
  bool weights_set;

 public:

  PhaseDiverseCDI(double beta=1.0, 
		  double gamma=1.0,
		  bool parallel=false,
		  int nx=0,
		  int ny=0,
		  int min_x=0,
		  int min_y=0,
		  int granularity=1);
  
  ~PhaseDiverseCDI();

  void iterate();
  
  void add_new_position(PlanarCDI * local, 
			double x=0, double y=0, 
			double alpha=1);
 
  void set_iterations_per_cycle(int iterations){
    iterations_per_cycle = iterations;
  }

  void set_feedback_parameter(double beta){
    this->beta = beta;    
    weights_set = false;

  };

  void set_amplification_factor(double gamma){
    this->gamma = gamma;
    weights_set = false;
 
  };

  void set_probe_scaling(int n_probe, double alpha){
    if(this->alpha.size() > n_probe)
      this->alpha.at(n_probe)=alpha;
    else{
      std::cout << "In PhaseDiverseCDI::set_probe_scaling, "
	   << "the probe number given is too large. "
	   << "Are you trying to set the value of alpha "
		<< "before calling add_new_position?" <<std::endl;
    }
    weights_set = false;

  };

  void initialise_estimate();
  Complex_2D * get_transmission();
  void set_transmission(Complex_2D & new_transmission);
  void adjust_positions(double step_size=4, bool forwards=true);
  double get_final_x_position(int n_probe){
    return x_position.at(n_probe);}
  double get_final_y_position(int n_probe){
    return y_position.at(n_probe);}



 private:

  void add_to_object(int n_probe);
  int check_position(int n_probe, double shift=4, int tries=0);
  void update_from_object(int n_probe);
  void get_result(PlanarCDI * local, Complex_2D & result);
  void set_result(PlanarCDI * local, Complex_2D & result);
  void update_to_object_sub_grid(int i, int j, 
				 double real_value, 
				 double imag_value);
  void update_from_object_sub_grid(int i, int j, 
				 double & real_value, 
				 double & imag_value);

  //void get_weights(int n_probe, Double_2D & weights);

  //void set_up_weight_norm();

  void reallocate_object_memory(int new_nx,int new_ny);

  void set_up_weights();
};


#endif
