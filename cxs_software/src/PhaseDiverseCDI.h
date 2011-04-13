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
  std::vector<Double_2D *> beam_shape;
  std::vector<double> x_position;
  std::vector<double> y_position;
  
  //parameters controlling the feedback
  double beta;
  double gamma;
  std::vector<double> alpha;
  
  Complex_2D & object;
  Double_2D weight_norm;

  int iterations_per_cycle;
  int scale;

 public:

  PhaseDiverseCDI(Complex_2D & object, int granularity=1,
		  double beta=1.0, double gamma=1.0);
  ~PhaseDiverseCDI();

  void iterate();
  
  void add_new_position(PlanarCDI * local, 
			Double_2D & beam_shape,
			double x=0, double y=0, 
			double alpha=1);
 
  void set_iterations_per_cycle(int iterations){
    iterations_per_cycle = iterations;
  }

  /**  void set_feedback_parameter(double beta){
    this->beta = beta;
  };

  void set_amplification_factor(double gamma){
    this->gamma = gamma;
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
    
    };**/

  void initialise_estimate();

 private:

  void add_to_object(int n_probe, double weight, double old_weight);
  int update_from_object_with_shift(int n_probe, double shift=4, int tries=0);
  void update_from_object(int n_probe);
  void get_result(PlanarCDI * local, Complex_2D & result);
  void set_result(PlanarCDI * local, Complex_2D & result);
  void update_to_object_sub_grid(int i, int j, 
				 double real_value, 
				 double imag_value);
  void update_from_object_sub_grid(int i, int j, 
				 double & real_value, 
				 double & imag_value);

  void get_weights(int n_probe, Double_2D & weights);


};


#endif
