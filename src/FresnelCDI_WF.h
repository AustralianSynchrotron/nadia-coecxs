#ifndef FCDI_WF_H
#define FCDI_WF_H

#include "PlanarCDI.h"

//forward declarations
class Complex_2D;

/**
 * @file FresnelCDI_WF.h
 * @class FresnelCDI_WF
 * @author  Nadia Davidson <> 
 *
 * @brief  The class which performs Fresnel CDI white 
 * field reconstuction.
 *
 *
 */

class FresnelCDI_WF: public PlanarCDI{

 private:

  double wavelength;
  double zone_to_focal_length;
  double focal_to_detector_length;
  double pixel_length;
  Complex_2D forward_coefficient;
  Complex_2D backward_coefficient;
  
 public:
  
  FresnelCDI_WF(Complex_2D & initial_guess,
		double beam_wavelength,
		double zone_focal_length,
		double focal_detector_length,
		double pixel_size);
  
  int iterate();  
  void initialise_estimate(int seed=0);
  void set_algorithm(int alg){
    std::cout << "WARNING: Algorithm can not be set when performing Frenel "
	      << "CDI white field recovery" << std::endl;
  }
  void get_algorithm(){
    std::cout << "Using the default Frenel CDI "
	      << "white field recovery algorithm" << std::endl;
  }

};

#endif
