#ifndef FCDI_H
#define FCDI_H

#include "PlanarCDI.h"

//forward declarations
class Complex_2D;

/**
 * @file FresnelCDI.h
 * @class FresnelCDI
 * @author  Nadia Davidson <> 
 *
 * @brief  The class which performs planar CDI reconstuction.
 *
 * The fundamental class used for performing planar CDI reconstuction.
 * A projection object should be created with Complex_2D object which
 * you wish the projection to act upon. The Complex_2D should already
 * be filled with an inital guess, or the
 * FresnelCDI::initialise_estimate() method can be used to set up an
 * inital guess with random numbers inside the support and zero
 * outside.
 *
 * The support and intensity in the detector plane must be set using
 * the approproate setter methods, and the algorithm and relaxation
 * parameters may also be set if something other than the default it
 * required. The FresnelCDI::set_custom_algorithm allows the
 * recostruction algorithm to be customised by the user. See page 22
 * of the H.M. Quiney review: TUTORIAL REVIEW, Coherent diffractive
 * imaging using short wavelength light sources, Journal of Modern
 * Optics, 2010, DOI: 10.1080/09500340.2010.495459 for the formalism.
 *
 * The reconstruction can be the performed with the
 * FresnelCDI::iterate() method.  This is automatically update the
 * values in the Complex_2D estimate and the current error. The error
 * can be retrieved using FresnelCDI::get_current_error(). The
 * algorithm, relaxation parameter, support and intensity can be reset
 * at any point during reconstuction.
 */

class FresnelCDI: public PlanarCDI{

 protected:

  Complex_2D & illumination;
  double wavelength;
  double sample_to_detector_length;
  double pixel_length;
  double norm;
  Complex_2D B_s;
  Complex_2D B_d;
  
 public:
  
  FresnelCDI(Complex_2D & initial_guess,
	     Complex_2D & white_field,
	     double beam_wavelength,
	     double focal_detector_length,
	     double focal_sample_length,
	     double pixel_size,
	     double normalisation=1.0
	     );

  virtual void project_intensity(Complex_2D & c); 
  virtual void initialise_estimate(int seed=0); //FIXME: seed doesn't so anything
  virtual void get_transmission_function(Complex_2D & result);
  
 protected:
  virtual void propagate_to_sample(Complex_2D & c);
  virtual void propagate_to_detector(Complex_2D & c);

};

#endif
