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
 * @brief  The class which performs planar CDI reconstuction.
 *
 * The fundamental class used for performing planar CDI reconstuction.
 * A projection object should be created with Complex_2D object which
 * you wish the projection to act upon. The Complex_2D should already
 * be filled with an inital guess, or the
 * FresnelCDI_WF::initialise_estimate() method can be used to set up an
 * inital guess with random numbers inside the support and zero
 * outside.
 *
 * The support and intensity in the detector plane must be set using
 * the approproate setter methods, and the algorithm and relaxation
 * parameters may also be set if something other than the default it
 * required. The FresnelCDI_WF::set_custom_algorithm allows the
 * recostruction algorithm to be customised by the user. See page 22
 * of the H.M. Quiney review: TUTORIAL REVIEW, Coherent diffractive
 * imaging using short wavelength light sources, Journal of Modern
 * Optics, 2010, DOI: 10.1080/09500340.2010.495459 for the formalism.
 *
 * The reconstruction can be the performed with the
 * FresnelCDI_WF::iterate() method.  This is automatically update the
 * values in the Complex_2D estimate and the current error. The error
 * can be retrieved using FresnelCDI_WF::get_current_error(). The
 * algorithm, relaxation parameter, support and intensity can be reset
 * at any point during reconstuction.
 */

class FresnelCDI_WF: public PlanarCDI{

 private:

  double zone_to_focal_length;

  double wavelength;

  double focal_to_detector_length;

  double pixel_length;

  
  Complex_2D * forward_coefficients_const;
  Complex_2D * backward_coefficients_const;
  
 public:
  
  FresnelCDI_WF(Complex_2D & initial_guess,
		double beam_wavelength,
		double zone_focal_length,
		double focal_detector_length,
		double pixel_size);
  ~FresnelCDI_WF();
  
  int iterate();
  
  void initialise_estimate(int seed=0);
  
};

#endif
