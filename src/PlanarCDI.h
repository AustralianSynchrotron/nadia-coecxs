#ifndef PLANARCDI_H
#define PLANARCDI_H

//#include <map>
//#include <string>
#include "PhaseRetrievalBase.h"

//forward declarations
class Complex_2D;

/**
 * @file PlanarCDI.h
 * @class PlanarCDI
 * @author  Nadia Davidson <> 
 *
 * @brief  The class which performs planar CDI reconstuction.
 *
 * The fundamental class used for performing planar CDI reconstuction.
 * A projection object should be created with Complex_2D object which
 * you wish the projection to act upon. The Complex_2D should already
 * be filled with an inital guess, or the
 * PlanarCDI::initialise_estimate() method can be used to set up an
 * inital guess with random numbers inside the support and zero
 * outside.
 *
 * The support and intensity in the detector plane must be set using
 * the approproate setter methods, and the algorithm and relaxation
 * parameters may also be set if something other than the default it
 * required. The PlanarCDI::set_custom_algorithm allows the
 * recostruction algorithm to be customised by the user. See page 22
 * of the H.M. Quiney review: TUTORIAL REVIEW, Coherent diffractive
 * imaging using short wavelength light sources, Journal of Modern
 * Optics, 2010, DOI: 10.1080/09500340.2010.495459 for the formalism.
 *
 * The reconstruction can be the performed with the
 * PlanarCDI::iterate() method.  This is automatically update the
 * values in the Complex_2D estimate and the current error. The error
 * can be retrieved using PlanarCDI::get_current_error(). The
 * algorithm, relaxation parameter, support and intensity can be reset
 * at any point during reconstuction.
 */
class PlanarCDI : public PhaseRetrievalBase{

 public:
  
  PlanarCDI(Complex_2D * complex);
  //  ~PlanarCDI();
  
  int iterate();
  
};

#endif
