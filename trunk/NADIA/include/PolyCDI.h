// Copyright 2012 T'Mir Julius for The ARC Centre of Excellence in 
// Coherent X-ray Science. This program is distributed under the GNU  
// General Public License. We also ask that you cite this software in 
// publications where you made use of it for any part of the data     
// analysis. 

/**
 * @file PolyCDI.h
 * @class PolyCDI
 * @author T'Mir Julius <cxs-softwarse@physics.unimelb.edu.au> 
 *
 * @brief  The class which performs partially coherent 
 * reconstruction.
 *
 * This class can be used to perform partially coherent reconstruction
 * of plane-wave partialy CDI data. The number of modes and legendre 
 * polynomials used is defined by the user. The number of Legendre polynomial * must exceed the number of modes.
 */

#ifndef POLYCDI_H
#define POLYCDI_H

#include <BaseCDI.h>
#include <vector>
#include <cstring>

#define WL 0
#define WEIGHT 1

//forward declarations
class Complex_2D;
class PolyCDI:public BaseCDI {

protected:


  /** A vector holding a pointer to each of the FresnelCDI/PlanarCDI objects */
  std::vector<Complex_2D> singleCDI;

  /** A vector holding the initial wave described as a series of modes.*/
  std::vector<Complex_2D> singlemode;

  /** A vector holding the weighting function for each frame */  
  std::vector<Double_2D * > weights; 

  /** A vector of the transverse positions in x */
  std::vector<double> x_position;

  /** A vector of the transverse positions in y */
  std::vector<double> y_position;

  /** parameters controlling the feedback */
  double beta;

  /** An array of the eigenvectors for the system of JC=nSC */
  std::vector<double> eigen;

  /** The number of iterations to perform for each 'local' frame */
  int iterations_per_cycle;

  /** The number of orthogonal components for the light source.*/
  int nleg;

  /** The number of component modes for the light source. nleg must 
    be bigger than nmode*/
  int nmode;

  /** Number of wavelengths in the spectrum */
  double nlambda;

  /** Padding required to expand the wavelengths */
  int paddingx;
  int paddingy;

  /** Some features of the spectrum */
  Double_2D  spectrum;

  /** The calculated intensity at the detector */
  Double_2D intensity_sqrt_calc;

  double lambdac;

  /** The matrices describing the source properties.*/
  Double_2D * magnitude;

  /** a flag for running in either series or parallel mode */
  bool parallel; 

public:

  PolyCDI(Complex_2D & initial_guess,
      double beta=1.0,
      int n_best=1,
      bool parallel=0
      );


  enum {CROSS_CORRELATION,MINIMUM_ERROR};

  /** 
   * Destructor for PhaseDiverseCDI
   */
  ~PolyCDI();


  /**
   * Initialise the estimate of the 'global' object function. The
   * initialisation is performed using the current estimate for each
   * 'local' frame. For this reason you must first initialise each
   * FresnelCDI/PlanarCDI object individually before calling this
   * function.
   */
  void initialise_estimate();

  void initialise_estimate(int seed);

  /** scale the highest occupancy mode 
   * this overwrites the function of the same
   * name in BaseCDI
   */
  void scale_intensity(Complex_2D & c);

  /**
    expand wavelengths from central wavelength 
   */
  void expand_wl(Complex_2D & c);


  /**
   * for scaling the transmission
   */
  //  void scale_intensity(Complex_2D & c);

  /**
   * add the intensities across all modes 
   */
  Double_2D sum_intensity(std::vector<Complex_2D> & c);

  /**
   * The iterate the algorithm. This overwrites the
   * the class of the same name in BaseCDI.
   */ 

  int iterate();


  /**
   * calculate the transmission function by dividing 
   * the highest occupacy mode at the source by the 
   * highest occupancy mode at the detector
   */
  void update_transmission();

  /////////////////////////////////
  // Get and setter methods
  /////////////////////////////////

  /**
   * Set the number of 'local' frame iterations before updating the
   * result to the 'global' object function.
   *
   * @param iterations The number of 'local' iterations.
   */
  void set_iterations_per_cycle(int iterations){
    iterations_per_cycle = iterations;
  }

  /**
   * This function allows you to access the 'global' sample function.
   *
   * @return The current estimate of either the transmission (for
   * FresnelCDI) or exit-surface-wave (for PlanarCDI)
   */
  Complex_2D  get_transmission();

  /**
   * Set the 'global' sample function. This method could be used, for
   * example, instead of 'initialise_estimate' to initialise the
   * results to those from a previous reconstruction.
   *
   * @param new_transmission The transmission function 
   * exit-surface-wave to copy.
   */
  void set_transmission(Complex_2D & new_transmission);

  /**
   * Set the polychromatic spectrum wavelengths and weights
   * from a Double_2D()
   *
   * @param spectrum_in A Double_2D() containing the wavelength
   * and weights 
   */
  void set_spectrum(Double_2D spectrum_in);

  /**
   * Set the polychromatic spectrum wavelengths and weights
   * from a file
   *
   * @param file_in a string containing the filename of a SPECTRA
   * test file
   */
  void set_spectrum(std::string file_name);


  /**
   * calculate and return the current intensity of the modes multiplied
   * by the transmissoion fnction
   */
  Double_2D get_intensity(){
    return(intensity_sqrt_calc);
  };

  /**
   * Propagates the modes to the detector. Specifically for use with 
   * the simulations
   */
  Double_2D propagate_modes_to_detector();

  /**
   * Returns a given mode. If the mode requested is too big, it returns
   * the final mode
   */
  Complex_2D get_mode(int mode);

  /**
   *Propagate the CDI's 
   */
  void propagate_from_detector(Complex_2D & c);
  void propagate_to_detector(Complex_2D & c);

  /**
   * Set the minimum contribution of a mode as a proportion of the 
   * dominant mode for it to be included in the reconstruction 
   */
private:

  /**
   * This function is used to reallocate memory for the 'global'
   * sample object. It is only used when add_new_position is called,
   * in the case that the frame does not fit within the bounds of the
   * current object array.
   */
  void reallocate_transmission_memory(int new_nx,int new_ny);

  /**
   * This function is used during reconstruction in parallel mode. It
   * is similar to simply scaling all the elements of the 'global'
   * sample array, (e.g. with object.scale(factor). However, this
   * method preserves the value of the elements outside the boundary
   * of the sample. ie. it checks the weights of all the local frames
   * to work out which area to scale.
   * 
   * @param factor The scaled elements go from 'a' to 'factor*a'.
   *
   */
  void scale_object(double factor);

  /** These functions pad and unpad the support in order to allow for
   * larger supports when scaling the intensity of low energies.
   */

  void pad_support();
  void unpad_support();

};


#endif
