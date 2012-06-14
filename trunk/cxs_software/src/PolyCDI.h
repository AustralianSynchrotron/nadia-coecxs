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

#ifndef PCDI_H
#define PCDI_H

#include "BaseCDI.h"
#include <vector>

//forward declarations
class Complex_2D;
class PolyCDI:public BaseCDI {

protected:


  /* A vector holding a pointer to each of the FresnelCDI/PlanarCDI objects */
  std::vector<Complex_2D> singleCDI;

  /*A vector holding the initial wave described as a series of modes.*/
  std::vector<Complex_2D> singlemode;

  /** A vector holding the weighting function for each frame */  
  std::vector<Double_2D * > weights; 

  /* A vector of the transverse positions in x */
  std::vector<double> x_position;

  /* A vector of the transverse positions in y */
  std::vector<double> y_position;

  //parameters controlling the feedback
  double beta;

  /** The current estimate of the transmission function */
  Complex_2D transmission;

  /** An array of the eigenvectors for the system of JC=nSC */
  std::vector<double> eigen;

  /** The number of iterations to perform for each 'local' frame */
  int iterations_per_cycle;

  /** The number of orthogonal components for the light source.*/
  int nleg;

  /** The number of component modes for the light source. nleg must 
    be bigger than nmode*/
  int nmode;

  /** Coherence lengths */
  double lcy;
  double lcx;

  /** Length of area being imaged **/
  double lx;
  double ly;

  /** The size of a pixel */
  double pxsize;
  double pysize;

  /** The minimum value of the contribution of a mode 
    * as a proportion of the dominant mode
    */
  double threshold;

  /** The matrices describing the source properties.*/
  Complex_2D * jmatrix;
  Complex_2D * hmatrix;
  Complex_2D * smatrix;
  
  Double_2D * magnitude;

  /** a flag for running in either series or parallel mode */
  bool parallel; 

public:

  PolyCDI(Complex_2D & initial_guess,
      double beta=1.0,
      double lcx=0,
      double lcy=0,
      double lx=0,
      double ly=0,
      int n_best=1.0,
      bool parallel=0
      );


  enum {CROSS_CORRELATION,MINIMUM_ERROR};

  /** 
   * Construct a PhaseDiverseCDI object for phase diverse or
   * ptychographic reconstruction. The data should be entered later
   * using the 'add_new_position' function.
   *
   * @param beta The feedback parameter. By default this is 1 (no feedback).
   * @param gamma The amplification factor. By default this is 1 (no
   * amplification).
   * @param parallel true - run in parallel mode, false - run in series
   *        mode. By default series mode is set.
   * @param granularity  A factor which controls sub-pixel alignment. 
   *   1 = regular, 2 = 2 'global' pixels for every 1 'local' pixel.
   *   i.e. 2 allows alignment to within half a pixel. 
   *   NOTE: This is not currently working properly!
   */

  /*  PolyCDI(double beta=1.0, 
      double gamma=1.0,
      bool parallel=false,
      int granularity=1
      );
   */
  /** 
   * Destructor for PhaseDiverseCDI
   */
  ~PolyCDI();


  /**
   *  
   *
   */
  void add_new_position(BaseCDI * local, 
      double x=0, double y=0, 
      double alpha=1);

  /**
   * Initialise the estimate of the 'global' object function. The
   * initialisation is performed using the current estimate for each
   * 'local' frame. For this reason you must first initialise each
   * FresnelCDI/PlanarCDI object individually before calling this
   * function.
   */
  void initialise_estimate();

  void initialise_estimate(int seed);

  /** Initialise the wave matrices */
  void initialise_matrices(int leg, int modes);

  /**
   *uses the complex_2d multiply function to apply 
   *the transmission function
   */
  void apply_transmission(Complex_2D & c);

  /* scale the highest occupancy mode 
   * this overwrites the function of the same
   * name in BaseCDI
   */
  void scale_intensity(std::vector<Complex_2D> & c);

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

  /*
   * generate the S and J matrices for the decomposition 
   * of the partially coherent wave where JC=nSC where
   * H = integral(P*l(r1)J(r1, r2)Pm(r2)) dr1 dr2 and 
   * S=integral(P*l(r)pm(r))dr where Pl is an orhtonormal
   * basis set, in this case, the Legendre polynomials
   */
  void initialse_matrices(int leg, int modes);

  /**
   * the J matrix where J = integral(P*l(r1)J(r1, r2)Pm(r2))dr1dr2 
   * the x and y are computed seperately, then multiplied together.
   * The result is a matrix of xn+y by in+j where 
   */
  void fill_jmatrix(Double_2D legmatrix, Double_2D roots);

  /**
   * the S matrix = integral(P*l(r)pm(r))dr = 2/(2n+1)
   * from the orthogonality requirments of Legendre 
   * Polynomials. We then turn it in to a 2D matrix
   * for the x and y dimensions
   */
  void fill_smatrix(Double_2D legmatrix, Double_2D roots);

  /**
   * fill a vector of Complex_2D for single modes. These 
   * modes do not evolve over time, and so are not BaseCDI's
   */
  void fill_modes(Complex_2D & c);

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
   * calculate and return the current intensity of the modes multiplied
   * by the transmissoion fnction
   */
  Double_2D get_intensity();

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
    * Set the minimum contribution of a mode as a proportion of the 
    * dominant mode for it to be included in the reconstruction 
    */
  void set_threshold(double d){
    threshold=d;
  }

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

  /**
   *Propagate the CDI's 
   */
  void propagate_from_detector(Complex_2D & c);
  void propagate_to_detector(Complex_2D & c);
};


#endif
