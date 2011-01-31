;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Author: Nadia Davidson (nadiamd@unimelb.edu.au)
; Date: 24th January 2011
; 
; This code provides wrappers to the functions in the COECXS C++
; library. It calls the methods which are defined in the
; IDL_interface.c file (and compiled into libIDLCOECXS.so), thus
; allowing CDI reconstruction to be performed in IDL. Some examples
; are provided in this directory, showing how you can use these
; methods.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Couldn't work out how to make global variables
; so this is my way around it.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function lib_name
return, 'libIDLCOECXS.so'
end


function pixels
return, 512
end

function nx
return, call_external(lib_name(),'IDL_get_array_x_size')
end

function ny
return, call_external(lib_name(),'IDL_get_array_y_size')
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Here are the real wrappers to the C++ code
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;+
; NAME:
;       CXS_INIT_PLANAR
;
; PURPOSE:
;       Set-up a planar CDI reconstruction. This will
;       initialise the reconstruction with the data
;       and support. Some defaults will be set and
;       memory will be allocated ready for reconstruction.
;       It is necessary to call this procedure before
;       attempting to call any of the reconstruction methods
;       (e.g. before setting the algorithm or 
;       calling CXS_ITERATE).
;
;       Calling this procedure will initialise the reconstruction 
;       algorithm to hybrid-input-output with a relaxation parameter of 0.9.
;
;
; CALLING SEQUENCE:
;
;	CXS_INIT_PLANAR, Data, Support [,Starting_point]
;
;
; INPUTS:
;
;	Data: 
;             The detector illumination. It should be
;             in the form of a 2D array
;
;       Support: 
;             A 2D Array giving the sample support.
;             Values or 1 or greater are considered inside
;             the support. All others are considered to be
;             outside the support.
;
;       Starting_point: 
;             As an option you may supply an initial 
;             guess of the exit-surface-wave for the sample. 
;             This maybe useful, for example, if you wish to 
;             start from the end point of a previous run. The
;             format of this parameter much be a 2D array of
;             COMPLEX variables. If this parameter is not supplied,
;             the starting point is initialised to be zero outside
;             the support and a random number inside the support, 
;             for both the magnitude and phase.
;
; EXAMPLE:
;        An example of loading two 2D arrays from file and using
;        them to initialise the planar reconstruction:
;
;        my_support = cxs_read_tiff(1024,1024,'planar_support.tiff')
;        my_data = cxs_read_tiff(1024,1024,'planar_data.tiff')
;        CXS_INIT_PLANAR, my_data, my_support
;
;-
pro cxs_init_planar, data, support, complex_array
n = size(data)
nx = n[2]
ny = n[1]
IF N_Params() EQ 3 THEN $
  b = call_external(lib_name(),'IDL_planar_init',nx,ny,complex_array) $
ELSE $
  b = call_external(lib_name(),'IDL_planar_init',nx,ny)
cxs_set_support, support
cxs_set_intensity, data
IF N_Params() EQ 2 THEN $
  cxs_initialise_esw
end


;+
; NAME:
;       CXS_INIT_FRESNEL_WF
;
; PURPOSE:
;       Set-up a Fresnel white-field CDI reconstruction. This will
;       initialise the reconstruction with the white-field intensity, 
;       zone-plate support and experimental parameters. Some defaults 
;       will be set and memory will be allocated ready for 
;       reconstructing the white-field (phase and magnitude) in the
;       detector plane. It is necessary to call this procedure before
;       attempting to call any of the reconstruction methods
;       (e.g. before setting the algorithm or calling CXS_ITERATE).
;
; CALLING SEQUENCE:
;
;	CXS_INIT_FRESNEL_WF, data, support, beam_wavelength,
;                            zone_focal_length, focal_detector_length,
;                            pixel_size [,starting_point]
;
; INPUTS:
;
;	data: 
;             The white-field illumination. It should be
;             in the form of a 2D array.
;
;       support: 
;             A 2D array of integers or doubles giving the zone-plate support.
;             Values or 1 or greater are considered inside
;             the support. All others are considered to be
;             outside the support.
;
;       beam_wavelength:
;             The beam wavelength.
;
;       zone_focal_length:
;             The distance between the zone plate and the focal point.
;
;       focal_detector_length:
;             The distance between the focal point and the detector.
;
;       pixel_size:
;             The side length of one detector pixel.
;
;       starting_point: 
;             As an option you may supply an initial guess of the 
;             white-field. This may be useful, for example, if you 
;             wish to start from the end point of a previous run. The
;             format of this parameter must be a 2D array of
;             COMPLEX variables. If this parameter is not supplied,
;             the starting point is initialised to be zero outside
;             the support, a random number inside the support 
;             for the magnitude and zero for the phase.
;
; EXAMPLE:
;
;        An example of loading two 2D arrays from file and using
;        them to initialise the white-field reconstruction for FCDI:
;
;        my_support = cxs_read_tiff(1024,1024,'support.tiff')
;        my_data = cxs_read_tiff(1024,1024,'data.tiff')
;        cxs_init_fresnel_wf, my_data, my_support, 4.892e-10, 16.353e-3, 0.9078777,13.5e-6
;-
pro cxs_init_fresnel_wf, data, $
                         support, $
                         beam_wavelength, $
                         zone_focal_length, $
                         focal_detector_length, $
                         pixel_size, $
                         complex_array
n = size(data)
nx = n[2]
ny = n[1]
IF N_Params() EQ 7 THEN $
  b = call_external(lib_name(),'IDL_fresnel_wf_init',nx,ny, $
                    double(beam_wavelength), $
                    double(zone_focal_length), $
                    double(focal_detector_length), $
                    double(pixel_size), $
                    complex_array) $
ELSE $
  b = call_external(lib_name() ,'IDL_fresnel_wf_init',nx,ny, $
                    double(beam_wavelength), $
                    double(zone_focal_length), $
                    double(focal_detector_length), $
                    double(pixel_size))

cxs_set_support, support
cxs_set_intensity, data

IF N_Params() EQ 6 THEN $
  cxs_initialise_esw

end

;+
; NAME:
;       CXS_INIT_FRESNEL
;
; PURPOSE:
;       Set-up a Fresnel CDI reconstruction. This will
;       initialise the reconstruction using a previously reconstructed
;       white-field, detector data, sample support and experimental 
;       parameters. Some defaults will be set and memory will be
;       allocated ready for reconstructing the sample
;       exit-surface-wave. It is necessary to call this procedure before
;       attempting to call any of the reconstruction methods
;       (e.g. before setting the algorithm or calling CXS_ITERATE).
;
;       Calling this procedure will initialise the reconstruction algorithm
;       to the error-reduction with a relaxation parameter of 0.9.
;
; CALLING SEQUENCE:
;
;	CXS_INIT_FRESNEL, data, support, white-field, beam_wavelength,
;	                  focal_detector_length, focal_sample_length, 
;                         pixel_size [, normalisation, starting_point ]
;
; INPUTS:
;
;	data: 
;             The detector data with the sample in place. It should be
;             in the form of a 2D array.
;
;       support: 
;             A 2D array of integers or doubles which give the sample support.
;             Values or 1 or greater are considered inside
;             the support. All others are considered to be
;             outside the support.
;
;       white-field:
;             A COMPLEX 2D array of the reconstructed white-field in
;             the detector plane. This can be recovered using 
;             CXS_INIT_FRESNEL_WF followed by CXS_ITERATE.
;
;       beam_wavelength:
;             The beam wavelength.
;
;       focal_detector_length:
;             The distance between the focal point and the detector.
;
;       focal_sample_length:
;             The distance between the focal point and the sample.
;
;       pixel_size:
;             The side length of one detector pixel.
;
;       normalisation: 
;             The factor to scale the white-field before
;             performing FCDI. If this parameter is excluded, the
;             ratio of the square-root of the intensity data and the
;             white-field magnitude is used as the normalisation.
;
;       starting_point: 
;             As an option you may supply an initial 
;             guess of the exit-surface-wave for the sample. 
;             This maybe useful, for example, if you wish to 
;             start from the end point of a previous run. The
;             format of this parameter much be a 2D array of
;             COMPLEX variables. If this parameter is not supplied,
;             the initialisation described in Harry's review paper: 
;             page 29. (in particular e.q. 137) is used.
;
; EXAMPLE:
;
;        cxs_init_fresnel, my_data, my_supports, my_white-field, $
;                          4.892e-10, 0.9078777, 2.16e-3, $
;                          13.5e-6
;-
pro cxs_init_fresnel, data, support, $
                      white_field, $
                      beam_wavelength, $
                      focal_detector_length, $
                      focal_sample_length, $
                      pixel_size, $
                      normalisation,$
                      complex_array
n = size(data)
nx = n[2]
ny = n[1]
IF N_Params() EQ 7 THEN BEGIN
  mag2_wf = abs(white_field)^2
  normalisation = total(data) / total(mag2_wf)
  print, normalisation
ENDIF
IF N_Params() EQ 9 THEN $
  b = call_external(lib_name(),'IDL_fresnel_init',nx,ny, $
                    white_field, $
                    double(beam_wavelength), $
                    double(focal_detector_length), $
                    double(focal_sample_length), $
                    double(pixel_size), $
                    double(normalisation),$
                    complex_array) $
ELSE $
  b = call_external(lib_name() ,'IDL_fresnel_init',nx,ny, $
                    white_field, $
                    double(beam_wavelength), $
                    double(focal_detector_length), $
                    double(focal_sample_length), $
                    double(pixel_size), $
                    double(normalisation))

cxs_set_support, support
cxs_set_intensity, data

IF N_Params() LT 9 THEN $
  cxs_initialise_esw
end


;+
; NAME:
;       CXS_SET_SUPPORT
;
; PURPOSE:
;       Set the support shape to be used in reconstruction. This with override 
;       the support given to any of the CXS_INIT methods and maybe called
;       at any time during the reconstruction.  
;
; CALLING SEQUENCE:
;
;	CXS_SET_SUPPORT, support
;
; INPUTS:
;
;       support: 
;             A 2D array of doubles or integers giving the sample's 
;             (or zone-plate's) support. Values of 1 or greater are 
;             considered inside the support. All others are considered 
;             to be outside the support.
;
; EXAMPLE:
;
;       cxs_set_support, my_support
;
;-
pro cxs_set_support, array
n = size(array)
a = call_external(lib_name(),'IDL_set_support',n[1],n[2],double(array))
end

;+
; NAME:
;       CXS_SET_INTENSITY
;
; PURPOSE:
;       Set the detector intensity data. This will override 
;       the intensity given to any of the CXS_INIT methods.
;       In general, users should not need to call this method.
;
; CALLING SEQUENCE:
;
;	CXS_SET_INTENSITY, data
;
; INPUTS:
;
;       data: 
;             The detector data. It should be in the form of a 2D
;             array or doubles or integers.
;
;
; EXAMPLE:
;
;       cxs_set_intensity, my_data
;
;-
pro cxs_set_intensity, array
n = size(array)
a= call_external(lib_name(),'IDL_set_intensity',n[1],n[2],double(array))
end

;+
; NAME:
;       CXS_INITIALISE_ESW
;
; PURPOSE:
;       Initialise the exit-surface-wave guess. The initialisation
;       will depend on the reconstruction type (see the procedures
;       which begin "CXS_INIT_" for a description). The "INIT" procedure  
;       will call this procedure if no starting guess is provided. 
;       It is useful if you wish to run the same reconstruction
;       several time with a different random starting point, or if 
;       you wish to reset the reconstruction back to the original guess.
;
;
; CALLING SEQUENCE:
;
;	CXS_INITIALISE_ESW, seed
;
; INPUTS:
;
;       seed: 
;             Seed for the random number generator used to initialise
;             the guess. It should be an integer. This is ignored in 
;             the case of Fresnel CDI.
;
;
; EXAMPLE:
;
;       cxs_initialise_esw, 6
;-
pro cxs_initialise_esw, seed
IF N_Params() EQ 0 THEN $
  seed = 0
b = call_external(lib_name() ,'IDL_initialise_esw',long(seed)) 
end

;+
; NAME:
;       CXS_ITERATE
;
; PURPOSE:
;       Perform the iterative reconstruction. The number of iterations
;       to perform should be given and the result of the final
;       iteration is returned. For planar and Fresnel CDI this will be
;       the exit surface wave of the sample. For Fresnel white-field
;       reconstruction it will be the white-field at the detector surface.
;       The magnitude of the result is also displayed on the screen as
;       a 512x512 pixel image. The iteration number and corresponding
;       error (see CXS_GET_ERROR) will be printed on the screen.
;
;       CXS_ITERATE may be called several time and the reconstruction
;       will start from where is ended. i.e. calling CXS_ITERATE(50),
;       followed by a second CXS_ITERATE(50) is equivalent to calling
;       CXS_ITERATE(100).
;
;       During the reconstruction, the best (lowest error) result is
;       also stored. It maybe retrieved by calling CXS_GET_BEST_RESULT.
; 
;
; CALLING SEQUENCE:
;
;	result = CXS_ITERATE([iterations])
;
; INPUTS:
;
;       iterations: 
;             The number of iterations to perform (an integer). 
;             If left empty, one iteration is performed.
;
; OUTPUTS:
;
;       result:
;             A COMPLEX 2D array. For planar and Fresnel CDI this will be
;             the exit surface wave of the sample. For Fresnel white-field
;             reconstruction it will be the white-field at the detector surface.
;
;
; EXAMPLE:
;
;       my_result = cxs_iterate(100)
;-
function cxs_iterate, iterations
nx = call_external(lib_name(),'IDL_get_array_x_size')
ny = call_external(lib_name(),'IDL_get_array_y_size')
result = make_array(nx,ny,/COMPLEX)
IF N_Params() EQ 0 THEN $
   iterations = 1
b = call_external(lib_name() ,'IDL_iterate',long(iterations),result) 
window, XSIZE=pixels(), YSIZE=pixels()
TVSCL, rebin(ABS(result),pixels(),pixels())
return, result
end

;+
; NAME:
;       CXS_SET_RELAXATION_PARAMETER
;
; PURPOSE:
;       Set the relaxation parameter. The default relaxation
;       parameter used in Planar and Fresnel reconstruction is 0.9.
;       In Fresnel white-field reconstruction, this parameter is 
;       not used.
;
; CALLING SEQUENCE:
;
;	CXS_SET_RELAXATION_PARAMETER, beta
;
; INPUTS:
;
;       beta: 
;             The relaxation parameter.
;
; EXAMPLE:
;
;       cxs_set_relaxation_parameter, 0.9
;-
pro cxs_set_relaxation_parameter, beta
b = call_external(lib_name() ,'IDL_set_relaxation_parameter',double(beta)) 
end

;+
; NAME:
;       CXS_APPLY_SHRINKWRAP
;
; PURPOSE:
;       Apply the shrinkwrap algorithm. The current exit-surface-wave
;       magnitude is used to update the support; it is convoluted with
;       a Gaussian and then a threshold is applied. You can use the
;       cxs_get_support function to see how the support have been
;       modified after calling this procedure.
;
; CALLING SEQUENCE:
;
;	CXS_APPLY_SHRINKWRAP [,gauss_width, threshold ]
;
; INPUTS:
;
;       gauss_width: 
;             The width (1-standard deviation. in pixels) of the Gaussian
;             used for smearing. If this parameter is not passed, 
;             a width of 1.5 pixels is used.
;
;       threshold:
;             All pixels which are below the threshold are set to zero 
;             (outside the support). The threshold should be given as
;             a fraction of the maximum pixel value (in double format). 
;             If this parameter is not passed a threshold of 0.1 (10%)
;             is used.
;
; EXAMPLE:
;       Perform 1000 iterations in total, applying shrink-wrap at the 400th iteration:
;
;       ....
;       a = CXS_ITERATE(400)
;       CXS_APPLY_SHRINKWRAP
;       a = CXS_ITERATE(600)
;-
pro cxs_apply_shrinkwrap, gauss_width, threshold
if N_Params() lt 2 then threshold = 0.1 
if N_Params() lt 1 then gauss_width = 1.5
b = call_external(lib_name() ,'IDL_apply_shrinkwrap',double(gauss_width),double(threshold))
end

;+
; NAME:
;       CXS_SET_ALGORITHM
;
; PURPOSE:
;       Select the reconstruction algorithm to use. Options are:
;       'ER' - error reduction 
;       'BIO' - basic input-output 
;       'BOO' - basic output-output 
;       'HIO' - hybrid input-output 
;       'DM' - difference map 
;       'SF' - solvent-flipping 
;       'ASR' - averaged successive reflections 
;       'HPR' - hybrid projection reflection 
;       'RAAR' - relaxed averaged alternating reflectors
;
; CALLING SEQUENCE:
;
;       CXS_SET_ALGORITHM, algorithm
;
; INPUTS:
;
;       algorithm:
;             Set the reconstruction algorithm. It should be
;             one of the options listed above. Note that it is
;             passed as a string.
;
; EXAMPLE:
;       Perform 1000 iterations in total, changing to error-reduction at the 400th iteration:
;       ....
;       a = CXS_ITERATE(400)
;       CXS_SET_ALGORITHM, 'ER'
;       a = CXS_ITERATE(600)
;-
pro cxs_set_algorithm, algorithm
b = call_external(lib_name() ,'IDL_set_algorithm',algorithm)
end

;+
; NAME:
;       CXS_SET_CUSTOM_ALGORITHM
;
; PURPOSE:
;       Set a custom reconstruction algorithm.
;       
;       Iterative reconstruction algorithms can be expressed as a
;       combination of several operators. The parameters to this
;       procedure set the coefficients for combinations of these
;       operators. For a description, please see page 22 of the
;       H.M. Quiney review: TUTORIAL REVIEW, Coherent diffractive
;       imaging using short wavelength light sources, Journal of
;       Modern Optics, 2010, DOI: 10.1080/09500340.2010.495459.  Some
;       description is also given in the C++ doxygen documentation
;       for PlanarCDI::set_custom_algorithm.
;
;
; CALLING SEQUENCE:
;
;       CXS_SET_CUSTOM_ALGORITHM, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10
;
; INPUTS:
;
;       m1-m10:
;             Coefficients to the operator combinations.
;
; EXAMPLE:
;       Perform 1000 iterations in total, changing to a custom
;       algorithm at the 400th iteration and print the algorithm to screen:

;       ....
;       a = CXS_ITERATE(400)
;       CXS_SET_ALGORITHM, 0.5, 0.5, 0.5, 0.5, 0.5, $ 
;                          0.5, 0.5, 0.5, 0.5, 0.5
;       CXS_PRINT_ALGORITHM                    
;       a = CXS_ITERATE(600)
;-
pro cxs_set_custom_algorithm, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10
b = call_external(lib_name() ,'IDL_set_custom_algorithm', $
                  double(m1), $
                  double(m2), $
                  double(m3), $ 
                  double(m4), $ 
                  double(m5), $
                  double(m6), $
                  double(m7), $
                  double(m8), $
                  double(m9), $
                  double(m10) )
end


;+
; NAME:
;       CXS_GET_BEST_RESULT
;
; PURPOSE:
;       Get the best (lowest error) result found during the reconstruction.
;       The result is returned and the magnitude of the result is
;       displayed as a 512x512 pixel image on the screen.
;
; CALLING SEQUENCE:
;
;       result = CXS_GET_BEST_RESULT()
;
; OUTPUTS:
;
;       result:
;             A 2D array of COMPLEX variables. For planar and Fresnel
;             reconstruction this will be the  exit-surface-wave for
;             the sample. For Fresnel white-field reconstruction, this
;             will be the complex white-field in the detector plane.
;
;
; EXAMPLE:
;       Perform 1000 iterations and get the lowest error result.
;       ....
;       a = CXS_ITERATE(1000)
;       a = CXS_GET_BEST_RESULT()
;-
function cxs_get_best_result 
result = make_array(nx(),ny(),/COMPLEX)
b = call_external(lib_name() ,'IDL_get_best_result',result)
window, XSIZE=pixels(), YSIZE=pixels()
TVSCL, rebin(ABS(result),pixels(),pixels())
return, result
end

;+
; NAME:
;       CXS_GET_SUPPORT
;
; PURPOSE:
;       Get the support. This maybe useful to see how shrinkwrap has
;       effected the support. The result is displayed as a 512x512 
;       pixel image on the screen.
;
; CALLING SEQUENCE:
;
;       support = CXS_GET_SUPPORT()
;
; OUTPUTS:
;
;       support:
;             The support. A 2D array of doubles. A pixel with value below 1
;             is outside the support and 1 and above is inside the support.
;
; EXAMPLE:
;       View the support after applying shrink-wrap.
; 
;       CXS_APPLY_SHRINKWRAP
;       a = CXS_GET_SUPPORT()
;-
function cxs_get_support
result = make_array(nx(),ny(),/DOUBLE)
b = call_external(lib_name() ,'IDL_get_support',result) 
window, XSIZE=pixels(), YSIZE=pixels()
TVSCL, rebin(result,pixels(),pixels())
return, result
end

;+
; NAME:
;       CXS_GET_ERROR
;
; PURPOSE:
;       Get the error metric. This is defined as the difference
;       between the estimated diffraction and the actual diffraction
;       pattern and is calculated as:
;           sum over pixels of ( M - sqrt(I) ) ^2 / sum(I) 
;       Where I is the detector data and M is the magnitude of the
;       estimate in the detector plane. Note that due to the way this
;       quantity is calculated, it actually corresponds to the
;       previous estimate rather than the current iteration.
;
;
; CALLING SEQUENCE:
;
;       error = CXS_GET_ERROR()
;
; OUTPUTS:
;
;       error:
;             The error. 
;
; EXAMPLE:
;
;       a = CXS_GET_ERROR()
;- 
function cxs_get_error
result = double(0.0)
b = call_external(lib_name(),'IDL_get_error',result)
return, result
end

;+
; NAME:
;       CXS_GET_TRANSMISSION_FUNCTION
;
; PURPOSE:
;       Get the transmission function from the current estimate of the
;       exit-surface-wave. The magnitude of the result will be
;       displayed on the screen. Note that this functions
;       is only available for Fresnel CDI reconstruction.
;
; CALLING SEQUENCE:
;
;       result = CXS_GET_TRANSMISSION_FUNCTION()
;
; OUTPUTS:
;
;       result:
;             The transmission function for the sample. A 2D array of
;             COMPLEX variables is returned. 
;
; EXAMPLE:
;       a = CXS_GET_TRANSMISSION_FUNCTION()
;-
function cxs_get_transmission_function
result = make_array(nx(),ny(),/COMPLEX)
b = call_external(lib_name() ,'IDL_get_transmission_function',result) 
window, XSIZE=pixels(), YSIZE=pixels()
TVSCL, rebin(abs(result),pixels(),pixels())
return, result
end

;+
; NAME:
;       CXS_CLEAR_MEMORY
;
; PURPOSE:
;       Clean-up after a reconstruction has been performed. This 
;       procedure should be called at the very end of a program.
;       It will free up the memory that was allocated when one of 
;       the "CXS_INIT_.." methods was called.
;
; CALLING SEQUENCE:
;
;       CXS_CLEAR_MEMORY
;
; EXAMPLE:
;       CXS_CLEAR_MEMORY
;-
pro cxs_clear_memory
b = call_external(lib_name() ,'IDL_deallocate_memory')
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;function get_result_at_sample


;function get_result_at_detector


;function get_result_at_zone_plate



;+
; NAME:
;       CXS_GET_INTENSITY_AUTOCORRELATION
;
; PURPOSE:
;       Get the autocorrelation function from the intensity data.
;       This method is only useful for planar CDI reconstruction.
;
; CALLING SEQUENCE:
;
;       a = CXS_GET_INTENSITY_AUTOCORRELATION()
;
; OUTPUTS:
;       a:
;             The autocorrelation function (a 2D array of doubles).
;
; EXAMPLE:
;       a = CXS_GET_INTENSITY_AUTOCORRELATION()
;-
function cxs_get_intensity_autocorrelation
result = make_array(nx(),ny(),/DOUBLE)
b = call_external(lib_name() ,'IDL_get_intensity_autocorrelation',result) 
window, XSIZE=pixels(), YSIZE=pixels()
TVSCL, rebin(result,pixels(),pixels())
return, result
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; io functions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;ok
function cxs_read_dbin, nx, ny, filename
result = make_array(nx,ny, /DOUBLE)
b = call_external(lib_name() ,'IDL_read_dbin',nx,ny,filename,result)
window, XSIZE=pixels(), YSIZE=pixels()
TVSCL, rebin(result,pixels(),pixels())
return, result
end

;ok
function cxs_read_tiff, nx, ny, filename
result = make_array(nx,ny, /DOUBLE)
b = call_external(lib_name() ,'IDL_read_tiff',nx,ny,filename,result)
window, XSIZE=pixels(), YSIZE=pixels()
TVSCL, rebin(result,pixels(),pixels())
return, result
end

;ok
function cxs_read_cplx, nx, ny, filename
result = make_array(nx,ny, /COMPLEX)
b = call_external(lib_name() ,'IDL_read_cplx',nx,ny,filename,result)
window, XSIZE=pixels(), YSIZE=pixels()
TVSCL, rebin(abs(result),pixels(),pixels())
return, result
end

;ok
pro cxs_write_cplx, complex_array, filename
n = size(complex_array)
b = call_external(lib_name() ,'IDL_write_cplx',n[1],n[2],complex_array, filename)
end

;ok
pro cxs_write_dbin, array, filename
n = size(array)
b = call_external(lib_name() ,'IDL_write_dbin',n[1],n[2],array, filename)
end

;ok
pro cxs_print_algorithm
b = call_external(lib_name() ,'IDL_print_algorithm')
end