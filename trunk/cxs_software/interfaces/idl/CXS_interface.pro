;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Author: Nadia Davidson (nadiamd@unimelb.edu.au)
; Date: 24th January 2011
; 
; This code provides wrappers to the functions in the COECXS C++
; librarie. It calls the methods which are defined in the
; IDL_interace.c file (and compiled into libIDLCOECXS.so), thus
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
;       CXS_PLANAR_INIT
;
; PURPOSE:
;       Set-up a planar CDI reconstuction. This will
;       initialise the reconstruction with the data
;       and support. Some defaults will be set and
;       memory will be allocated ready for reconstruction.
;       It is necessary to call this procedure before
;       attempting to call any of the reconstruction methods
;       (e.g. before setting the algorithm or 
;       calling CXS_ITERATE).
;
; CALLING SEQUENCE:
;
;	CXS_PLANAR_INIT, Data, Support [,Starting_point]
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
;             for both the mangitude and phase.
;
; EXAMPLE:
;        An example of loading two 2D arrays from file and using
;        them to initialise the planar reconstruction:
;
;        my_support = cxs_read_tiff(1024,1024,'planar_support.tiff')
;        my_data = cxs_read_tiff(1024,1024,'planar_data.tiff')
;        CXS_PLANAR_INIT, my_data, my_support
;
;-
pro cxs_planar_init, data, support, complex_array
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
;       CXS_FRESNEL_WF_INIT
;
; PURPOSE:
;       Set-up a Fresnel white-field CDI reconstuction. This will
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
;	CXS_FRESNEL_WF_INIT, data, support, beam_wavelength,
;                            zone_focal_length, focal_detector_length,
;                            pixel_size [,starting_point]
;
; INPUTS:
;
;	data: 
;             The detector illumination. It should be
;             in the form of a 2D array.
;
;       support: 
;             A 2D Array giving the sample support.
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
;             As an option you may supply an initial 
;             guess of the exit-surface-wave for the sample. 
;             This maybe useful, for example, if you wish to 
;             start from the end point of a previous run. The
;             format of this parameter much be a 2D array of
;             COMPLEX variables. If this parameter is not supplied,
;             the starting point is initialised to be zero outside
;             the support and a random number inside the support, 
;             for both the mangitude and phase.
;
; EXAMPLE:
;
;        An example of loading two 2D arrays from file and using
;        them to initialise the planar reconstruction:
;
;        my_support = cxs_read_tiff(1024,1024,'support.tiff')
;        my_data = cxs_read_tiff(1024,1024,'data.tiff')
;        cxs_fresnel_wf_init, my_data, my_support, 4.892e-10, 16.353e-3, 0.9078777,13.5e-6
;-
pro cxs_fresnel_wf_init, data, $
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
;       CXS_FRESNEL_INIT
;
; PURPOSE:
;       Set-up a Fresnel CDI reconstuction. This will
;       initialise the reconstruction using a previously reconstructed
;       white-field, detector data, sample support and experimental 
;       parameters. Some defaults will be set and memory will be
;       allocated ready for reconstructing the sample
;       exit-surface-wave. It is necessary to call this procedure before
;       attempting to call any of the reconstruction methods
;       (e.g. before setting the algorithm or calling CXS_ITERATE).
;
; CALLING SEQUENCE:
;
;	CXS_FRESNEL_INIT, data, support, white-field, beam_wavelength,
;	                  focal_detector_length, focal_sample_length, 
;                         pixel_size [, normalisation, starting_point ]
;
; INPUTS:
;
;	data: 
;             The detector illumination. It should be
;             in the form of a 2D array.
;
;       support: 
;             A 2D Array giving the sample support.
;             Values or 1 or greater are considered inside
;             the support. All others are considered to be
;             outside the support.
;
;       white-field:
;
;
;
;       beam_wavelength:
;             The beam wavelength.
;
;
;       focal_detector_length:
;             The distance between the focal point and the detector.
;
;       pixel_size:
;             The side length of one detector pixel.
;
;       starting_point: 
;             As an option you may supply an initial 
;             guess of the exit-surface-wave for the sample. 
;             This maybe useful, for example, if you wish to 
;             start from the end point of a previous run. The
;             format of this parameter much be a 2D array of
;             COMPLEX variables. If this parameter is not supplied,
;             the starting point is initialised to be zero outside
;             the support and a random number inside the support, 
;             for both the mangitude and phase.
;
; EXAMPLE:
;
;        An example of loading two 2D arrays from file and using
;        them to initialise the planar reconstruction:
;
;        my_support = cxs_read_tiff(1024,1024,'support.tiff')
;        my_data = cxs_read_tiff(1024,1024,'data.tiff')
;        cxs_fresnel_wf_init, my_data, my_support, 4.892e-10, 16.353e-3, 0.9078777,13.5e-6
;-
pro cxs_fresnel_init, data, support, $
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


;ok
pro cxs_set_support, array
n = size(array)
a = call_external(lib_name(),'IDL_set_support',n[1],n[2],double(array))
end

;ok
pro cxs_set_intensity, array
n = size(array)
a= call_external(lib_name(),'IDL_set_intensity',n[1],n[2],double(array))
end

;ok
pro cxs_initialise_esw, seed
IF N_Params() EQ 0 THEN $
  seed = 0
b = call_external(lib_name() ,'IDL_initialise_esw',long(seed)) 
end

;ok
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

;ok
pro cxs_set_relaxation_parameter, beta
b = call_external(lib_name() ,'IDL_set_relaxation_parameter',double(beta)) 
end

;ok
pro cxs_apply_shrinkwrap, gauss_width, threshold
if N_Params() lt 2 then threshold = 0.1 
if N_Params() lt 1 then gauss_width = 1.5
b = call_external(lib_name() ,'IDL_apply_shrinkwrap',double(gauss_width),double(threshold))
end

;ok
pro cxs_set_algorithm, algorithm
b = call_external(lib_name() ,'IDL_set_algorithm',algorithm)
end

;ok
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


;ok
function cxs_get_best_result 
result = make_array(nx(),ny(),/COMPLEX)
b = call_external(lib_name() ,'IDL_get_best_result',result)
window, XSIZE=pixels(), YSIZE=pixels()
TVSCL, rebin(ABS(result),pixels(),pixels())
return, result
end

;ok
function cxs_get_support
result = make_array(nx(),ny(),/DOUBLE)
b = call_external(lib_name() ,'IDL_get_support',result) 
window, XSIZE=pixels(), YSIZE=pixels()
TVSCL, rebin(result,pixels(),pixels())
return, result
end

;ok (but need to check weather HIO is giving the correct
;error
function cxs_get_error
result = double(0.0)
b = call_external(lib_name(),'IDL_get_error',result)
return, result
end

;ok
function cxs_get_transmission_function
result = make_array(nx(),ny(),/COMPLEX)
b = call_external(lib_name() ,'IDL_get_transmission_function',result) 
window, XSIZE=pixels(), YSIZE=pixels()
TVSCL, rebin(abs(result),pixels(),pixels())
return, result
end

;ok
pro cxs_clear_memory
b = call_external(lib_name() ,'IDL_deallocate_memory')
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;function get_result_at_sample


;function get_result_at_detector


;function get_result_at_zone_plate




;ok
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
