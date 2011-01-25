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
;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;ok. complex_array option also tested
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

;ok
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

;ok
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
