.Compile CXS_interface.pro

support = cxs_read_tiff(1024,1024,'../../examples/image_files/planar_support.tiff')
data = cxs_read_tiff(1024,1024,'../../examples/image_files/poly_sim_intensity.tiff')


cxs_init_poly, data, support 

cxs_set_support, support

cxs_set_spectrum, '../../examples/image_files/spectrum.txt'

cxs_set_algorithm, 'HIO'

FOR I=0,4 DO BEGIN a = cxs_iterate(10) & cxs_apply_shrinkwrap & ENDFOR


cxs_set_algorithm, 'ER'

; Do another 150 iterations with shrink-wrap applied every 50 iterations
FOR I=0,2 DO BEGIN a = cxs_iterate(50) & cxs_apply_shrinkwrap & ENDFOR

; Do one last iteration to get the final result
a = cxs_iterate()

; Apply shrink-wrap again. Lets be less tight this time
; cxs_apply_shrinkwrap, 1, 0.05 

; Do one more iteration
;a = cxs_iterate()


; Lets get the result with the lowest error from the
; last 100 iterations:
;a = cxs_get_best_result()

; and use this to restart the reconstruction starting with this
; result:
;cxs_init_planar, data, support, a

; We have finished with the reconstruction now, so free-up the
; memory we allocated earlier (note this does not effect "a").
; You will not longer be able to call cxs_iterate or any of the
; cxs get and set methods.
cxs_clear_memory

; Now you can play with "a" however you like in IDL.

; e.g. get the phase and display it:
; phase = ATAN(a, /PHASE)
; window, XSIZE=512, YSIZE=512
; TVSCL, rebin(phase,512,512)

; or the magnitude:
; TVSCL, rebin(abs(a),512,512)

; or save the result to a file:
; write_cplx(a , 'result_of_my_planar_CDI.cplx')

