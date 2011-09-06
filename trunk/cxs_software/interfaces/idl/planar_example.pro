; Copyright 2011 Nadia Davidson 
; for The ARC Centre of Excellence in Coherent X-ray Science. 
;
; This program is distributed under the GNU General Public License. 
; We also ask that you cite this software in publications where you made 
; use of it for any part of the data analysis.

; This file demonstrates the steps in performing planar CDI
; reconstruction using the library in IDL. 

; Load up the module containing the wrapper code
.Compile CXS_interface.pro

; Load some image files of the data and the support. 
; A 2D array of doubles is returned.
; Please replace the file-name with your 
; "cxs_software/example/image_files" directory if you are not
; running this example on osiris.
support = cxs_read_tiff(1024,1024,'../../examples/image_files/planar_support.tiff')
data = cxs_read_tiff(1024,1024,'../../examples/image_files/planar_data.tif')

; Set-up everything ready for planar reconstruction.
; You need to pass the image data and support.
; By default the HIO algorithm is used with a relaxation
; parameter of beta=0.9
cxs_init_planar, data, support

;Set the algorithm to hybrid input-out
cxs_set_algorithm, 'HIO'

; Perform 50 iterations and then 
; apply shrink wrap (default Gaussian width = 1.5 pixels,
; threshold = 0.1 time the maximum pixel value). 
; Do this 5 times (250 iterations in total).
FOR I=0,4 DO BEGIN a = cxs_iterate(50) & cxs_apply_shrinkwrap & ENDFOR

; Change to the error-reduction algorithm 
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

