; This file demonstrates the steps in performing Phase Diverse
; FCDI reconstruction using the library in IDL. 
; Thanks for Corey for the data which this example replies on.

; Load up the module containing the wrapper code
.Compile CXS_interface.pro

; Load some image files of the white-field data and 
; zone-plate support. A 2D array of doubles is returned.
; Please replace the file-name with your 
; "cxs_software/example/image_files" directory if you are not
; running this example on osiris.

white_field = ["wf_A_1024.cplx", $
               "wf_B_1024.cplx", $
               "wf_C_1024.cplx", $
               "wf_D_1024.cplx", $
               "wf_E_1024.cplx", $
               "wf_F_1024.cplx", $
               "wf_G_1024.cplx"]

data = ["A.dbin", $
        "B.dbin", $
        "C.dbin", $
        "D.dbin", $
        "E.dbin", $
        "F.dbin", $
        "G.dbin" ]

wavelength = 4.892e-10 ; wavelength
fz = 16.353e-3  ; zone plate to focal distance
ps = 13.5e-6 ; pixel size

;zone plate to detector distances
fd = [0.909513, 0.909388, 0.909263, 0.909213, $
      0.909088, 0.909088, 0.909088] 

;sample to detector distances
fs = [18.513e-3, 18.388e-3, 18.263e-3, 18.213e-3, $
      18.088e-3, 18.088e-3, 18.088e-3]

;white-field normalisation
norm = [0.984729833,0.97700431,0.986270638, 0.967487825,$
        0.980945916, 0.97628279, 0.963066039 ]

;relative positions in pixels 
x_pos = [0,23,3,35,30,37,60]
y_pos = [-150,-145,-134,-123,-150,-93,210]


;make a support (use the same for all frames)
s = CXS_GET_ROUND_SUPPORT(1024,1024,0.5)

;Make a new PhaseDiverseCDI object
; CXS_INIT_FRESNEL
CXS_INIT_PHASE_DIVERSE, parallel=1


FOR I=0,1 DO BEGIN $
   d = CXS_READ_DBIN(1024,1024,data[I]) & $
   w = CXS_READ_CPLX(1024,1024,white_field[I]) & $

   CXS_INIT_FRESNEL, d, s, w, wavelength, fd[I]-fz, fs[I]-fz, ps, 1 & $
   ;norm[I] & $
   
   CXS_SET_CHARGE_FLIPPING, 1 & $
   CXS_SET_TRANS_UNITY_CONSTRAINT, 1 & $
   
   CXS_PHASE_DIVERSE_ADD_POSITION, x_pos[I], y_pos[I] & $

ENDFOR


a = CXS_PHASE_DIVERSE_INIT_ESTIMATE()

CXS_PHASE_DIVERSE_ADJUST_POSITIONS

a = CXS_PHASE_DIVERSE_ITERATE(10) 

CXS_PHASE_DIVERSE_ADJUST_POSITIONS, TYPE=1

b = CXS_PHASE_DIVERSE_ITERATE(5) 


; Now get the transmission function based on the result of the
; final iteration.
; a = cxs_get_transmission_function()
; cxs_clear_memory

; We are finished with the white-field reconstruction now, so
; let free some memory

; Now you can play with "a" however you like in IDL.

; e.g. get the phase and display it:
;phase = ATAN(b, /PHASE)
;window, XSIZE=512, YSIZE=512
;TVSCL, rebin(phase,512,512)

; or the magnitude:
; TVSCL, rebin(abs(a),512,512)

; or save the result to a file:
; cxs_write_cplx(a , 'result_of_my_FCDI_CDI.cplx')

