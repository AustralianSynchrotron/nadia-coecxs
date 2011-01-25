

.Compile CXS_interface.pro


s = cxs_read_tiff(1024,1024,'/data/nadia/cxs_software_rel_0/cxs_software/examples/image_files/FCDI_wf_support.tiff')

d = cxs_read_dbin(1024,1024,'/data/nadia/cxs_software_rel_0/cxs_software/examples/image_files/FCDI_wf_data.dbin')

cxs_fresnel_wf_init, d, s, 4.892e-10, 16.353e-3, 0.9078777, 13.5e-6

a = cxs_iterate(20)

s = cxs_read_tiff(1024,1024,'/data/nadia/cxs_software_rel_0/cxs_software/examples/image_files/FCDI_support.tiff')

d = cxs_read_dbin(1024,1024,'/data/nadia/cxs_software_rel_0/cxs_software/examples/image_files/FCDI_data.dbin')

cxs_clear_memory


i = cxs_read_cplx(1024,1024,'my_fresnel.cplx')

cxs_fresnel_init, d, s, a, 4.892e-10, 0.9078777, 2.16e-3, 13.5e-6, 0.95, i
cxs_set_algorithm, 'ER'
cxs_print_algorithm

a = cxs_iterate(20)
a = cxs_get_transmission_function()

phase = ATAN(a, /PHASE)
window, XSIZE=pixels(), YSIZE=pixels()
TVscl, rebin(phase,pixels(),pixels())

;mag = abs(a)
;mag = (1-mag)*1000
;window, XSIZE=pixels(), YSIZE=pixels()
;TV, rebin(mag,pixels(),pixels())


