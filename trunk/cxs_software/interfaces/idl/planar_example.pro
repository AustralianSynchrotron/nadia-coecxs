

.Compile CXS_interface.pro

s = cxs_read_tiff(1024,1024,'/data/nadia/cxs_software_rel_0/cxs_software/examples/image_files/planar_support.tiff')
d = cxs_read_dbin(1024,1024,'/data/nadia/cxs_software_rel_0/cxs_software/examples/image_files/planar_data.dbin')

cxs_planar_init, d, s

;cxs_set_algorithm, 'HIO'
a = cxs_iterate(5)
cxs_apply_shrinkwrap
a = cxs_iterate(50)
cxs_apply_shrinkwrap
cxs_set_algorithm, 'ER'
a = cxs_iterate(50)
cxs_apply_shrinkwrap
a = cxs_iterate(50)

