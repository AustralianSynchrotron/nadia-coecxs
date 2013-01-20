#!/bin/bash
#
# Copyright 2011 Nadia Davidson 
# for The ARC Centre of Excellence in Coherent X-ray Science. 
#
# This program is distributed under the GNU General Public License. 
# We also ask that you cite this software in publications where you made 
# use of it for any part of the data analysis.
#
# This shows how to run the planar reconstuction
# using the command line tool. 

PATH=$PATH:../bin

#run the planar example:
CDI_reconstruction.exe planar_example.config

#########################################################

#to run fresnel reconstruction:
#reconstruct the white-field then the sample
#uncomment the lines below.

#CDI_reconstruction.exe fresnel_example.config fresnel_wf
#CDI_reconstruction.exe fresnel_example.config fresnel

#########################################################

#to run multiple times with a different starting seed do
#something like the following:

#for a in `seq 10`
#do
#  CDI_reconstruction.exe planar_example.config planar $a &> log_$a
#  mv planar.cplx planar_result_${a}_.cplx
#done  

#a tool for merging the output can then be used (not written yet).

