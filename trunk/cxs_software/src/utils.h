// Copyright 2011 Nadia Davidson for The ARC Centre of Excellence in 
// Coherent X-ray Science. This program is distributed under the GNU  
// General Public License. We also ask that you cite this software in 
// publications where you made use of it for any part of the data     
// analysis. 

#ifndef UTILS_H
#define UTILS_H

#include "Double_2D.h"
//class Double_2D;
class Complex_2D;

void crop(Double_2D & image, Double_2D & new_image, int x_start, int y_start);
void rescale(Double_2D & image, double scale);
void slow_align(Double_2D & image1, Double_2D & image2, int & offset_x, int & offset_y, int step_size=8.0, int min_x=0, int max_x=0, int min_y=0, int max_y=0);

void align(Double_2D & first_image, Double_2D & second_image,
	   int & offset_x, int & offset_y,
	   int min_x=0, int max_x=0,
	   int min_y=0, int max_y=0,
	   Double_2D * first_image_weights = 0,
	   Double_2D * second_image_weights = 0,
	   double overlap_fraction = 0.2);

double edges(Double_2D & image);
double line_out(Double_2D & image);
double calculate_high_frequency_ratio(Double_2D & image);
double calculate_image_entropy(Double_2D & image);
double calculate_gradients(Double_2D & image, double threshold=0.03);
double calculate_mean_difference(Double_2D & image);
double calculate_average_energy_density(Double_2D & image);
double vollaths_4(Double_2D & image);
double vollaths_5(Double_2D & image);
double deviation_from_zero(Double_2D & image);
double count_pixels(Double_2D & image, double threshold);
double diff_of_squares(Double_2D & image1, Double_2D & image2);
double simple(Double_2D & image, double scale);
double sobel_gradient(Double_2D & image);
double laplace_gradient(Double_2D & image);
double edge_grad(Double_2D & image, Double_2D & mask);
double calculate_image_entropy_2(Double_2D & image);

void interpolate( const Complex_2D & original, Complex_2D & big);
void interpolate( const Double_2D & original, Double_2D & big);
void shrink( const Complex_2D & original, Complex_2D & small);
void shrink( const Double_2D & original, Double_2D & small);


#endif

