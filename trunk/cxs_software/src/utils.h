#ifndef UTILS_H
#define UTILS_H

class Double_2D;
class Complex_2D;

void crop(Double_2D & image, Double_2D & new_image, int x_start, int y_start);
void rescale(Double_2D & image, double scale);

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

#endif

