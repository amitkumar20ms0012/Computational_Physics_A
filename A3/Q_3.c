#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

// Macro functions to extract real and imaginary parts
#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

// Function definition for the given function f(x)
double f(double x) {
    if (x == 0.0) {
        return 1.0;
    } else {
        return (sin(x) / x);
    }
}

int main() {
    int n = 1024; // Sampling points    
    double x_min = -50.0;
    double x_max = 50.0;
    double x_arr[n], k_arr[n], fact_q[2 * n], f_data[n], ft_f_data[2 * n];
    FILE *ft_data;
    double d = (x_max - x_min) / (n - 1);

    // Open file for writing
    ft_data = fopen("Q_3_fft_gsl_data.txt", "w");

    // Generate sample points and initialize arrays
    for (int i = 0; i < n; i++) {
        x_arr[i] = x_min + i * d;
        if (i < n / 2)
            k_arr[i] = 2 * M_PI * (i / (n * d));
        else
            k_arr[i] = 2 * M_PI * ((i - n) / (n * d));
        f_data[i] = f(x_arr[i]);
        REAL(ft_f_data, i) = f(x_arr[i]);
        IMAG(ft_f_data, i) = 0.0;
        REAL(fact_q, i) = cos(k_arr[i] * x_min);
        IMAG(fact_q, i) = -sin(k_arr[i] * x_min);
    }

    // Perform FFT
    gsl_fft_complex_radix2_transform(ft_f_data, 1, n, +1);

    // Normalize
    for (int i = 0; i < n; i++) {
        REAL(ft_f_data, i) = (1.0 / sqrt(n)) * REAL(ft_f_data, i);
        IMAG(ft_f_data, i) = (1.0 / sqrt(n)) * IMAG(ft_f_data, i);
    }

    // Compute product
    double p[2 * n];
    for (int i = 0; i < n; i++) {
        REAL(p, i) = REAL(fact_q, i) * REAL(ft_f_data, i) - IMAG(fact_q, i) * IMAG(ft_f_data, i);
        IMAG(p, i) = REAL(fact_q, i) * IMAG(ft_f_data, i) + IMAG(fact_q, i) * REAL(ft_f_data, i);
    }

    // Compute the Fourier Transform
    for (int i = 0; i < n; i++) {
        REAL(p, i) = d * sqrt(n / (2 * M_PI)) * REAL(p, i);
    }

    // Write data to file
    fprintf(ft_data, "# X_val f(x) k_val FT(f(x))\n");
    for (int i = 0; i < n; i++) {
        fprintf(ft_data, "%f %f %f %f\n", x_arr[i], f_data[i], k_arr[i], REAL(p, i));
    }
    fclose(ft_data);
    return 0;
}
