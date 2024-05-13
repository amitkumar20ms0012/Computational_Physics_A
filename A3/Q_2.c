#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

double f(double x) {
    if (x == 0.0) {
        return 1.0;
    } else {
        return (sin(x) / x);
    }
}

int main() {
    int n = 512; // Number of sample points
    float x_min = -50.0, x_max = 50.0, d = 0.0, *k_arr, *x_arr; // declaring x_min, x_max and delta(d)
    fftw_complex *in, *out, *ft_factors, *prod;
    FILE *ft_data;
    fftw_plan p;

    // Allocate memory for arrays
    x_arr = calloc(n, sizeof(float));
    k_arr = calloc(n, sizeof(float));
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
    ft_factors = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
    prod = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);

    // Create FFTW plan
    p = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    // Open file for writing
    ft_data = fopen("fft_data.txt", "w");

    // Compute delta (spacing between sample points)
    d = (x_max - x_min) / (n - 1);

    // Generate sample points and initialize input array
    for (int i = 0; i < n; i++) {
        x_arr[i] = x_min + i * d;
        if (i < n / 2)
            k_arr[i] = 2 * M_PI * (i / (n * d));
        else
            k_arr[i] = 2 * M_PI * ((i - n) / (n * d));

        in[i][0] = f(x_arr[i]);
        in[i][1] = 0.0;
        ft_factors[i][0] = cos(k_arr[i] * x_min);
        ft_factors[i][1] = -sin(k_arr[i] * x_min);
    }

    // Execute FFT
    fftw_execute(p);

    // Print DFT
    printf("DFT Printing\n");
    for (int i = 0; i < n; i++) {
        printf("%f\t%f\n", out[i][0], out[i][1]);
    }

    // Normalize
    for (int i = 0; i < n; i++) {
        out[i][0] = (1.0 / sqrt(n)) * out[i][0];
        out[i][1] = (1.0 / sqrt(n)) * out[i][1];
    }

    // Complex multiplication by factors
    for (int i = 0; i < n; i++) {
        prod[i][0] = ft_factors[i][0] * out[i][0] - ft_factors[i][1] * out[i][1];
        prod[i][1] = ft_factors[i][0] * out[i][1] + ft_factors[i][1] * out[i][0];
    }

    // Constructing FT
    for (int i = 0; i < n; i++) {
        out[i][0] = sqrt(n / (2 * M_PI)) * d * prod[i][0];
        out[i][1] = sqrt(n / (2 * M_PI)) * d * prod[i][1];
    }

    // Write data to file
    fprintf(ft_data, "# X_val f(x) k_val FT(f(x))\n");
    for (int i = 0; i < n; i++) {
        if (i == 0) {
            printf("Creating file\n");
        }
        fprintf(ft_data, "%f\t%f\t%f\t%f\n", x_arr[i], in[i][0], k_arr[i], out[i][0]);
    }
    fclose(ft_data);

    // Free allocated memory
    free(x_arr);
    free(k_arr);
    fftw_free(in);
    fftw_free(out);
    fftw_free(ft_factors);
    fftw_free(prod);
    fftw_destroy_plan(p);

    return 0;
}
