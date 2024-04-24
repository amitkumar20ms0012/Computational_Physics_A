#include <stdio.h>
#include <math.h>

// Function to compute the exact solution
double exact_solution(double t) {
    return pow(t + 1, 2) - 0.5 * exp(t);
}

// Function to implement Euler's method and write data to file
void euler_method(double h, double y0) {
    double t = 0;
    double y = y0;
    double error_bound, error;

    FILE *fp;
    fp = fopen("Q_16_data.txt", "w");
    if (fp == NULL) {
        printf("Error opening file.\n");
        return;
    }

    fprintf(fp, "t\t\tEuler's Method\t\tExact Solution\t\tError\t\tError Bound\n");

    while (t <= 2) {
        error_bound = 0.2 * (exp(t) - 1) * fabs((pow(t + 1, 2) + 0.5 * exp(t))); // Upper bound for the error
        error = fabs(exact_solution(t) - y); // Error between exact and numerical solutions
        fprintf(fp, "%.4f\t\t%.6f\t\t%.6f\t\t%.6f\t\t%.6f\n", t, y, exact_solution(t), error, error_bound);

        // Euler's method iteration
        y += h * (y - pow(t, 2) + 1);
        t += h;
    }

    fclose(fp);
}

int main() {
    double h = 0.2;
    double y0 = 0.5;

    // Call Euler's method function
    euler_method(h, y0);

    return 0;
}
