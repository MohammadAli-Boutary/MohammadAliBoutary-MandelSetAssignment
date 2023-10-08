#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N 2
#define NPixels 800

int main() {
    int maxiter = 1000; 
    double real_min = -N;
    double real_max = N;
    double imag_min = -N;
    double imag_max = N;
    int width = NPixels;
    int height = NPixels;
    double scale_real = (real_max - real_min) / (double)width;
    double scale_imag = (imag_max - imag_min) / (double)height;
    long min_color = 0;
    long max_color = 0;
    double scale_color = (double)(max_color - min_color) / (double)(maxiter - 1);
    clock_t start_time = clock();

    for (int row = 0; row < height; ++row) {
        for (int col = 0; col < width; ++col) {
            double c_real = real_min + (col * scale_real);
            double c_imag = imag_min + ((height - 1 - row) * scale_imag);
            double z_real = 0.0;
            double z_imag = 0.0;
            int k = 0;
            double lengthsq, temp;

            while (k < maxiter) {
                temp = z_real * z_real - z_imag * z_imag + c_real;
                z_imag = 2.0 * z_real * z_imag + c_imag;
                z_real = temp;
                lengthsq = z_real * z_real + z_imag * z_imag;
                if (lengthsq >= (N * N)) {
                    break;
                }
                ++k;
            }
            long color = (long)((k - 1) * scale_color) + min_color;
        }
    }
    clock_t end_time = clock();
    double execution_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("Execution Time: %f seconds\n", execution_time);
    return 0;
}
