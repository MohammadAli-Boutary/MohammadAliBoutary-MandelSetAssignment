#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define WIDTH 800
#define HEIGHT 800
#define MAX_ITER 1000

int mandelbrot(double real, double imag) {
    int n;
    double r = 0.0;
    double i = 0.0;
    
    for (n = 0; n < MAX_ITER; n++) {
        double r2 = r * r;
        double i2 = i * i;
        if (r2 + i2 > 4.0) {
            return n;
        }
        i = 2 * r * i + imag;
        r = r2 - i2 + real;
    }
    
    return MAX_ITER;
}

int main(int argc, char** argv) {
    int rank, size;
    double start_time, end_time;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    start_time = MPI_Wtime(); 
    
    // Calculate the height chunk for each process
    int chunk_height = HEIGHT / size;
    int start_row = rank * chunk_height;
    int end_row = start_row + chunk_height;
    int image[WIDTH][chunk_height];
    
    // Calculate the Mandelbrot set for the assigned rows
    for (int y = start_row; y < end_row; y++) {
        for (int x = 0; x < WIDTH; x++) {
            double real = (x - WIDTH / 2.0) * 4.0 / WIDTH;
            double imag = (y - HEIGHT / 2.0) * 4.0 / HEIGHT;
            
            image[x][y - start_row] = mandelbrot(real, imag);
        }
    }
    
    // Gather the data from all processes
    int* gathered_image = NULL;
    if (rank == 0) {
        gathered_image = (int*)malloc(sizeof(int) * WIDTH * HEIGHT);
    }
    MPI_Gather(image, WIDTH * chunk_height, MPI_INT, gathered_image, WIDTH * chunk_height, MPI_INT, 0, MPI_COMM_WORLD);

    
    end_time = MPI_Wtime(); // End measuring time
    
    if (rank == 0) {
        printf("Execution time: %lf seconds\n", end_time - start_time);
    }
    
    MPI_Finalize();
    return 0;
}
