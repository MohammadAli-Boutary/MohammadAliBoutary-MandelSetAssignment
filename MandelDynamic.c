
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <gd.h>
#include <unistd.h>
#include <float.h>

#define N 2
#define NPixels 800

#define WORK_TAG 1
#define DATA_TAG 2
#define STOP_TAG 3

typedef struct {
    double real, imag;
} complex;

int main(int argc, char *argv[]) {
    int num_procs; // number of processors
    int rank;      // current process id
    int return_status; // return status of master and slave work
    int maxiter = 1000; // number of iterations to calculate if the pixel goes out of bounds or not
    double real_min = -N;
    double real_max = N;
    double imag_min = -N;
    double imag_max = N;
    int width = NPixels;
    int height = NPixels;

    if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
        exit(EXIT_FAILURE);
    }

    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (num_procs < 2) { // number of processors must be at least 2 in order to parallelize
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    if (rank == 0) {
        long min_color = 0;
        long max_color = 256;
        int this_row, next_row;
        double start_time, end_time;
        long *data_msg = malloc((width + 1) * sizeof(*data_msg)); // the Mandelbrot set strip processors will send to the master after computing its color
        MPI_Status status;
        int tasks_not_done;
        int id;

        start_time = MPI_Wtime();
        MPI_Bcast(&min_color, 1, MPI_LONG, 0, MPI_COMM_WORLD);
        MPI_Bcast(&max_color, 1, MPI_LONG, 0, MPI_COMM_WORLD);

        next_row = 0;
        tasks_not_done = 0;

        for (int i = 1; i < num_procs; ++i) { // Start from rank 1
            MPI_Send(&next_row, 1, MPI_INT, i, WORK_TAG, MPI_COMM_WORLD);
            ++next_row;
            ++tasks_not_done;
        }

        while (tasks_not_done > 0) {
            MPI_Recv(data_msg, width + 1, MPI_LONG, MPI_ANY_SOURCE, DATA_TAG, MPI_COMM_WORLD, &status);
            --tasks_not_done;
            id = status.MPI_SOURCE;

            if (next_row < height) {
                MPI_Send(&next_row, 1, MPI_INT, id, WORK_TAG, MPI_COMM_WORLD);
                ++next_row;
                ++tasks_not_done;
            } else {
                MPI_Send(&next_row, 0, MPI_INT, id, STOP_TAG, MPI_COMM_WORLD);
            }
        }

        end_time = MPI_Wtime();
        fprintf(stdout, "Execution Time = %f\n", end_time - start_time);
    } else {
        MPI_Status status;
        int current_row;
        long min_color, max_color;
        double scale_real, scale_imag, scale_color;
        long *data_msg = malloc((width + 1) * sizeof(*data_msg));

        MPI_Bcast(&min_color, 1, MPI_LONG, 0, MPI_COMM_WORLD);
        MPI_Bcast(&max_color, 1, MPI_LONG, 0, MPI_COMM_WORLD);

        scale_real = (double) (real_max - real_min) / (double) width;
        scale_imag = (double) (imag_max - imag_min) / (double) height;

        scale_color = (double) (max_color - min_color) / (double) (maxiter - 1);

        while (((MPI_Recv(&current_row, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status)) == MPI_SUCCESS) && (status.MPI_TAG == WORK_TAG)) {
            data_msg[0] = current_row;

            for (int col = 0; col < width; ++col) {
                complex z, c;

                z.real = z.imag = 0;

                c.real = real_min + ((double) col * scale_real);
                c.imag = imag_min + ((double) (height - 1 - current_row) * scale_imag);

                int k = 0;
                double lengthsq, temp;

                do {
                    temp = z.real * z.real - z.imag * z.imag + c.real;
                    z.imag = 2 * z.real * z.imag + c.imag;
                    z.real = temp;
                    lengthsq = z.real * z.real + z.imag * z.imag;
                    ++k;
                } while (lengthsq < (N * N) && k < maxiter);

                long color = (long) ((k - 1) * scale_color) + min_color;
                
                data_msg[col + 1] = color;
            }

            MPI_Send(data_msg, width + 1, MPI_LONG, 0, DATA_TAG, MPI_COMM_WORLD);
        }
    }

    MPI_Finalize();
    return 0;
}
