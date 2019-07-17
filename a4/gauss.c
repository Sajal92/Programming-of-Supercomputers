#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

int main(int argc, char** argv) {
    char matrix_name[200], vector_name[200], solution_name[200];
    int dims[2], rhs_rows;
    int rows, columns, size, rank;
    double *matrix_local_block, *rhs_local_block;
    double total_time, io_time = 0, kernel_time, mpi_time = 0;
    double total_start, io_start, kernel_start, mpi_start;
    MPI_File matrix_file, vector_file, solution_file;
    MPI_Request reqs[2];
    MPI_Status stts[2], status;
    int ret;

    if (argc != 2) { 
        perror("The base name of the input matrix and vector files must be given\n"); 
        exit(-1);
    }

    int print_a = 0;
    int print_b = 0;
    int print_x = 0;

    sprintf(matrix_name,   "%s.mat.bin", argv[1]);
    sprintf(vector_name,   "%s.vec.bin", argv[1]);
    sprintf(solution_name, "%s.sol.bin", argv[1]);

    MPI_Init(&argc, &argv); 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        printf("Solving the Ax=b system with Gaussian Elimination:\n");
        printf("READ:  Matrix   (A) file name: \"%s\"\n", matrix_name);
        printf("READ:  RHS      (b) file name: \"%s\"\n", vector_name);
        printf("WRITE: Solution (x) file name: \"%s\"\n", solution_name);
    }

    total_start = MPI_Wtime();

    int row, column, index, i;
    io_start = MPI_Wtime();
    
    char error[200];
    int errLen;
    ret = MPI_File_open(MPI_COMM_WORLD, matrix_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &matrix_file);
    if (ret != 0) {
        MPI_Error_string(ret, error, &errLen);
        printf("%s\n", error);
        perror("Could not open the specified matrix file");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    
    MPI_File_read_at_all(matrix_file, 0, dims, 2, MPI_INT, &status);
    rows = dims[0]; columns = dims[1];
    if (rows != columns) {
        perror("Only square matrices are allowed");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    if (rows % size != 0) {
        perror("The matrix should be divisible by the number of processes");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    
    ret = MPI_File_open(MPI_COMM_WORLD, vector_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &vector_file);
    if (ret != 0) {
        MPI_Error_string(ret, error, &errLen);
        printf("%s\n", error);
        perror("Could not open the specified vector file");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    MPI_File_read_at_all(vector_file, 0, &rhs_rows, 1, MPI_INT, &status);
    if (rhs_rows != rows){
        perror("RHS rows must match the sizes of A");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    
    int local_block_size = rows / size;
    
    MPI_File_set_view(matrix_file, 2 * sizeof(int) + rank * local_block_size * columns * sizeof(double), MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
    MPI_File_set_view(vector_file, sizeof(int) + rank * local_block_size * sizeof(double), MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
    matrix_local_block = (double *) malloc(local_block_size * columns * sizeof(double));
    rhs_local_block = (double *) malloc(local_block_size * sizeof(double));
    MPI_File_iread_all(matrix_file, matrix_local_block, local_block_size * columns, MPI_DOUBLE, &reqs[0]);
    MPI_File_iread_all(vector_file, rhs_local_block, local_block_size, MPI_DOUBLE, &reqs[1]);
    MPI_Waitall(2, reqs, stts);
    
    MPI_File_close(&matrix_file);
    MPI_File_close(&vector_file);
    
    io_time += MPI_Wtime() - io_start;
    
    int process, column_pivot;

	double tmp, pivot;
	double *pivots = (double *) malloc((local_block_size + (rows * local_block_size) + 1) * sizeof(double));
	double *local_work_buffer = (double *) malloc(local_block_size * sizeof(double));
	double *accumulation_buffer = (double *) malloc(local_block_size * 2 * sizeof(double));
	double *solution_local_block = (double *) malloc(local_block_size * sizeof(double));

	kernel_start = MPI_Wtime();

	for(process = 0; process < rank; process++) {
		mpi_start = MPI_Wtime();
		MPI_Recv(pivots, (local_block_size * rows + local_block_size + 1), MPI_DOUBLE, process, process, MPI_COMM_WORLD, &status);
		mpi_time += MPI_Wtime() - mpi_start;

		for(row = 0; row < local_block_size; row++){
			column_pivot = ((int) pivots[0]) * local_block_size + row;
			for (i = 0; i < local_block_size; i++){
				index = i * rows;
				tmp = matrix_local_block[index + column_pivot];
				for (column = column_pivot; column < columns; column++){
					matrix_local_block[index + column] -=  tmp * pivots[(row * rows) + (column + local_block_size + 1)];
				}
				rhs_local_block[i] -= tmp * pivots[row + 1];
				matrix_local_block[index + column_pivot] = 0.0;
			}
		}
	}

	for(row = 0; row < local_block_size; row++){
		column_pivot = (rank * local_block_size) + row;
		index = row * rows;
		pivot = matrix_local_block[index + column_pivot];
		assert(pivot!= 0);

		for (column = column_pivot; column < columns; column++){
			matrix_local_block[index + column] = matrix_local_block[index + column]/pivot; 
			pivots[index + column + local_block_size + 1] = matrix_local_block[index + column];
		}

		local_work_buffer[row] = (rhs_local_block[row])/pivot;
		pivots[row+1] =  local_work_buffer[row];

		for (i = (row + 1); i < local_block_size; i++) {
			tmp = matrix_local_block[i*rows + column_pivot];
			for (column = column_pivot+1; column < columns; column++){
				matrix_local_block[i*rows+column] -=  tmp * pivots[index + column + local_block_size + 1];
			}
			rhs_local_block[i] -= tmp * local_work_buffer[row];
			matrix_local_block[i * rows + row] = 0;
		}
	}

	for (process = (rank + 1); process < size; process++) {
		pivots[0] = (double) rank;
		mpi_start = MPI_Wtime();
		MPI_Send( pivots, (local_block_size * rows + local_block_size + 1), MPI_DOUBLE, process, rank, MPI_COMM_WORLD);
		mpi_time += MPI_Wtime() - mpi_start;
	} 

	for (process = (rank + 1); process<size; ++process) {
		mpi_start = MPI_Wtime();
		MPI_Recv( accumulation_buffer, (2 * local_block_size), MPI_DOUBLE, process, process, MPI_COMM_WORLD, &status); 
		mpi_time += MPI_Wtime() - mpi_start;

		for (row  = (local_block_size - 1); row >= 0; row--) {
			for (column = (local_block_size - 1);column >= 0; column--) {
				index = (int) accumulation_buffer[column];
				local_work_buffer[row] -= accumulation_buffer[local_block_size + column] * matrix_local_block[row * rows + index];
			}
		}
	}

	for (row = (local_block_size - 1); row >= 0; row--) {
		index = rank * local_block_size + row;
		accumulation_buffer[row] = (double) index;
		accumulation_buffer[local_block_size+row] = solution_local_block[row] = local_work_buffer[row];
		for (i = (row - 1); i >= 0; i--){
			local_work_buffer[i] -= solution_local_block[row] * matrix_local_block[ (i * rows) + index];
		}
	}

	for (process = 0; process < rank; process++){
		mpi_start = MPI_Wtime();
		MPI_Send( accumulation_buffer, (2 * local_block_size), MPI_DOUBLE, process, rank, MPI_COMM_WORLD); 
		mpi_time += MPI_Wtime() - mpi_start;
	}

	kernel_time = MPI_Wtime() - kernel_start;

    io_start = MPI_Wtime();
    ret = MPI_File_open(MPI_COMM_WORLD, solution_name, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &solution_file);
    if (ret != 0) {
        MPI_Error_string(ret, error, &errLen);
        printf("%s\n", error);
        perror("Could not open the specified solution file");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    
    if (rank == 0)
        MPI_File_write_at(solution_file, 0, &rows, 1, MPI_INT, &status);
    MPI_File_set_view(solution_file, sizeof(int) + rank * local_block_size * sizeof(double), MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
    MPI_File_write_all(solution_file, solution_local_block, local_block_size, MPI_DOUBLE, &status);
    
    MPI_File_close(&solution_file);
    io_time += MPI_Wtime() - io_start;

    /*if(print_a){
        printf("\nSystem Matrix (A):\n");
        for (row = 0; row < rows; row++) {
            for (column = 0; column < columns; column++){
                printf("%4.1f ", matrix_2d_mapped[row][column]);
            }
            printf("\n");
        }
    }

    if(print_b){
        printf("\nRHS Vector (b):\n");
        for (row = 0; row < rows; row++) {
            printf("%4.1f\n", rhs[row]);
        }
    }

    if(print_x){
        printf("\n\nSolution Vector (x):\n");
        for(row = 0; row < rows; row++){
            printf("%4.4f\n",solution[row]);
        }
    }*/

	total_time = MPI_Wtime() - total_start;

	printf("[R%02d] Times: IO: %f; Compute: %f; MPI: %f; Total: %f;\n", 
			rank, io_time, kernel_time, mpi_time, total_time);
    
    free(matrix_local_block);
	free(rhs_local_block);
	free(pivots);
	free(local_work_buffer);
	free(accumulation_buffer);
	free(solution_local_block);

	MPI_Finalize(); 
	return 0;
}

