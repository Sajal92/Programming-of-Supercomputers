#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

int main(int argc, char** argv) {
	char matrix_name[200], vector_name[200], solution_name[200];
	int rows, columns, size, rank;
	double **matrix_2d_mapped, *matrix_1D_mapped, *rhs, *solution;
	double total_time, io_time = 0, setup_time, kernel_time, mpi_time = 0;
	double total_start, io_start, setup_start, kernel_start, mpi_start;
	FILE *matrix_file, *vector_file, *solution_file;
	MPI_Status status;     

	if (argc != 2) { 
		perror("The base name of the input matrix and vector files must be given\n"); 
		exit(-1);
	}

	int print_a = 0;
	int print_b = 0;
	int print_x = 0;

	sprintf(matrix_name,   "%s.mat", argv[1]);
	sprintf(vector_name,   "%s.vec", argv[1]);
	sprintf(solution_name, "%s.sol", argv[1]);

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

	int row, column, index;
	if (rank == 0) {
		io_start = MPI_Wtime();
		if ((matrix_file = fopen(matrix_name, "r")) == NULL) {
			perror("Could not open the specified matrix file");
			MPI_Abort(MPI_COMM_WORLD, -1);
		}

		fscanf(matrix_file, "%d %d", &rows, &columns);     
		if (rows != columns) {
			perror("Only square matrices are allowed\n");
			MPI_Abort(MPI_COMM_WORLD, -1);
		}  	
		if (rows % size != 0) {
			perror("The matrix should be divisible by the number of processes\n");
			MPI_Abort(MPI_COMM_WORLD, -1);
		}  	

		matrix_2d_mapped = (double **) malloc(rows * sizeof(double *));
		for (row = 0; row < rows; row++) {
			matrix_2d_mapped[row] = (double *) malloc(rows * sizeof(double));
			for (column = 0; column < columns; column++) {
				fscanf(matrix_file, "%lf", &matrix_2d_mapped[row][column]);
			}
		}
		fclose(matrix_file);

		if ((vector_file = fopen(vector_name, "r")) == NULL) {
			perror("Could not open the specified vector file");
			MPI_Abort(MPI_COMM_WORLD, -1);
		}

		int rhs_rows;
		fscanf(vector_file, "%d", &rhs_rows);     
		if (rhs_rows != rows) {
			perror("RHS rows must match the sizes of A");
			MPI_Abort(MPI_COMM_WORLD, -1);
		}

		rhs = (double *) malloc(rows * sizeof(double));
		for (row = 0; row < rows; row++) {
			fscanf(vector_file, "%lf", &rhs[row]);
		}
		fclose(vector_file); 
		io_time += MPI_Wtime() - io_start;

		matrix_1D_mapped = (double *) malloc(rows * rows * sizeof(double));
		index = 0;
		for (row = 0; row < rows; row++) {
			for (column = 0; column < columns; column++) {
				matrix_1D_mapped[index++] = matrix_2d_mapped[row][column];
			}
		}
		solution = (double *) malloc(rows * sizeof(double));
	}

	setup_start = MPI_Wtime();
    
    int i;
    MPI_Group group_slaves, group_all;
    MPI_Group *groups;
    MPI_Alloc_mem(size * sizeof(MPI_Group), MPI_INFO_NULL, &groups);
    MPI_Comm_group(MPI_COMM_WORLD, &group_all);
    for (i = 0; i < size; i++) {
        const int i_rank[1] = {i};
        MPI_Group_incl(group_all, 1, i_rank, &groups[i]);
    }
    MPI_Group_difference(group_all, groups[0], &group_slaves);
    
    MPI_Win win;
	int *shared_mem;
    MPI_Win_allocate(2 * sizeof(int), sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &shared_mem, &win);
    //MPI_Alloc_mem(2 * sizeof(int), MPI_INFO_NULL, &shared_mem);
    //MPI_Win_create(shared_mem, 2 * sizeof(int), sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
	if (rank == 0) {
        shared_mem[0] = rows; shared_mem[1] = columns;
        
        MPI_Win_post(group_slaves, 0, win);
        MPI_Win_wait(win);
        //MPI_Win_fence(0, win);
	} else {
        int matrix_dim[2];
        
        MPI_Win_start(groups[0], 0, win);
        MPI_Get(matrix_dim, 2, MPI_INT, 0, 0, 2, MPI_INT, win);
        MPI_Win_complete(win);
        rows = matrix_dim[0]; columns = matrix_dim[1];
	}
	MPI_Win_free(&win);

	int local_block_size = rows / size;
	int process, column_pivot;

	double tmp, pivot;
	double *matrix_local_block = (double *) malloc(local_block_size * rows * sizeof(double));
	double *rhs_local_block = (double *) malloc(local_block_size * sizeof(double));
	double *pivots = (double *) malloc((local_block_size + rows * local_block_size + 1) * sizeof(double));
	double *local_work_buffer = (double *) malloc(local_block_size * sizeof(double));
	double *accumulation_buffer = (double *) malloc(local_block_size * 2 * sizeof(double));
	double *solution_local_block = (double *) malloc(local_block_size * sizeof(double));

    MPI_Win win_mtx, win_rhs;
    double *shared_mem_mtx, *shared_mem_rhs;
    MPI_Win_allocate(rows * rows * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &shared_mem_mtx, &win_mtx);
    MPI_Win_allocate(rows * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &shared_mem_rhs, &win_rhs);
	if (rank == 0) {
        memcpy(shared_mem_mtx, matrix_1D_mapped, rows * rows * sizeof(double));
        MPI_Win_post(group_slaves, 0, win_mtx);
        MPI_Win_wait(win_mtx);
        
        memcpy(shared_mem_rhs, rhs, rows * sizeof(double));
        MPI_Win_post(group_slaves, 0, win_rhs);
        MPI_Win_wait(win_rhs);
		
		for (i = 0; i < local_block_size * rows; i++) {
			matrix_local_block[i] = matrix_1D_mapped[i];
		}
		for (i = 0; i < local_block_size; i++) {
			rhs_local_block[i] = rhs[i];
		}
	} else {
        MPI_Win_start(groups[0], 0, win_mtx);
        MPI_Get(matrix_local_block, local_block_size * rows, MPI_DOUBLE, 0, rank * local_block_size * rows, local_block_size * rows, MPI_DOUBLE, win_mtx);
        MPI_Win_complete(win_mtx);
        
        MPI_Win_start(groups[0], 0, win_rhs);
        MPI_Get(rhs_local_block, local_block_size, MPI_DOUBLE, 0, rank * local_block_size, local_block_size, MPI_DOUBLE, win_rhs);
        MPI_Win_complete(win_rhs);
	}
	MPI_Win_free(&win_mtx);
    MPI_Win_free(&win_rhs);

	setup_time = MPI_Wtime() - setup_start;
	kernel_start = MPI_Wtime();

    MPI_Win win_pivot;
    double *shared_mem_pivot;
    MPI_Win_allocate((local_block_size * rows + local_block_size + 1) * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &shared_mem_pivot, &win_pivot);
	for (process = 0; process < rank; process++) {
		mpi_start = MPI_Wtime();
        MPI_Win_start(groups[process], 0, win_pivot);
        MPI_Get(pivots, local_block_size * rows + local_block_size + 1, MPI_DOUBLE, process, 0, local_block_size * rows + local_block_size + 1, MPI_DOUBLE, win_pivot);
        MPI_Win_complete(win_pivot);
		mpi_time += MPI_Wtime() - mpi_start;

		for (row = 0; row < local_block_size; row++) {
			column_pivot = ((int) pivots[0]) * local_block_size + row;
			for (i = 0; i < local_block_size; i++) {
				index = i * rows;
				tmp = matrix_local_block[index + column_pivot];
				for (column = column_pivot; column < columns; column++) {
					matrix_local_block[index + column] -=  tmp * pivots[(row * rows) + (column + local_block_size + 1)];
				}
				rhs_local_block[i] -= tmp * pivots[row + 1];
				matrix_local_block[index + column_pivot] = 0.0;
			}
		}
	}

	for (row = 0; row < local_block_size; row++) {
		column_pivot = rank * local_block_size + row;
		index = row * rows;
		pivot = matrix_local_block[index + column_pivot];
		assert(pivot!= 0);

		for (column = column_pivot; column < columns; column++) {
			matrix_local_block[index + column] = matrix_local_block[index + column] / pivot; 
			pivots[index + column + local_block_size + 1] = matrix_local_block[index + column];
		}

		local_work_buffer[row] = rhs_local_block[row] / pivot;
		pivots[row + 1] = local_work_buffer[row];

		for (i = (row + 1); i < local_block_size; i++) {
			tmp = matrix_local_block[i*rows + column_pivot];
			for (column = column_pivot+1; column < columns; column++) {
				matrix_local_block[i*rows+column] -=  tmp * pivots[index + column + local_block_size + 1];
			}
			rhs_local_block[i] -= tmp * local_work_buffer[row];
			matrix_local_block[i * rows + row] = 0;
		}
	}

	pivots[0] = (double) rank;
    memcpy(shared_mem_pivot, pivots, (local_block_size * rows + local_block_size + 1) * sizeof(double));
    for (process = rank + 1; process < size; process++) {
		mpi_start = MPI_Wtime();
        MPI_Win_post(groups[process], 0, win_pivot);
        MPI_Win_wait(win_pivot);
		mpi_time += MPI_Wtime() - mpi_start;
	}
	MPI_Win_free(&win_pivot);

//     MPI_Win win_acc;
//     double *shared_mem_acc;
//     MPI_Win_allocate(2 * local_block_size * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &shared_mem_acc, &win_acc);
	for (process = rank + 1; process < size; process++) {
		mpi_start = MPI_Wtime();
		MPI_Recv(accumulation_buffer, 2 * local_block_size, MPI_DOUBLE, process, process, MPI_COMM_WORLD, &status);
//         MPI_Win_start(groups[process], 0, win_acc);
//         MPI_Get(accumulation_buffer, 2 * local_block_size, MPI_DOUBLE, process, 0, 2 * local_block_size, MPI_DOUBLE, win_acc);
//         MPI_Win_complete(win_acc);
		mpi_time += MPI_Wtime() - mpi_start;

		for (row = local_block_size - 1; row >= 0; row--) {
			for (column = local_block_size - 1; column >= 0; column--) {
				index = (int) accumulation_buffer[column];
				local_work_buffer[row] -= accumulation_buffer[local_block_size + column] * matrix_local_block[row * rows + index];
			}
		}
	}

	for (row = local_block_size - 1; row >= 0; row--) {
		index = rank * local_block_size + row;
		accumulation_buffer[row] = (double) index;
		accumulation_buffer[local_block_size+row] = solution_local_block[row] = local_work_buffer[row];
		for (i = row - 1; i >= 0; i--) {
			local_work_buffer[i] -= solution_local_block[row] * matrix_local_block[i * rows + index];
		}
	}

// 	memcpy(shared_mem_acc, accumulation_buffer, 2 * local_block_size * sizeof(double));
	for (process = 0; process < rank; process++) {
		mpi_start = MPI_Wtime();
		MPI_Send(accumulation_buffer, 2 * local_block_size, MPI_DOUBLE, process, rank, MPI_COMM_WORLD);
//         MPI_Win_post(groups[process], 0, win_acc);
//         MPI_Win_wait(win_acc);
		mpi_time += MPI_Wtime() - mpi_start;
	}
// 	MPI_Win_free(&win_acc);

    MPI_Win win_gather;
    double *shared_mem_gather;
    MPI_Win_allocate(local_block_size * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &shared_mem_gather, &win_gather);
	if (rank == 0) {
		for (i = 0; i < local_block_size; i++) {
			solution[i] = solution_local_block[i];
		}
		mpi_start = MPI_Wtime();
		for (i = 1; i < size; i++) {
            MPI_Win_start(groups[i], 0, win_gather);
            MPI_Get(solution + i * local_block_size, local_block_size, MPI_DOUBLE, i, 0, local_block_size, MPI_DOUBLE, win_gather);
            MPI_Win_complete(win_gather);
		}
		mpi_time += MPI_Wtime() - mpi_start;
	} else {
        memcpy(shared_mem_gather, solution_local_block, local_block_size * sizeof(double));
		mpi_start = MPI_Wtime();
        MPI_Win_post(groups[0], 0, win_gather);
        MPI_Win_wait(win_gather);
		mpi_time += MPI_Wtime() - mpi_start;
	}

	kernel_time = MPI_Wtime() - kernel_start;

	if (rank == 0) {
		io_start = MPI_Wtime();
		if ((solution_file = fopen(solution_name, "w+")) == NULL) {
			perror("Could not open the solution file");
			MPI_Abort(MPI_COMM_WORLD, -1);
		}

		fprintf(solution_file, "%d\n", rows);
		for (i = 0; i < rows; i++) {
			fprintf(solution_file, "%f ", solution[i]);
		}
		fprintf(solution_file, "\n");
		fclose(solution_file);
		io_time += MPI_Wtime() - io_start;

		if (print_a) {
			printf("\nSystem Matrix (A):\n");
			for (row = 0; row < rows; row++) {
				for (column = 0; column < columns; column++) {
					printf("%4.1f ", matrix_2d_mapped[row][column]);
				}
				printf("\n");
			}
		}

		if (print_b) {
			printf("\nRHS Vector (b):\n");
			for (row = 0; row < rows; row++) {
				printf("%4.1f\n", rhs[row]);
			}
		}

		if (print_x) {
			printf("\n\nSolution Vector (x):\n");
			for (row = 0; row < rows; row++) {
				printf("%4.4f\n", solution[row]);
			}
		}
	}

	total_time = MPI_Wtime() - total_start;

	printf("[R%02d] Times: IO: %f; Setup: %f; Compute: %f; MPI: %f; Total: %f;\n", 
			rank, io_time, setup_time, kernel_time, mpi_time, total_time);

	if (rank == 0) {
		for (i = 0; i < rows; i++) {
			free(matrix_2d_mapped[i]);
		}
		free(matrix_2d_mapped);
		free(rhs);
		free(solution);
	}
	free(matrix_local_block);
	free(rhs_local_block);
	free(pivots);
	free(local_work_buffer);
	free(accumulation_buffer);
	free(solution_local_block);

    for (i = 0; i < size; i++)
        MPI_Group_free(&groups[i]);
    MPI_Free_mem(groups);
    MPI_Group_free(&group_all);
    MPI_Group_free(&group_slaves);
    MPI_Finalize(); 
	return 0;
}

