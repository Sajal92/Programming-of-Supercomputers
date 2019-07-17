#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <errno.h>

int main(int argc, char** argv) {

	int i, j;
	int rows, columns;
	char matrix_name[200];
	char vector_name[200];
	FILE *matrix_file;
	FILE *vector_file;

	if( argc == 3 ) { printf("Creating a matrix of size: %s\n", argv[1]); }
	else if( argc > 3 ) { printf("Only the size of the square matrix should be supplied.\n"); return 0; }
	else { printf("The size of the square matrix and its name must be suplied (eg. ./mgen 16 size16x16 ).\n"); return 0; }

	srand(time(0));

	const char *p = argv[1];
	char *end;
	rows = columns = strtol(p, &end, 10);
	if (errno == ERANGE){
		perror("range error, got ");
		exit(-1);
	}

	sprintf(matrix_name, "%s.mat", argv[2]);
	sprintf(vector_name, "%s.vec", argv[2]);
	printf("Matrix file name: \"%s\"; vector file name: \"%s\";\n", matrix_name, vector_name);

	matrix_file = fopen(matrix_name, "w+");
	fprintf(matrix_file, "%d %d\n", rows, columns);
	for(i=0; i<rows; i++) {
		for(j=0; j<columns; j++) {
			if(j == i){
				fprintf(matrix_file,"%d ", (rand()%100 + 1));
			} else {
				fprintf(matrix_file,"%d ", rand()%10);
			}
		}
		fprintf(matrix_file,"\n");
	}
	fclose(matrix_file);

	vector_file = fopen(vector_name, "w+");
	fprintf(vector_file, "%d\n", columns);
	for(j=0; j<columns; j++) {
		fprintf(vector_file,"%d ", (rand()%100 + 1));
	}
	fclose(vector_file);

	printf("Matrix and vector files written.\n");

	return 0;
}
