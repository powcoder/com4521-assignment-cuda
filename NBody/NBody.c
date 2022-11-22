https://powcoder.com
代写代考加微信 powcoder
Assignment Project Exam Help
Add WeChat powcoder
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <omp.h>

#include "NBody.h"
#include "NBodyVisualiser.h"

#define USER_NAME "acu19jz"		//replace with your username
#define COMMENT_CHAR         '#' //comment char which will be skip when read from input
#define LINE_LEN   1000


void print_help();
void step(void);


int N; // Number of Nbody
int D; // dimension
int Iter = 1; // iteration number
nbody* input_data;
float* den_arr; //dencity array
MODE M; // operation mode, either  'CPU' or 'OPENMP'




int main(int argc, char* argv[]) {

	print_help();


	//TODO: Processes the command line arguments
	//Check input 1 and 2 exceptions
	if (!atoi(argv[1]) || !atoi(argv[2])) {
		fprintf(stderr, "Invalid input\n");
		exit(1);
	}

	//Check input 3 exception
	if (strcmp(argv[3], "CPU") != 0 && strcmp(argv[3], "OPENMP") != 0) {
		fprintf(stderr, "Invalid input\n");
		exit(1);
	}

	//Get the value of N
	N = atoi(argv[1]);
	input_data = (nbody*)malloc(N * sizeof(nbody));//Allocate memmory to the array to store input

	//Get the value of D
	D = atoi(argv[2]);
	den_arr = (float*)malloc(D * D * sizeof(float));//Allocate memmory to the array to store density

	//Select run model for step()


	if (strcmp(argv[3], "CPU") == 0) M = 0;
	if (strcmp(argv[3], "OPENMP") == 0) M = 1;


	//If there are no more inputs, run visualisation model
	if (argc == 4)
	{
		//If there is no input file, get random input 
		for (int i = 0; i < N; i++) {
			input_data[i].m = 1.0 / N;
			input_data[i].x = (float)rand() / (float)RAND_MAX;
			input_data[i].y = (float)rand() / (float)RAND_MAX;
			input_data[i].vx = 0;
			input_data[i].vy = 0;
		}

		//Call visualisation model
		printf("Visualisation model\n");
		initViewer(N, D, M, step);
		setNBodyPositions(input_data);
		setHistogramData(den_arr);
		startVisualisationLoop();
	}


	//Check if there are more more inputs
	if (argc > 4)
	{
		//Exceptions
		if (strcmp(argv[4], "-i") != 0 && strcmp(argv[4], "-f") != 0)
		{
			fprintf(stderr, "Ivalid input, should be -f or -i");
			exit(1);
		}
		else
		{
			//Check -i or -f
			for (int i = 4; i < argc; i++)
			{
				if (strcmp(argv[i], "-i") == 0)
				{
					if (!atoi(argv[i + 1])) {
						fprintf(stderr, "Invalid input, num of iterations should be an int\n");
						exit(1);
					}
					else
					{
						Iter = atoi(argv[i + 1]); // Iteration

						//Call step to run console model and calculate time
						clock_t begin = clock();
						step();
						clock_t end = clock();
						double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
						printf("Console model, time spent: %lf\n", time_spent);

					}

				}
				if (strcmp(argv[i], "-f") == 0)
				{
					//Read input file, should be a txt file
					read_input(argv[i + 1]);

					//Call console model
					clock_t begin = clock();
					step();
					clock_t end = clock();
					double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
					printf("Console model, time spent: %lf \n", time_spent);

					//Call visualisation model
					printf("Visualisation model");
					initViewer(N, D, M, step);
					setNBodyPositions(input_data);
					setHistogramData(den_arr);
					startVisualisationLoop();



				}
			}
		}

	}



	//Free memmory
	free(input_data);
	free(den_arr);

	return 0;
}




void step(void)
{
	float* ax_arr = NULL; // Array to store a
	float* ay_arr = NULL;

	//Initialize density array
	for (int i = 0; i < D * D; i++)
	{
		den_arr[i] = 0;
	}

	//Allocate memory to acceleration arrary
	ax_arr = (float*)malloc(N * sizeof(float));
	ay_arr = (float*)malloc(N * sizeof(float));


	int i;


#pragma omp parallel for schedule(dynamic) if(M==1)
	for (i = 0; i < Iter; i++)
	{
		int k;
		//#pragma omp parallel for if(M==1)
		for (k = 0; k < N; k++)
		{

			float x_total = 0;
			float y_total = 0; // Use for calculate a

			for (int j = 0; j < N; j++)
			{
				//printf("j %d\n", j);
				if (k == j) continue;
				float x_below, y_below;
				float x_up, y_up;
				x_below = pow(pow(input_data[j].x - input_data[k].x, 2) + pow(SOFTENING, 2), 3.0 / 2);

				x_up = input_data[j].m * (input_data[j].x - input_data[k].x);

				x_total += x_up / x_below;

				y_below = pow(pow(input_data[j].y - input_data[k].y, 2) + pow(SOFTENING, 2), 3.0 / 2);
				y_up = input_data[j].m * (input_data[j].y - input_data[k].y);
				y_total += y_up / y_below;

			}
			ax_arr[k] = G * x_total;
			ay_arr[k] = G * y_total;

			input_data[k].vx += ax_arr[k] * dt;  //velocity
			input_data[k].vy += ay_arr[k] * dt;

			input_data[k].x += input_data[k].vx * dt; // location
			input_data[k].y += input_data[k].vy * dt;

			if (input_data[k].x < 1 && input_data[k].y < 1 && input_data[k].x>0 && input_data[k].y>0)
			{

				int x_index = floor(input_data[k].x / (1.0 / D));
				int y_index = floor(input_data[k].y / (1.0 / D));
#pragma omp atomic
				den_arr[x_index * D + y_index] += (float)1.0 / N;

			}
		}
	}
	free(ax_arr);
	free(ay_arr);

}

int read_input(char* filename) {

	char line[LINE_LEN];
	FILE* f;
	char* token;
	char* pr;//  pointer to turn string to float using strod
	float num;


	if ((f = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "Can't open file %c", filename);
		exit(1);
	}

	int counter = 0;
	rewind(f); // set to the beginning of the f
	while (fgets(line, LINE_LEN, f) != NULL)
	{

		if (line[0] != COMMENT_CHAR)
		{

			token = strtok(line, ",");
			//
			int k = 0;
			float arr[5];
			int i = 0;
			while (token != NULL)
			{
				num = strtod(token, &pr); arr[i] = num;
				i++;
				token = strtok(NULL, ",:");
			}
			if (arr[4] == 0)
			{
				arr[4] = 1.0 / N;
			}
			input_data[counter].x = arr[0];
			input_data[counter].y = arr[1];
			input_data[counter].vx = arr[2];
			input_data[counter].vy = arr[3];
			input_data[counter].m = arr[4];
			counter++;
		}

	}
	if (counter != N)
	{
		fprintf(stderr, "Input file N does not match N\n");
		exit(1);
	}

}

void print_help() {
	printf("nbody_%s N D M [-i I] [-i input_file]\n", USER_NAME);

	printf("where:\n");
	printf("\tN                Is the number of bodies to simulate.\n");
	printf("\tD                Is the integer dimension of the activity grid. The Grid has D*D locations.\n");
	printf("\tM                Is the operation mode, either  'CPU' or 'OPENMP'\n");
	printf("\t[-i I]           Optionally specifies the number of simulation iterations 'I' to perform. Specifying no value will use visualisation mode. \n");
	printf("\t[-f input_file]  Optionally specifies an input file with an initial N bodies of data. If not specified random data will be created.\n");
}
