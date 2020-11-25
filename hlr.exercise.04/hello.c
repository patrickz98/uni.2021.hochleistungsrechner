#include <omp.h>

#include <stdio.h>
#include <stdlib.h>

int main (int argc, char** argv)
{
	int thread_id, nthreads;
	
	#pragma omp parallel private(thread_id)
	{
		thread_id = omp_get_thread_num();

		printf("Hello World from thread %d\n", thread_id);
		
		#pragma omp barrier
		
		if (thread_id == 0)
		{
			nthreads = omp_get_num_threads();
			printf("There are %d threads\n",nthreads);
		}
	}

	return EXIT_SUCCESS;
}
