#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

int isprime(int n) {
  int i, sq;

  sq = (int) sqrt(n);
  for (i=3; i<=sq; i+=2)
    if ((n%i) == 0) return 0;

  return 1;
}

int main(int argc, char *argv[])
{
  int i, p=0, max;
  //max=100000000;   // 100 million (primes=5761455)
  max=2100000000;  // 2.1 billion (primes=102886526)

  int ntasks, rank, size, length, start, h, psum;
  double start_time, end_time, wct;
  char node[MPI_MAX_PROCESSOR_NAME];
  
  MPI_Init(&argc, &argv);     //MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Rank 
  MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
  MPI_Get_processor_name(node, &length); //Name

  start_time = MPI_Wtime(); // Timer
  start = (rank*2)+1; //
  h = ntasks*2;

  MPI_Bcast(&max,1,MPI_INT,0,MPI_COMM_WORLD);  //broadcasting max variable to all ranks

  if (rank==0){
  printf("Numbers to be scanned = %d\n", max);
    for (i=start; i<=max; i+=h){
      if (isprime(i)) p++;
    }
    MPI_Reduce(&p,&psum,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    printf("Rank = %d/%d , Node = %s , Primes found = %d\n", rank, ntasks, node, p);
  }
  
  if(rank>0){
  for (i=start; i<=max; i+=h) {
    if (isprime(i)) p++;
  }
  printf("Rank = %d/%d , Node = %s , Primes found = %d\n", rank, ntasks, node, p);
  MPI_Reduce(&p,&psum,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if(rank==0){
  printf("Total primes = %d\n", psum);
  end_time = MPI_Wtime();
  wct=end_time-start_time;
  printf("Wallclock time elapsed: %.2lf seconds\n", wct);
  }

  MPI_Finalize();
}
