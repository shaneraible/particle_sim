#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <cuda.h>
#include "common.h"

#define NUM_THREADS 256
#define MAX_BIN_SIZE 15
extern double size;

typedef struct {
    particle_t particles[MAX_BIN_SIZE];
    int size;
} bin_t;

//
//  benchmarking program
//

__device__ void apply_force_gpu(particle_t &particle, particle_t &neighbor)
{
  double dx = neighbor.x - particle.x;
  double dy = neighbor.y - particle.y;
  double r2 = dx * dx + dy * dy;
  if( r2 > cutoff*cutoff )
      return;
  //r2 = fmax( r2, min_r*min_r );
  r2 = (r2 > min_r*min_r) ? r2 : min_r*min_r;
  double r = sqrt( r2 );

  //
  //  very simple short-range repulsive force
  //
  double coef = ( 1 - cutoff / r ) / r2 / mass;
  particle.ax += coef * dx;
  particle.ay += coef * dy;

}

__device__ void rebin_gpu(particle_t * particles, bin_t * bins, int n, int bins_per_dim){

}

__global__ void compute_forces_gpu(particle_t * particles, bin_t * bins, int n, int bins_per_dim, int size)
{
    // Get thread (particle) ID
    int thread = threadIdx.x + blockIdx.x * blockDim.x;
    if(thread >= n) return;

    int i = thread / bins_per_dim;
    int j = thread % bins_per_dim;

    int block_num = i*bins_per_dim+j;
    bin_t &current_bin = bins[block_num];

    for (int k = 0; k < current_bin.size; k++){
        current_bin.particles[k].ax = current_bin.particles[k].ay = 0;
    }

    
    // zero out acceleration
    // current_bin.particles[particle_idx].ax = current_bin.particles[particle_idx].ay = 0;
    // for (int i = 0; i < current_bin.size; i++){
    //     int bin_r = round(double(particles[i].y)/size*(bins_per_dim-1));
    //     int bin_c = round(double(particles[i].x)/size*(bins_per_dim-1));

    //     for(int r = max(bin_r - 1, 0); r <= min(bin_r+1, bins_per_dim - 1); r ++)
    //     {
    //         for(int c = max(bin_c - 1, 0); c <= min(bin_c+1, bins_per_dim - 1); c++)
    //         {
    //             bin_t &neighbor = bins[r*bins_per_dim + c];
    //             //forces within this bin
    //             for (int j = 0; j < neighbor.size; j++){
    //                 if(threadIdx.x!=j) apply_force_gpu(current_bin.particles[i], neighbor.particles[j]);
    //             }
    //         }
    //     }
    // } 

    // //forces within this bin
    for (int i = 0; i < current_bin.size; i++){
      for (int j = 0; j < current_bin.size; j++){
        if(threadIdx.x!=j) apply_force_gpu(current_bin.particles[i], current_bin.particles[j]);
      }
    }   

    // for(int j = 0 ; j < n ; j++)
    //     apply_force_gpu(particles[tid], particles[j]);

}

__global__ void move_gpu (particle_t * particles, int n, double size)
{

    // Get thread (particle) ID
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if(tid >= n) return;

    particle_t * p = &particles[tid];
    //
    //  slightly simplified Velocity Verlet integration
    //  conserves energy better than explicit Euler method
    //
    p->vx += p->ax * dt;
    p->vy += p->ay * dt;
    p->x  += p->vx * dt;
    p->y  += p->vy * dt;

    //
    //  bounce from walls
    //
    while( p->x < 0 || p->x > size )
    {
        p->x  = p->x < 0 ? -(p->x) : 2*size-p->x;
        p->vx = -(p->vx);
    }
    while( p->y < 0 || p->y > size )
    {
        p->y  = p->y < 0 ? -(p->y) : 2*size-p->y;
        p->vy = -(p->vy);
    }

}

void init_bins(bin_t *bins, particle_t* particles, int n, int bins_per_dim, double grid_size){
    for(int i=0; i<bins_per_dim*bins_per_dim; i++){

        bins[i].size = 0;
    }

        //store the point into the grid
    for (int i = 0; i < n; i++){
        int x = min((int)(particles[i].x / (grid_size/bins_per_dim)), bins_per_dim - 1);
        int y = min((int)(particles[i].y / (grid_size/bins_per_dim)), bins_per_dim - 1);
        int index = y*bins_per_dim + x;
    
        if(bins[index].size >= MAX_BIN_SIZE){
            printf("BIN SIZE TOO SMALL\n");
            return;
        }
        bins[index].particles[bins[index].size++] = particles[i];

    }
}

int main( int argc, char **argv )
{    

    // This takes a few seconds to initialize the runtime
    cudaThreadSynchronize(); 

    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        return 0;
    }
    
    int n = read_int( argc, argv, "-n", 1000 );

    char *savename = read_string( argc, argv, "-o", NULL );
    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );

    set_size( n );
    init_particles( n, particles );

    double size = get_size();
    int bins_per_dim = int(sqrt(n/4));

    bin_t *bins = (bin_t*) malloc(bins_per_dim*bins_per_dim*sizeof(bin_t));
    init_bins(bins, particles, n, bins_per_dim, size);
    
    bin_t *bins_gpu; 
    cudaMalloc((void **) &bins_gpu, bins_per_dim * bins_per_dim * sizeof(bin_t));

    // GPU particle data structure
    particle_t * d_particles;
    cudaMalloc((void **) &d_particles, n * sizeof(particle_t));


    cudaThreadSynchronize();
    double copy_time = read_timer( );

    // Copy the particles to the GPU
    cudaMemcpy(d_particles, particles, n * sizeof(particle_t), cudaMemcpyHostToDevice);
    cudaMemcpy(bins_gpu, bins, bins_per_dim * bins_per_dim * sizeof(bin_t), cudaMemcpyHostToDevice);

    cudaThreadSynchronize();
    copy_time = read_timer( ) - copy_time;
    
    //
    //  simulate a number of time steps
    //
    cudaThreadSynchronize();
    double simulation_time = read_timer( );

    for( int step = 0; step < NSTEPS; step++ )
    {
        //
        //  compute forces
        //

        int blks = (bins_per_dim*bins_per_dim + NUM_THREADS - 1)/NUM_THREADS;
        compute_forces_gpu <<< blks, NUM_THREADS >>> (d_particles, bins_gpu, n, bins_per_dim, size);
        
        //
        //  move particles
        //
	    move_gpu <<< blks, NUM_THREADS >>> (d_particles, n, size);
        
        //
        //  save if necessary
        //
        if( fsave && (step%SAVEFREQ) == 0 ) {
	        // Copy the particles back to the CPU
            cudaMemcpy(particles, d_particles, n * sizeof(particle_t), cudaMemcpyDeviceToHost);
            
            cudaMemcpy(bins, bins_gpu, bins_per_dim*bins_per_dim*sizeof(bin_t), cudaMemcpyDeviceToHost);
            int count = 0;
            for (int p = 0; p < bins_per_dim; p++){
                for (int q = 0; q < bins_per_dim; q++){
                    bin_t& curr = bins[p*bins_per_dim + q];
                    for (int k = 0; k < curr.size; k++){
                        particles[count++] = curr.particles[k];
                    }
                }
            }

            save( fsave, n, particles);
            cudaMemset(bins_gpu, 0, bins_per_dim * bins_per_dim * sizeof(bin_t));
            init_bins(bins, particles, n, bins_per_dim, size);
            cudaMemcpy(bins_gpu, bins, bins_per_dim * bins_per_dim * sizeof(bin_t), cudaMemcpyHostToDevice);

	    }
    }
    cudaThreadSynchronize();
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "CPU-GPU copy time = %g seconds\n", copy_time);
    printf( "n = %d, simulation time = %g seconds\n", n, simulation_time );
    
    free( particles );
    cudaFree(d_particles);
    if( fsave )
        fclose( fsave );
    
    return 0;
}
