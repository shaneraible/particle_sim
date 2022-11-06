#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include "bin_grid.h"

//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    int navg,nabsavg=0;
    double davg,dmin, absmin=1.0, absavg=0.0;

    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-no turns off all correctness checks and particle output\n");
        return 0;
    }
    
    int n = read_int( argc, argv, "-n", 1000 );

    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );
    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );
    double s = get_size();
    int sn = int(sqrt(n/4));
    // printf("%d\n", sn);
    // BinGrid bingrid(sn, sn, s);
    // bingrid.add_particles_to_bins(n, particles);

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
	
    for( int step = 0; step < NSTEPS; step++ )
    {
	    navg = 0;
        davg = 0.0;
	    dmin = 1.0;

        BinGrid bingrid(sn, sn, s);
        bingrid.add_particles_to_bins(n, particles);

        //
        //  compute forces
        //
        // for( int i = 0; i < n; i++ )
        // {
        //     particles[i].ax = particles[i].ay = 0;
        //     for (int j = 0; j < n; j++ )
		// 		apply_force( particles[i], particles[j],&dmin,&davg,&navg);
        // }

        // for(int i = 0; i < n; i++)
        // {
        //     particles[i].ax = particles[i].ay = 0;
        //     int bin_r = particles[i].y / bin_y, bin_c = particles[i].x / bin_x;
        //     // Traversing the neighbors
        //     for(int r = max(bin_r - 1, 0); r <= min(bin_r+1, bin_j - 1); r ++)
        //     {
        //         for(int c = max(bin_c - 1, 0); c <= min(bin_c+1, bin_i - 1); c++)
        //         {
        //             bin_t neighbor = bin_list[r + c*bin_j];
        //             //printf("Neighbor index = %d with size: %d\n", r+c*bin_j, neighbor.bin_size);
        //             for(int j = 0; j < neighbor.bin_size; j ++)
        //                 apply_force(particles[i], particles[neighbor.indeces[j]], &dmin, &davg, &navg);    
        //         }
        //     }
        // }
        int min_r=-1, min_c=-1, min_n=-1, min_p=-1;
        int sum = (n-1)*(n)/2;

        // int i = 0;
        // for(auto bin: bingrid.bins){
        //     for(auto idx: bin){
        //         sum-=idx;
        //         // printf("%d\n", idx);

        //     }
        //     ++i;
        //     // count += bin.size();
        //     if(bin.size() == 0) printf("HERE %d\n", i);
        // }

        // if(sum!=0) printf("%d\n", sum);

        for(int i=0; i<n; i++){
            particles[i].ax = particles[i].ay = 0;
            int bin_r = round(double(particles[i].y)/s*(sn-1));
            int bin_c = round(double(particles[i].x)/s*(sn-1));

            for(int r = max(bin_r - 1, 0); r <= min(bin_r+1, sn - 1); r ++)
            {
                for(int c = max(bin_c - 1, 0); c <= min(bin_c+1, sn - 1); c++)
                {
                    auto neighbor = bingrid.bins[r*sn + c];
                    // printf("Neighbor index = %d with size: %d\n", r+c*sn, neighbor.size());
                    for(int j = 0; j < neighbor.size(); j ++){
                        int k = dmin;
                        apply_force(particles[i], particles[neighbor[j]], &dmin, &davg, &navg);
                        if(dmin<k){
                            min_r = r;
                            min_c = c;
                            min_n = neighbor[j];
                            min_p = i;
                        }
                    }
                }
            }
        }

        // if(dmin<.4)
        //     printf("step: %d -- r: %d, c: %d, dmin: %f, particle: %d, neighbor: %d\n", step, min_r, min_c, dmin, min_p, min_n);

        

        // for(int b=0; b<bingrid.num_bins(); b++){
            // auto bin = *bingrid.bins[b];
            // printf("%ld\n", bin.size());
            // for(int x=0; x<sn; x++){
            //     for(int y = 0; y<sn; y++){
            //         int idx = y*sn + x;

            //         for(int x_n=max(0, x-1); x_n<=min(sn-1, x+1); x_n++){
            //         for(int y_n =max(0, y-1); y_n<=min(sn-1, y+1); y_n++){
            //             int idx_n = y*sn + x;
                        
            //             for( int i = 0; i < bingrid.bins[idx].size(); i++ )
            //             {
            //                 particles[bingrid.bins[idx][i]].ax = particles[bingrid.bins[idx][i]].ay = 0;

            //                 for (int j = 0; j < bingrid.bins[idx_n].size(); j++ )
            //                     if(i!=j) apply_force(particles[bingrid.bins[idx][i]], particles[bingrid.bins[idx_n][j]],&dmin,&davg,&navg);
            //             }
                        
            //         }
            //         }

            //         //apply to the neighbors
            //     }
            // }
        // }

        //
        //  move particles
        //
        for( int i = 0; i < n; i++ ) 
            move( particles[i] );		

        if( find_option( argc, argv, "-no" ) == -1 )
        {
          //
          // Computing statistical data
          //
          if (navg) {
            absavg +=  davg/navg;
            nabsavg++;
          }
          if (dmin < absmin) absmin = dmin;
		
          //
          //  save if necessary
          //
          if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
        }
    }
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d, simulation time = %g seconds", n, simulation_time);

    if( find_option( argc, argv, "-no" ) == -1 )
    {
      if (nabsavg) absavg /= nabsavg;
    // 
    //  -The minimum distance absmin between 2 particles during the run of the simulation
    //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
    //  -A simulation where particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
    //
    //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
    //
    printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
    if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
    if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");     

    //
    // Printing summary data
    //
    if( fsum) 
        fprintf(fsum,"%d %g\n",n,simulation_time);
 
    //
    // Clearing space
    //
    if( fsum )
        fclose( fsum );    
    free( particles );
    if( fsave )
        fclose( fsave );
    
    return 0;
}
