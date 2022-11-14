// #include "bin.h"
#include "common.h"
#include <cstdio>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <vector>

class BinGrid{
    int x_num_bins, y_num_bins;
    double size;

    public:
    std::vector<std::vector<int> > bins;

    BinGrid(int x, int y, double size);
    void add_particles_to_bins(int num_particles, particle_t *particles);
    void clear_grid();
    int num_bins();
    void add_to_grid(particle_t *particle, int index);
};
