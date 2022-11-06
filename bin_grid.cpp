#include "bin_grid.h"


BinGrid::BinGrid(int x, int y, double size): bins(x*y, std::vector<int>()){
    this->x_num_bins = x;
    this->y_num_bins = y;
    this->size = size;
    // printf("%d x %d\n", x, y); 

}

void BinGrid::add_particles_to_bins(int num_particles, particle_t *particles){
    for (int i=0; i<num_particles; i++){
        add_to_grid(particles+i, i);
    }
}

int BinGrid::num_bins(){
    return this->x_num_bins*this->y_num_bins;
}

void BinGrid::clear_grid(){

}

void BinGrid::add_to_grid(particle_t *particle, int index_in_list){
    int idx_y = round(double(particle->y)/this->size*(this->y_num_bins-1));
    int idx_x = round(double(particle->x)/this->size*(this->x_num_bins-1));

    int single_dim_idx = idx_y*x_num_bins + idx_x;
    // printf("%f, %f -> %d, %d = %d\n", particle->x, particle->y, idx_x, idx_y, single_dim_idx); 
    // printf("%x\n", &this->bins[single_dim_idx]);

    this->bins[single_dim_idx].push_back(index_in_list);

}
