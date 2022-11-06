#include "bin.h"

int Bin::particle_at(int i){
    if(i >= this->particle_indeces.size() || i<0)
    {
        return -1;
    }

    return this->particle_indeces.size();
}

int Bin::get_size(){ return this->particle_indeces.size()}

void Bin::add_particle(int idx){
    this->particle_indeces.push_back(idx);
}

void Bin::remove_particle(int idx){
    this->particle_indeces.push_back(idx);
}