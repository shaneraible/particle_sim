
#include <vector>


class Bin{
    std::vector<int> particle_indeces;
    
    public:
    int get_size();
    int particle_at(int i);
    void add_particle(int idx);
};