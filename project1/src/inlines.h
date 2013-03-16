inline int calculate_cell_index(const int i, const int j, const int k, const int c_x, const int c_y, const int c_z) {
    return i*c_y*c_z + j*c_y + k;
}

inline int sign(double d) {
    return d<0?-1:d>0;
}

