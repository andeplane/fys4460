#include "threadcontrol.h"
#include <Cell.h>
#include <inlines.h>

ThreadControl::ThreadControl()
{

}

void ThreadControl::setup(int nodes_, int cells_x, int cells_y, int cells_z, vector<Cell*> &cells) {
    for(int i=0;i<nodes_;i++) {
        ThreadNode node;
        nodes.push_back(node);
    }

    //    cell_indices_per_thread
    int cell_rows_per_node = cells_z/(nodes_-1);
    int current_node = 1;
    int cell_rows_so_far = 0;
    int cell_index;
    int max_node_index = nodes_ - 1;

    for(int k=0;k<cells_z;k++) {
        cell_rows_so_far++;
        for(int i=0;i<cells_x;i++) {
            for(int j=0;j<cells_y;j++) {
                cell_index = calculate_cell_index(i,j,k,cells_x,cells_y,cells_z);
                nodes[current_node].owned_cells.insert(cell_index);
            }
        }

        if(cell_rows_so_far >= cell_rows_per_node && current_node < max_node_index) {
            current_node++;
            cell_rows_so_far = 0;
        }
    }

    for(int i=0;i<nodes_;i++) {
        for(set<int>::iterator it=nodes[i].owned_cells.begin(); it!= nodes[i].owned_cells.end();it++) {
            int cell_index = *it;
            Cell *cell = cells[cell_index];
            Cell *neighbour_above = cells[calculate_cell_index(cell->i,cell->j,(cell->k+1)%cells_z,cells_x,cells_y,cells_z)];

            nodes[i].connected_cells.insert(neighbour_above->index);
            nodes[i].connected_cells.insert(cell_index);
        }
    }
}
