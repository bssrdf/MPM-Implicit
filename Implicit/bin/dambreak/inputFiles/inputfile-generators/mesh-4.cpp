
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

int main () {

    double xmin = 0.0, xmax = 0.4, xspacing = 0.05;
    double ymin = 0.0, ymax = 0.4, yspacing = 0.05;

    int xnumNodes = (xmax - xmin)/xspacing + 1;
    int ynumNodes = (ymax - ymin)/yspacing + 1;
    int xnumElem = (xmax - xmin) / (xspacing);
    int ynumElem = (ymax - ymin) / (yspacing);

    int totalNodes = xnumNodes * ynumNodes;
    int totalElem = xnumElem * ynumElem;
    std::cout << "Total number of Nodes: " << totalNodes << "\n";
    std::cout << "Total number of Elements: " << totalElem << "\n";

    std::ofstream  fmesh, fcons;
    fmesh.open("mesh.dat", std::ios::app);
    fcons.open("mesh_constraints.dat", std::ios::app);

    double x, y;

    fmesh << "! elementShape quadrilateral" << "\n";
    fmesh << "! elemNumPoints 4" << "\n";
    fmesh << totalNodes << "\t" << totalElem << "\n";

    for (int j = 0; j < ynumNodes; j++) {
        for (int i = 0; i < xnumNodes; i++) {
            x = xmin + xspacing * i;
            y = ymin + yspacing * j;
            fmesh << x << "\t" << y << "\t" << 0. << "\n";
        }
    }

    for (int n = 0; n < ynumElem; n++) {
            std::vector<int> grid1, grid2;
            for (int p = 0; p < xnumNodes; p++) {
                int G1 = n*xnumNodes + p;
                int G2 = (n+1)*xnumNodes + p;
                grid1.push_back(G1);
                grid2.push_back(G2);
            } 
        for (int m = 0; m < xnumElem; m++) {
            int node1 = grid1[m];
            int node2 = grid1[m+1];
            int node3 = grid2[m+1];
            int node4 = grid2[m];
           fmesh << node1 << "\t" << node2 << "\t" << node3 << "\t" << node4 << "\n";
        }
    }
//    std::vector<int> FgridX, LgridX, FgridY, LgridY;
//    int gx = xnumNodes * (ynumNodes-1);
    int totalCons = 2*xnumNodes + 2*ynumNodes;

    fcons << totalCons << "\t" << 0 << "\n";
    for (unsigned a = 0; a < xnumNodes; a++) {
//        LgridX.push_back(gx+a);
 
        fcons << a << "\t" << 1 << "\n";
//        if (a = 0)
//            fcons << a << "\t" << 0 << "\n";
//        if (a = xnumNodes-1)
//            fcons << a << "\t" << 0 << "\n";
    }
    for (unsigned b = 1; b < ynumNodes; b++) {
        int gy1 = b*xnumNodes;
        int gy2 = gy1 + (xnumNodes-1);
        fcons << gy1 << "\t" << 0 << "\n";
        fcons << gy2 << "\t" << 0 << "\n";
    }

   for (unsigned a = 0; a < xnumNodes; a++) {
        unsigned xTop = totalNodes - 1 -a;
        fcons << xTop << "\t" << 1 << "\n";
   }

    fmesh.close();
    fcons.close();
    return 0;

}
