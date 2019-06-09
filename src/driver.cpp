
//======================================
// 1D Eiuler approximate Riemann sovler
#include "../include/EulerShockTube1D.h"

//======================================
// structured grids 2D class 
#include "../include/gridGen2D.h"


int main(){

    if (false){
        // Solver solver;
        // solver.Euler1D();
        // solver.output();
        driverEuler1D();
    }else{
        driverGrid2D();
    }
    return 1;
}