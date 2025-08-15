#include "Environment.h"

void Environment::printStep(double time) {
    int numM = 0;
    int numT8 = 0;
    int numT4 = 0;
    int numC = 0;

    for(auto &cell : cell_list){
        if(cell.type == 1){
            numM++;
        } else if(cell.type == 3){
            numT8++;
        } else if(cell.type == 2){
            numT4++;
        } else if(cell.type == 0){
            numC++;
        }
    }
    std::cout << std::fixed << std::setprecision(5);
    std::cout << "Time: " << std::setw(10) << (time / 24) << " | cancer: " << std::setw(10) << numC 
          << " | cd8: " << std::setw(10) << numT8 << " | cd4: " << std::setw(10) << numT4 
          << " | macrophage: " << std::setw(10) << numM << std::endl;
}

void Environment::updateTimeSeries() {
    /*
     * function to save time series of number of each state of cell present at every times step
     * also tracks total number of kills
     * these lists get saved as a csv after the end of the simulation
     */
    int numT8_activated = 0;
    int numT8_suppressed = 0;
    int numT4_reg = 0;
    int numT4_help = 0;
    int numC = 0;

    for(auto &cell : cell_list){
        if(cell.type == 3 && cell.state == 6){
            numT8_activated++;
        }else if(cell.type == 3 && cell.state == 7) {
            numT8_suppressed++;
        }
        else if(cell.type == 2 && cell.state == 4) {
            numT4_help++;
        }
        else if(cell.type == 2 && cell.state == 5) {
            numT4_reg++;
        }
        else if(cell.type == 0){
            numC++;
        }
    }

    cancerTS.push_back(numC);
    cd8TS.push_back(numT8_activated);
    cd8_suppTS.push_back(numT8_suppressed);
    cd4TS.push_back(numT4_help);
    cd4_regTS.push_back(numT4_reg);
    killCountTS.push_back(killCount);

    int m0 = 0;
    int m1 = 0;
    int m2 = 0;
    for(auto &c : cell_list){
        if(c.type == 1) {
            if (c.state == 0) { m0++; }
            if (c.state == 1) { m1++; }
            if (c.state == 2) { m2++; }
        }
    }

    m0TS.push_back(m0);
    m1TS.push_back(m1);
    m2TS.push_back(m2);

    radiusTS.push_back(tumorRadius);
}

int Environment::countCancerCells() {
    int numC = 0;
    for(auto &cell : cell_list){
        if(cell.type == 0){
            numC++;
        }
    }
    return numC;
}
