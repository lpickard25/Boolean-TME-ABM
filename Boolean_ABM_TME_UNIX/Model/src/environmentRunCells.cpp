#include "Environment.h"
#include "ModelUtil.h"

void Environment::neighborInfluenceInteractions(double tstep, size_t step_count) {

    /*
     * FIRST LOOP
     * - determine neighbors
     * - determine influences on a cell
     * - perform indirect interactions
     *
     * SECOND LOOP
     * - direct interactions
     *
     * THIRD LOOP
     * - differentiate
     */

#pragma omp parallel for
    for(int i=0; i<cell_list.size(); ++i){
        // reset neighborhood and influence
        if(std::isnan(cell_list[i].x[0]) || std::isnan(cell_list[i].x[1])) {
            std::cerr << "error in neighborhoodInfInt: " <<cell_list[i].type << cell_list[i].x[0] << cell_list[i].x[1] << std::endl;
            exit(-1);
        }
        cell_list[i].neighbors.clear();
        //reset cancer distance and presence of neighbors
        cell_list[i].cancerDistance = 100;
        cell_list[i].cancerNeighbor = false;
        cell_list[i].clearInfluence();

        for(auto &c : cell_list){
            // assume that a cell cannot influence itself
            if(cell_list[i].id != c.id){
                cell_list[i].neighboringCells(c.x, c.id, c.state);
                cell_list[i].addInfluence(c.x, c.influenceRadius, c.state);
                //cell_list[i].addChemotaxis(c.x, c.influenceRadius, c.type);
            }

        }

        cell_list[i].indirectInteractions(tstep, step_count);

    }

#pragma omp parallel for
    for(int i=0; i<cell_list.size(); ++i){
        for(auto &c : cell_list[i].neighbors){
            cell_list[i].directInteractions(cell_list[c].state,
                                            cell_list[c].x,
                                            cell_list[c].directInteractionProperties(cell_list[i].state, step_count),
                                            tstep);
        }
    }
#pragma omp parallel for
    for(int i=0; i<cell_list.size(); ++i){
        cell_list[i].differentiate(tstep);
    }


}


void Environment::calculateForces(double tstep) {
    /*
     * 1. Calculate total force vector for each cell
     * 2. Resolve forces on each cell
     * 3. Determine current overlap for each cell
     * 4. Determine if each cell is compressed
     */

    // divide tstep into smaller steps for solving
    // only solve forces between neighboring cells to improve computation time
    int Nsteps = static_cast<int>(tstep/dt);
    for(int q=0; q<Nsteps; ++q) {
        // migrate first
#pragma omp parallel for
        for(int i=0; i<cell_list.size(); ++i){
            cell_list[i].migrate(dt, tumorCenter);
        }

        // calc forces
#pragma omp parallel for
        for(int i=0; i<cell_list.size(); ++i){
            //std::cout << "cell: " << i << std::endl;

            for(auto &c : cell_list[i].neighbors) {
                //std::cout << "neighbor: " << c << std::endl;
                //if (cell_list[i].type != 0 and cell_list[c].type != 0){
                    cell_list[i].calculateForces(cell_list[c].x, cell_list[c].radius, cell_list[c].type, cell_list[c].id);
                //}
            }
        }

#pragma omp parallel for
        for(int i=0; i<cell_list.size(); ++i){
            cell_list[i].resolveForces(dt, tumorCenter, necroticRadius, necroticForce);
        }
    }

        // calculate overlaps and proliferation states
#pragma omp parallel for
        for(int i=0; i<cell_list.size(); ++i){
            for(auto &c : cell_list[i].neighbors){
                cell_list[i].calculateOverlap(cell_list[c].x, cell_list[c].radius);
            }
            cell_list[i].isCompressed();

        }
    }

void Environment::internalCellFunctions(double tstep, size_t step_count) {
    /*
     * cell death via aging
     * cell proliferation
     * remove cell if out of bounds
     */
    const int numCells = cell_list.size();
    for(int i=0; i<numCells; ++i){
        //if (cell_list[i].type != 0) {
            cell_list[i].age(tstep, step_count);
        //}
        // if in necrotic core, die
        if(cell_list[i].calcDistance(tumorCenter) < necroticRadius){
            cell_list[i].state = -1;
            std::cout << "cell died" << std::endl;
        }

        cell_list[i].prolifState();
        std::array<double, 3> newLoc = cell_list[i].proliferate(tstep);

        if(newLoc[2] == 1){
            if(std::isnan(newLoc[0]) || std::isnan(newLoc[1])) {
                std::cerr << cell_list[i].type << " in internal cell function" << newLoc[0] << newLoc[1] << std::endl;
                exit(-1);
            }
            if(cell_list[i].type == 3){
                int phenotypeIdx = getRandomNumber(tCellPhenotypeTrajectory.size());
                std::vector<std::string> trajec_phenotype = get2dvecrow(tCellPhenotypeTrajectory, phenotypeIdx);
                if(trajec_phenotype.empty() || trajec_phenotype.size() == 0){
                    std::cerr << "WARNING INTERNAL CELL FUNCTIONS: t_cell_phenotype_Trajectory is empty!" << std::endl;
                }
                cell_list.push_back(Cell({newLoc[0], newLoc[1]},
                                     cell_list.size(),
                                     cellParams,
                                     cell_list[i].type, trajec_phenotype, step_count));
            }
            else{
                cell_list.push_back(Cell({newLoc[0], newLoc[1]},
                                     cell_list.size(),
                                     cellParams,
                                     cell_list[i].type, tCellPhenotypeTrajectory_1));

            }
            cell_list[cell_list.size() - 1].inherit(cell_list[i].inheritanceProperties());
        }
    }

    // remove dead cells
    std::vector<int> dead;
    for(int i=0; i<cell_list.size(); ++i){
        if(cell_list[i].state == -1){
            dead.push_back(i);
        }
    }
    std::reverse(dead.begin(), dead.end());
    for(auto &i : dead){
        cell_list.erase(cell_list.begin()+i);
    }

    for(int i=0; i<cell_list.size(); ++i){
        cell_list[i].updateID(i);
        if(cell_list[i].state == -1){
            throw std::runtime_error("Environment::internalCellFunctions -> dead cell not removed");
        }
    }
}

void Environment::runVessels(size_t step_count) {
    /*
     * iterates over each vessel and determines its protumor neighbors, vessel neighbors,
     * and closest cancer neighbor (CCN)
     * determined whether vessel allows extravasation based on these
     */

#pragma omp parallel for
    for(int i=0; i<vessel_list.size(); ++i) {
        vessel_list[i].protumorneighbors.clear();
        vessel_list[i].vesselneighbors.clear();
        vessel_list[i].CCNdistance = 100;
        for(auto &c : cell_list) {
            if (c.state == 3) {
                vessel_list[i].neighboringCells(c.x, c.state);
            }
        }
        vessel_list[i].setExtravasationBool(step_count);
        for(auto &v : vessel_list) {
            if (vessel_list[i].id != v.id) {
                vessel_list[i].neighboringVessels(v.x);
            }
        }
        //std::cout << "Vessel neighbors: " << size(vessel_list[i].vesselneighbors) << std::endl;
    }
}

void Environment::runSprout(size_t step_count) {
    /*
     * function only runs if angiogensis is true
     * iterates over vessels and determines if/where it sprouts
     * potential for parallelization here
     */
    std::vector<std::array<double, 3>> newVessels(size(vessel_list));
    for(int i=0; i<vessel_list.size(); ++i) {
        std::array<double, 3> newVessel = vessel_list[i].sprout();
        newVessels[i] = newVessel;
    }

    for(auto &v :newVessels) {
        if (v[2] == 1) {
            vessel_list.push_back(Vessel({v[0],v[1]}, vessel_list.size(), false));
        }
    }
}



void Environment::runCells(double tstep, size_t step_count) {
    neighborInfluenceInteractions(tstep, step_count);

    calculateForces(tstep);

    internalCellFunctions(tstep, step_count);

    runVessels(step_count);

    if (angiogenesis) {
        runSprout(step_count);
    }

}