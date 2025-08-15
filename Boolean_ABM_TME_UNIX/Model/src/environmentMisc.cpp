#include "Environment.h"
#include "ModelUtil.h"

void Environment::initializeCells(std::array<double,2> coords, double clusterRadius, int celltype) {
    /*
     * places the initial tumor cluster(s)
     * cluster Radius is in # of cells
     * included celltype to allow initialization of other cell types -- used for testing/troubleshooting
     */


    if (celltype == 3) {
        // code needed for Tcell initialization
        int phenotypeIdx = getRandomNumber(tCellPhenotypeTrajectory.size());
        std::vector<std::string> trajec_phenotype = get2dvecrow(tCellPhenotypeTrajectory, phenotypeIdx);
        cell_list.push_back(Cell(coords, cell_list.size(), cellParams, 3, trajec_phenotype,0));
        for(int i=1; i<clusterRadius; ++i){
            double circumfrence = 2*i*cellParams[4][0]*3.1415;
            double nCells = circumfrence/cellParams[4][0];
            for(int j=0; j<nCells; ++j){
                double x = (i * cellParams[4][0] * cos(2 * 3.1415 * j / nCells))+coords[0];
                double y = (i * cellParams[4][0] * sin(2 * 3.1415 * j / nCells))+coords[1];
                cell_list.push_back(Cell({x,y}, cell_list.size(), cellParams, 3, trajec_phenotype,0));
            }
        }
    }
    else {
        cell_list.push_back(Cell(coords, cell_list.size(), cellParams, celltype, tCellPhenotypeTrajectory_1));

        double radiiCells = envParams[0];
        for(int i=1; i<clusterRadius; ++i){
            double circumfrence = 2*i*cellParams[4][0]*3.1415;
            double nCells = circumfrence/cellParams[4][0];
            for(int j=0; j<nCells; ++j){
                double x = (i * cellParams[4][0] * cos(2 * 3.1415 * j / nCells))+coords[0];
                double y = (i * cellParams[4][0] * sin(2 * 3.1415 * j / nCells))+coords[1];

                cell_list.push_back(Cell({x, y},cell_list.size(), cellParams, celltype, tCellPhenotypeTrajectory_1));
            }
        }
    }
}

void Environment::initializeVessels(int xDist, int yDist) {
    /*
     * function to randomly initialize vessel point source locations within a given area
     * initialized based on desired microvessel density (in vessel/mm^2)
     */
    const double vesselDensity = 14;  // in vessel/mm^2, 14 represents healthy breast tissue
    const int area = xDist*yDist*4; //in um squared
    double vesselNumber = area*vesselDensity/1000000;
    for(int i = 0; i < vesselNumber; ++i) {

        std::array<double, 2> location = makeVesselLoc(xDist, yDist);
        vessel_list.push_back(Vessel(location,vessel_list.size(),false));
    }
}

void Environment::checkAngiogenesis() {
    /*
     * function to determine if angiogenic switch has turned on
     * based on number of cancercells present in simulation
     */
    if (!angiogenesis) {
        int numCancer = countCancerCells();
        if (numCancer > 4000) {
            angiogenesis = true;
            std::cout << "Angiogenesis started on: " << steps/24 << std::endl;
        }
    }
}

std::array<double, 2> Environment::makeVesselLoc(int xDist, int yDist) {
    /*
     * function to create a random location within given boundaries
     */
    std::uniform_real_distribution<double> random(-1.0,1.0);

    double randomX = random(mt);
    double randomY = random(mt);
    double xloc = randomX*(xDist) + tumorCenter[0];
    double yloc = randomY*(yDist) + tumorCenter[1];

    return {xloc, yloc};
}



void Environment::recruitImmuneCells(double tstep,  size_t step_count, std::string mig_cells) {
    /*
     * function to initialize new immune cells based on the tumor burden and each immune cells
     * recruitment rate
     */
    // recruitment is scaled by number of cancer cells
    auto day = static_cast<double>(tstep*steps/24.0);
    if(day < recruitmentDelay){return;}

    int numC = cancerTS[steps - recruitmentDelay*24/tstep];

    // iterates through three recruitment rates
    for(int i=0; i<stod(mig_cells); ++i){
    //for(int i=0; i<1; ++i){
        // rate times the number of cancer cells
        double recRate = immuneCellRecRates[i]*static_cast<double>(numC);
        // number of cells to recruit is the rate times the tstep
        immuneCells2rec[i] += tstep*recRate;
        // iterates while there are still cells left to recruit
        while(immuneCells2rec[i] >= 1){
            // get recruitment location == the thing I need to change
            std::array<double, 2> recLoc = recruitmentLocation();
            if(std::isnan(recLoc[0]) || std::isnan(recLoc[1])) {
                std::cerr << " in recruitment" << recLoc[0] << recLoc[1] << std::endl;
                exit(-1);
            }
            
            //i==0 represents the idx associated with initializing a cd8 t cell
            if(i == 0){
                
                size_t phenotypeIdx = getRandomNumber(tCellPhenotypeTrajectory.size()); 

                std::vector<std::string> trajec_phenotype = get2dvecrow(tCellPhenotypeTrajectory, phenotypeIdx);
                if(trajec_phenotype.empty() || trajec_phenotype.size() == 0){
                    std::cerr << "WARNING recruitImmuneCells: t_cell_phenotype_Trajectory is empty!" << std::endl;
                    std::cerr << phenotypeIdx << std::endl; 
                    std::cerr << "ATTEMPTING FORCE RESOLVE... " << std::endl; 
                    phenotypeIdx = getRandomNumber(tCellPhenotypeTrajectory.size()); 
                    trajec_phenotype = get2dvecrow(tCellPhenotypeTrajectory, phenotypeIdx);
                    if(trajec_phenotype.empty() || trajec_phenotype.size() == 0){
                        std::cerr << "FORCE RESOLVE FAILED... " << std::endl; 
                    }
                    else{
                        std::cerr << "FORCE RESOLVE SUCCESSFUL... " << std::endl; 
                    }
                }
                

                cell_list.push_back(Cell(recLoc, 
                                    static_cast<int>(cell_list.size()), 
                                    cellParams, 
                                    immuneCellRecTypes[i],
                                    trajec_phenotype, 
                                    step_count)); 
            }
            else{
                cell_list.push_back(Cell(recLoc,
                                     static_cast<int>(cell_list.size()),
                                     cellParams,
                                     immuneCellRecTypes[i], tCellPhenotypeTrajectory_1));

            }
            
            immuneCells2rec[i] -= 1;
        }
    }
}

std::array<double, 2> Environment::recruitmentLocation() {
    /*
     * function to randomly select a viable point source location for an immune cell to be initialized at
     */
    // randomly selects a vessel from the list until it selects a vessel that allows extravasation
    std::uniform_int_distribution<int> random(0,vessel_list.size()-1);
    bool vesselFound = false;
    std::array<double, 2> location;
    while (!vesselFound) {

        int index = random(mt); //random number from 0-length of vessel list

        if (vessel_list[index].canExtravasate) {
            location = vessel_list[index].x;
            vesselFound = true;
        }
    }
    // adds slight random variability in the location to ensure cells do not have exact same locations
    // if cells have the exact same location can result in program erroring
    std::uniform_real_distribution<double> randomReal(-1.0,1.0);
    location[0] += randomReal(mt);
    location[1] += randomReal(mt);
    return location;


    /*
     * ORIGINAL FUNCTION
     * cells enter a random distance away from the tumor radius
     * cells enter at a random angle from the tumor center
     */
    //std::uniform_real_distribution<double> angle(-1.0, 1.0);
    // std::normal_distribution<double> angle(0.0,1.0);
    // std::array<double, 2> dx = {angle(mt),
    //                             angle(mt)};
    // double norm = sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
    //
    // std::uniform_real_distribution<double> loc(0.0, recDist);
    // double distance = loc(mt) + tumorRadius;
    //
    // return {distance*(dx[0]/norm), distance*(dx[1]/norm)};


    /*
     * FALL 2024 VERSION
     * Cells recruit randomly in a rectangle that is sized based on tumor and tumor center
     */

    // std::uniform_real_distribution<double> random(-1.0,1.0);
    // double randomX = random(mt);
    // double randomY = random(mt);
    // double xloc = randomX*(absMaxX+recDist) + tumorCenter[0];
    // double yloc = randomY*(absMaxY+recDist) + tumorCenter[1];
    // return {xloc, yloc};




    /*
     * recruit based on cytokine "concentrations". Recruit if a location is in a certain range
     */
    /*std::array<double, 2> x = {0.0, 0.0};
    double d = 0.0;
    std::normal_distribution<double> angle(0.0,1.0);
    std::uniform_real_distribution<double> distance_from_center(0.5*tumorRadius, 1.5*tumorRadius);
    std::uniform_real_distribution<double> recruitmentProbability(0.0, 1.0);
    std::array<double, 2> dx = {0.0, 0.0};
    double norm = 0;
    double distance = 0.0;
    bool recruited = false;
    while(!recruited){
        dx = {angle(mt),
              angle(mt)};
        norm = sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
        distance = distance_from_center(mt);
        x = {distance*dx[0]/norm, distance*dx[1]/norm};
        d = calculateDiffusibles(x);
        if(d <= maxRecCytoConc && recruitmentProbability(mt) < d){
            recruited = true;
        }
    }*/

    /*std::array<double, 2> x = {0.0, 0.0};
    double minD = 1e6;
    double maxD = 200.0;
    std::normal_distribution<double> angle(0.0,1.0);
    std::uniform_real_distribution<double> distance_from_center(0.5*tumorRadius, tumorRadius + 200);
    std::array<double, 2> dx = {0.0, 0.0};
    double norm = 0;
    double distance = 0.0;
    bool recruited = false;
    while(!recruited){
        dx = {angle(mt),
              angle(mt)};
        norm = sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
        distance = distance_from_center(mt);
        x = {distance*dx[0]/norm, distance*dx[1]/norm};
        for(auto & c : cell_list){
            if(c.state == 3){
                double dist = c.calcDistance(x);
                minD = std::min(dist, minD);
                if(minD < 50){
                    break;
                }

            }
        }
        if(minD >= 50.0 && minD < 200.0){
            recruited = true;
        } else{
            minD = 1e6;
        }
    }

    return x;*/
}

void Environment::tumorSize(int i) {
    tumorCenter = {0,0};
    // following values are needed for Fall 2024 version of immune cell recruitment
    double maxX = 0;
    double minX = 0;
    double maxY = 0;
    double minY = 0;
    double avgX = 0;
    double avgY = 0;
    // double numC = 0;
    for(auto &c : cell_list) {
        if(std::isnan(c.x[0]) || std::isnan(c.x[1])) {
            std::cerr << "error in tumorSize: " <<
                "type= " << c.type <<" state= " << c.state << c.x[0] << c.x[1] << std::endl;
            exit(-1);
        }
        if(c.type == 0) {
            maxX = std::max(maxX,c.x[0]);
            minX = std::min(minX,c.x[0]);
            maxY = std::max(maxY,c.x[1]);
            minY = std::min(minY,c.x[1]);
            // avgX += c.x[0];
            // avgY += c.x[1];
            // numC += 1;
        }
    }
    // avgX /= numC;
    // avgY /= numC;
    avgX = (maxX+minX)/2;
    avgY = (maxY+minY)/2;
    // absMaxX = std::max(abs(maxX),abs(minX));
    // absMaxY = std::max(abs(maxY),abs(minY));

    tumorCenter = {avgX, avgY};

    tumorRadius = 0;
    absMaxX = 0;
    absMaxY = 0;
    for(auto& c : cell_list){
        if(c.type == 0){
            tumorRadius = std::max(tumorRadius, c.calcDistance(tumorCenter));
            // following values are needed for fall 2024 version of immune cell recruitment
            absMaxX = std::max(absMaxX, std::fabs(c.x[0]-avgX));
            absMaxY = std::max(absMaxY, std::fabs(c.x[1]-avgY));
        }
    }

    //originally commented out

    // tumorCenter = {0.0, 0.0};
    // double avgX = 0;
    // double avgY = 0;
    // double numC = 0;
    // std::vector<int> mainTumorMass;
    // for(auto & c : cell_list){
    //     if(c.type == 0){
    //         int numNeighbors = 0;
    //         for(auto &c2 : cell_list){
    //             if(c2.type == 0 && c.calcDistance(c2.x) < 50.0){
    //                 ++numNeighbors;
    //             }
    //         }
    //         if(numNeighbors > 10){
    //             numC++;
    //             avgX += c.x[0];
    //             avgY += c.x[1];
    //             mainTumorMass.push_back(1);
    //         } else{
    //             mainTumorMass.push_back(0);
    //         }
    //     } else{
    //         mainTumorMass.push_back(0);
    //     }
    // }

    // avgX /= numC;
    // avgY /= numC;

    // tumorCenter = {avgX, avgY};

    //----END originally commented out

    // tumorRadius = 0.0;
    // for(int i=0; i<cell_list.size(); ++i){
    //     if(mainTumorMass[i] == 1){
    //         if(cell_list[i].type == 0){
    //             tumorRadius = std::max(tumorRadius, cell_list[i].calcDistance(tumorCenter));
    //         }
    //     }
    // }

}

void Environment::necrosis(double tstep) {
    int nCancer = 0;
    for(auto &c : cell_list){
        if(c.type == 0){
            ++nCancer;
        }
    }

    necroticRadius += tstep*nCancer*necroticGrowth;
    necroticRadius = std::max(std::min(necroticRadius, tumorRadius-necroticLimit), 0.0);
}

double Environment::calculateDiffusibles(std::array<double, 2> x) {
    double d = 0;
    for(auto & c: cell_list){
        if(c.state == 3){
            // is a cancer cell
            double alpha = -log2(c.probTh);
            double lambda = alpha*0.693/c.influenceRadius;
            double dist = c.calcDistance(x);
            d = 1 - (1 - d)*(1 - exp(-lambda*dist));
        }
    }

    return d;
}