#include <cstring>

#include "Environment.h"
#include "ModelUtil.h"
#include <GLFW/glfw3.h>



namespace fs = std::filesystem;
// initialize killCount so that it can be updated in other functions later
int Environment::killCount = 0;


void read_trajectory_csv(std::vector<std::vector<std::string>>& csv_data, std::string fpath){
    std::ifstream myfile(fpath);

    //try to read in the csv file
    if(!myfile.is_open()){
        std::cerr << "Error: Unable to open file at " << fpath << std::endl; 
        return; 
    }
    char ch;
    std::string cell;
    std::vector<std::string> current_row;
    bool inside_quotes = false;
    while (myfile.get(ch)) {
        if (ch == '"') {
            inside_quotes = !inside_quotes;
        } else if (ch == ',' && !inside_quotes) {
            current_row.push_back(cell);
            cell.clear();
        } else if (ch == '\n' && !inside_quotes) {
            current_row.push_back(cell);
            cell.clear();
            csv_data.push_back(current_row);
            current_row.clear();
        } else {
            cell += ch;
        }
    }
    // Close the file
    myfile.close(); 
}

std::vector<std::string> get_phenotype(std::vector<std::vector<std::string>>& trajectory_2d_vec){ 
    std::vector<std::string> res; 

    size_t max_col = trajectory_2d_vec[0].size(); 

    for(size_t i = 1; i < trajectory_2d_vec.size(); ++i){ 
        //start at i=1 because the first value is the title, change this if need be
        res.push_back(trajectory_2d_vec[i][max_col - 1]); 
    }

    return res; 

}


Environment::Environment(std::string folder, std::string set, std::string tCellTrajectoryPath): mt((std::random_device())()) {
    /*
     * initialize a simulation environment
     * -----------------------------------
     * loads the three parameter files
     *  cell parameters
     *  environment parameters
     *  recruitment parameters
     *
     * set environment variables to their respective values
     */

    saveDir = folder+"\\set_"+set;
    loadParams();

    for(int i=0; i<3; ++i){
        immuneCellRecRates.push_back(recParams[i]);
        immuneCells2rec.push_back(0.0);
    }
    std::cout << "parameters loaded" << std::endl;
    immuneCellRecTypes = {3, 1, 2}; // {CD8, macrophage, CD4} -> same order as RecRates

    recDist = recParams[3];
    //maxRecCytoConc = recParams[3];
    recruitmentDelay = recParams[4];

    simulationDuration = envParams[1];
    necroticGrowth = envParams[2];
    necroticForce = envParams[3];
    necroticLimit = envParams[4];

    tumorCenter = {0,0};
    tumorRadius = 0;
    necroticRadius = 0;

    steps = 0;
    day = 0;


    dt = 0.005;

    angiogenesis = false;


    // attempt to load in t cell trajectory files by iterating over files in tCellTrajectoryPath
    printf("Attempting to read trajectories located in %s...\n", tCellTrajectoryPath.c_str());

    
    for (const auto& entry : std::filesystem::directory_iterator(tCellTrajectoryPath)){
        //make sure we are only reading csv files 
        if(entry.path().extension().string() == ".csv"){
            
            std::vector<std::vector<std::string>> trajec_csv;

            read_trajectory_csv(trajec_csv, entry.path().string()); 
            std::vector<std::string> ptype = get_phenotype(trajec_csv); 
            tCellPhenotypeTrajectory.push_back(ptype); 
        } 
        else{
            std::cout << entry.path().string() << " cannot be read as a trajectory, if this is unexpected please confirm the integrity of the .csv file \n"; 
        }
    }

    printf("Finished reading trajectories located in %s...\n", tCellTrajectoryPath.c_str());
    
    //attempt to load in t cell trajectory file 
    std::string trajecPath =  "t_cell_trajectory/1Tcell_Sim_ABM.csv"; 
    // std::string trajecPath2 = tCellTrajectoryPath + "/2Tcell_Sim_ABM.csv"; 
    // std::string trajecPath3 = tCellTrajectoryPath + "/3Tcell_Sim_ABM.csv"; 

    std::vector<std::vector<std::string>> trajec_csv_1;
    
    // read_trajectory_csv(trajec_csv_1, trajecPath1); 
    // // if we only care about the phenotype we can simply create a std::vector<std::string or char> of phenotypes 

    std::cout << "new code start" << std::endl;
    std::vector<std::string> phenotype_trajec_1;
    tCellPhenotypeTrajectory_1 = phenotype_trajec_1; 
    
    std::cout << "Constructor done" << std::endl;
}

void Environment::initialize(std::string numClusters, std::string sizeClusters, std::string distance) {
    /*
     * function to intialize cancer clusters based on inputs from program arguments
     * can either initialize 1 cluster of a given size at (0,0)
     * or 2 clusters of given sizes at a given distance apart
     * or 3 clusters of given sizes at a given distance apart in an equilateral triangle
     * room here to make code capable of more complicated initializations, based on user input
     *
     * Eg: numClusters = 2 clusterSize = 2.7 distance = 1000
     * - two clusters, one with a radius of 2 cells, the other with a radius of 7 cells, 1000 um apart
     */

    double radius1;
    double radius2;
    double radius3;
    if (numClusters == "1") {
        radius1 = stod(sizeClusters);
        //std::cout << "radius: " << radius << std::endl;
        initializeCells({0.0,0.0},radius1,0);
        // special case for a radius of 1 -- initializes 3 cancer cells instead of just one
        // to help ensure cluster survives
        if (radius1 == 1.0) {
            std::array<double,2> center3 = {0.0, 0.0};
            center3[0] -= 20;
            initializeCells(center3,radius1,0);
            center3[1] += 20;
            initializeCells(center3,radius1,0);
        }
    }
    else if (numClusters == "2") {
        char *clusterSize = strdup(sizeClusters.c_str());
        char *r1 = strtok(clusterSize, ".");
        char *r2 = strtok(NULL, ".");
        radius1 = atof(r1);
        radius2 = atof(r2);
        //std::cout << sizeClusters  << std::endl;
        std::array<double,2> center1 = {stod(distance)/-2,0};
        std::array<double,2> center2 = {stod(distance)/2,0};
        //std::cout << "radii: " << radius1 << ", " << radius2 << std::endl;
        //std::cout << "xlocs: " << center1[0] << ", " << center2[0] << std::endl;
        initializeCells(center1, radius1, 0);

        initializeCells(center2, radius2, 0);
        if (radius2 == 1.0){
            std::array<double,2> center3 = center2;
            center3[0] += 20;
            initializeCells(center3,radius2,0);
            center3[1] += 20;
            initializeCells(center3,radius2,0);
        }
        if (radius1 == 1.0) {
            std::array<double,2> center3 = center1;
            center3[0] -= 20;
            initializeCells(center3,radius1,0);
            center3[1] += 20;
            initializeCells(center3,radius1,0);
        }
    }
    else if (numClusters == "3") {
        char *clusterSize = strdup(sizeClusters.c_str());
        char *r1 = strtok(clusterSize, ".");
        char *r2 = strtok(NULL, ".");
        char *r3 = strtok(NULL, ".");
        radius1 = atof(r1);
        radius2 = atof(r2);
        radius3 = atof(r3);
        std::array<double,2> center1 = {stod(distance)/-2,stod(distance)*sqrt(3)/4};
        std::array<double,2> center2 = {stod(distance)/2,stod(distance)*sqrt(3)/4};
        std::array<double,2> center3 = {0,stod(distance)*sqrt(3)/-4};
        initializeCells(center1, radius1, 0);
        initializeCells(center2, radius2, 0);
        initializeCells(center3, radius3, 0);
    }
    else {
        initializeCells({0.0,0.0},5 , 0);
    }
    //double cells = radius1*radius1*2.8+radius2*radius2*2.8;
    std::cout << "initialize cells finished" << std::endl;
    // follow values determine area in which vessels are initialized
    // I found that 4000 by 4000 contains a single centered cluster well, after 24 simulation days
    xdist = 4000+ stod(distance)/2;
    ydist = 4000;
    initializeVessels(xdist,ydist);
    std::cout << "initialize vessels finished" << std::endl;
}

void Environment::simulate(double tstep, std::string mig_cells) {
    /*
     * initializes and runs a simulation
     * ---------------------------------
     * place initial tumor and vessels
     * run simulation loop
     *  recruit immune cells
     *  run cell functions
     *  run vessel functions
     * ends once time limit is reached or there are no more cancer cells
     */

    // included for visualization
    int model_time = 0;
    int saveInterval = 24;  //24 to save every day, 1 to save every hour
    int scale = 2;
    int winWidth = 2000/scale;
    int winHeight = 2000/scale;
    GLFWwindow* win = createWindow(winWidth, winHeight, "ABM", true);


    //more code for visualization
    int dayNum = steps * tstep / 48;
    std::string filepath = saveDir +"\\images\\tumor_day_0.jpg";
    drawModel(visualize, win, scale, killCount);
    const char* cstr = filepath.c_str();
    saveToJPG(cstr, win);

    tumorSize(steps);
    //std::cout << "Initial Tumor Center: ";
    // for (int i = 0; i <2; i++) {
    //     std::cout << tumorCenter[i]<<",";
    // }
    //std::cout << "; Tumor Radius: "<<tumorRadius << " necrotic radius: " << necroticRadius<<std::endl;
    save(tstep, steps*tstep);
    updateTimeSeries();

    std::cout << "starting simulations...\n";

    while(tstep*steps/24 < simulationDuration) {
        recruitImmuneCells(tstep, tstep*steps, mig_cells);
        checkAngiogenesis();
        runCells(tstep, tstep*steps);
        tumorSize(steps);
        necrosis(tstep);

        steps += 1;
        //printStep(steps * tstep);
        updateTimeSeries();
        model_time = steps;
        if (fmod(steps * tstep, saveInterval) == 0) {

            // creat openGl window
            drawModel(visualize, win, scale, killCount);
            int dayNum = steps * tstep / saveInterval;
            std::string filepath = saveDir + "\\images\\tumor_day_"+std::to_string(dayNum)+".jpg";
            // convert std::string to const char*
            const char* cstr = filepath.c_str();
            saveToJPG(cstr, win);

        }
        if (fmod(steps * tstep, saveInterval) == 0) {
            // save every saveInterval
            save(tstep, steps*tstep);
            std::cout << "saved step = " << steps*tstep << std::endl;
            std::cout << "Vessel count = " << vessel_list.size() << std::endl;
            std::cout << "Cancer cells: " << countCancerCells() << std::endl;
        }

        // commented out following code so that simulation runs even if cancer is not present
        // int numC = 0;
        // for (auto &c: cell_list) {
        //     if (c.type == 0) {
        //         numC++;
        //     }
        // }
        // if (numC == 0) {
        //     save(tstep, steps*tstep);
        //     break;
        //}
    }
    // save csv after simulation is finished
    saveTimeSeries();
}

