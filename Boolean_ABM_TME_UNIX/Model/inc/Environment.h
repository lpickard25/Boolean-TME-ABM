#ifndef IMMUNE_MODEL_ENVIRONMENT_H
#define IMMUNE_MODEL_ENVIRONMENT_H

#include <vector>
#include <algorithm>
#include <random>
#include "Cell.h"
#include <iostream>
#include <iomanip>      
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <filesystem> 
#include <array>
#include <omp.h>
#include "stb_image.h"
#include "stb_image_write.h"
#include <GL/glew.h>
#include <GLFW/glfw3.h>




class Environment{

public:
    Environment(std::string folder, std::string set, std::string tCellTrajectoryPath);
    //destructor needed
    bool visualize;
    void simulate(double tstep);
    void simulate(double tstep, std::vector<std::array<double, 3>> clusters);
    void initialize(std::string size, std::string xloc, std::string yloc);
private:
    void runCells(double tstep, size_t step_count);
    void neighborInfluenceInteractions(double tstep, size_t step_count);
    void internalCellFunctions(double tstep, size_t step_count);
    void recruitImmuneCells(double tstep, size_t step_count);
    std::array<double, 2> recruitmentLocation();
    void tumorSize(int i);
    void necrosis(double tstep);
    double calculateDiffusibles(std::array<double, 2> x);

    void save(double tstep, double tstamp);

    void saveTimeSeries();

    void loadParams();

    void initializeCells(std::array<double, 2> coords, double clusterRadius, int celltype);
    void calculateForces(double tstep);
    void updateTimeSeries();

    void printStep(double time);

    void createSet(std::vector<std::array<double, 3>> clusters, std::string folder, int Runs);
    
    double dt;

    // cell lists
    std::vector<Cell> cell_list;

    // time courses
    std::vector<int> cancerTS;
    std::vector<int> cd8TS;
    std::vector<int> cd8_suppTS;
    std::vector<int> cd4TS;
    std::vector<int> cd4_regTS;
    std::vector<int> m0TS;
    std::vector<int> m1TS;
    std::vector<int> m2TS;
    std::vector<int> radiusTS;

    /*
    t cell trajectory matrix: we can either represent as a vector of chars 
    where a char maps to a phenotypic state or an int where the int maps to 
    a phenotypic state
    */
    std::vector<std::string> tCellPhenotypeTrajectory_1;

    std::vector<std::vector<std::string>> tCellPhenotypeTrajectory; 

    

    // parameter lists
    std::vector<std::vector<double>> cellParams;
    std::vector<double> recParams;
    std::vector<double> envParams;

    std::string saveDir;
    int steps;
    std::vector<double> immuneCellRecRates;
    std::vector<int> immuneCellRecTypes;
    std::vector<double> immuneCells2rec;
    double recDist;
    double maxRecCytoConc;
    double tumorRadius;
    double necroticGrowth;
    double necroticRadius;
    double necroticForce;
    double necroticLimit;
    std::array<double, 2> tumorCenter;
    double recruitmentDelay;

    // environment params
    double simulationDuration;
    int day;

    std::mt19937 mt;

    // simulation params
    GLFWwindow* createWindow(int width, int height, const char* title, bool active);
    void drawModel(bool active, GLFWwindow* win, int scale);
    void updateWindow(GLFWwindow* win);
    std::vector<int> setColor(int cellType, int cellState);
    std::vector<int> setColor(int cellState);
    void drawCircle(float x, float y, float radius, std::vector<int> color, int scale);
    void drawCellOutlines(float x, float y, float radius,int scale);
    void saveToJPG(const char* filepath, GLFWwindow* win);
    std::array<std::vector<Cell>,8> createSubLists();

};

#endif //IMMUNE_MODEL_ENVIRONMENT_H
