//
// Created by lgpda on 2/28/2025.
// The code contained herein is used to initialize, update, and interact with the vessel point sources of the simulation
//

#include "Vessel.h"
#include "Environment.h"
#include <cmath>


// INITIALIZE VESSEL
Vessel::Vessel(std::array<double, 2> loc, int idx, bool extraB): mt((std::random_device())()) {
    x=loc;
    id = idx;
    radius = 5;    //um
    sproutSpeed = 50;    //um/hr
    sproutBias = 0.65;
    canExtravasate = extraB;
    if (extraB) {vesselState = -2;}
    else {vesselState = -3;}
    canSprout = extraB;
    target = {0,0};
    CCNdistance = 100;     //um

}

// blank initializer
Vessel::Vessel(): mt((std::random_device())()) {
    x={0,0};
    radius = 10;    //um
    sproutSpeed = 50;    //um
    sproutBias = 0.65;
    canExtravasate = false;
    canSprout = false;
    target = {0,0};
    CCNdistance = 100;

}

void Vessel::neighboringCells(std::array<double, 2> otherX, int otherState){
    /*
     * code to determine the protumor neighbors (at this moment, only considering cancer cells)
     * as well as the closest cancer cell to vessel -- this determines direction of sprouting
     */
    double dis = calcDistance(otherX);
    if(dis <= 100){     // um
        protumorneighbors.push_back(otherX);
        if( otherState == 3 and dis < CCNdistance) {
            target = otherX;
            CCNdistance = dis;
        }
    }
}

void Vessel::neighboringVessels(std::array<double, 2> otherX) {
    /*
     *  code to determine vessel neighbors, this is used to determine if vessel can sprout
     */
    double dis = calcDistance(otherX);
    if(dis <= 100 and dis != 0.0){     // um
        vesselneighbors.push_back(otherX);
    }
}

void Vessel::setExtravasationBool(size_t step_count) {
    // lower threshold of cancer neighbors at which vessels allow extravasation
    // upper threshold of cancer neighbors at which vessel collapses
    int neighborThreshold = 90;
    if ((size(protumorneighbors) > 0)) {
        canExtravasate = true;
        vesselState = -2;
    }
    if ((size(protumorneighbors) > neighborThreshold) ) {
        canExtravasate = false;
        vesselState = -3;
    }
}

double Vessel::setSproutProb() {
    /*
     * function to set sprout probability based on the vessels distance from nearest cancer
     * exponential decay of sprout probability as distance increases
     */
    double const sproutGradient = 0.02;   //measure of how sharp the gradient is: 0.007 ~ half of max chance at 100 um
    double const beta = 0.0;       // distance from closest cancer cell in um where vessels are most likely to divide
    double const maxSproutProb = 0.1;     // max probability of dividing 1/hr
    double sproutProb;

    if (CCNdistance == 100) {
        sproutProb = 0.0;
    }
    else {
        sproutProb = exp(-abs(sproutGradient*CCNdistance-(sproutGradient*beta)))*maxSproutProb;
    }
    return sproutProb;
}

void Vessel::setCanSprout() {
    /*
     * determines boolean of whether vessel can sprout based on concentration of vessel neighbors
     */
    double maxVessels = 3.14*100*100*200/1000000;

    // std::cout << "max Vessels: " << maxVessels << std::endl;
    // where 80 or 161 is max MVD in vessel/mm2 and 100 is radius of neighborhood area
    if (size(vesselneighbors) > maxVessels or !canExtravasate) {
        canSprout = false;
    }
    else {canSprout = true;}
}

std::array<double, 3> Vessel::sprout() {
    /*
     * checks if vessel can sprout, sets sproutprob, tests sproutprob, and determines location of sprout if passed
     * vessels sprout in a random biased direction towards closest cancer neighbor
     * function returns location of new vessel {x,y,1}
     * OR {0,0,0} if no sprout occurred
     */

    std::array<double, 3> sprout {0,0,0};
    setCanSprout();
    if (!canSprout) {return sprout;}

    double sproutProb = setSproutProb();
    //double sproutProb = 1.1;
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    double test = dis(mt);
    if(test < sproutProb) {
        std::uniform_real_distribution<double> vect(-1.0, 1.0);
        std::array<double, 2> dx_random = {vect(mt), vect(mt)};
        //std::cout << "Random: (" << dx_random[0] << ", " << dx_random[1] << ")" << std::endl;
        dx_random = unitVector(dx_random);

        std::array<double, 2> target_direction = {target[0] - x[0],
                                      target[1] - x[1]};
        target_direction = unitVector(target_direction);
        std::array<double, 2> dx_movement = {0,0};

        for(int i=0; i<2; ++i){
            dx_movement[i] = sproutBias*target_direction[i] + (1- sproutBias)*dx_random[i];
        }
        dx_movement = unitVector(dx_movement);
        sprout[0] = x[0]+ sproutSpeed*dx_movement[0];
        sprout[1] = x[1]+ sproutSpeed*dx_movement[1];
        //std::cout << "Sprout: " << sprout[0] << ", " << sprout[1] << std::endl;
        sprout[2] = 1;
    }
    return sprout;
}




// OTHER FUNCTIONS
double Vessel::calcDistance(std::array<double, 2> otherX) {
    double d0 = (otherX[0] - x[0]);
    double d1 = (otherX[1] - x[1]);

    return sqrt(d0*d0 + d1*d1); // distance in um
}

std::array<double, 2> Vessel::unitVector(std::array<double, 2> v) {
    double norm = calcNorm(v);
    return {v[0]/norm, v[1]/norm};
}

double Vessel::calcNorm(std::array<double, 2> dx){
    return sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
}

