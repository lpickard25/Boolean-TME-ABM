//
// Created by lgpda on 2/28/2025.
//


#ifndef VESSEL_H
#define VESSEL_H

#include <array>
#include <vector>
#include <random>

class Vessel {
    public:
    /*
     * Functions
     */
    explicit Vessel(std::array<double, 2> loc, int idx, bool extraP);

    Vessel();

    void neighboringCells(std::array<double, 2> loc, int otherState);
    void neighboringVessels(std::array<double, 2> otherX);

    void setExtravasationBool(size_t step_count);

    double setSproutProb();
    void setCanSprout();
    std::array<double, 3> sprout();

    double calcDistance(std::array<double, 2> otherX);
    std::array<double, 2> unitVector(std::array<double, 2> v);
    double calcNorm(std::array<double, 2> dx);

    /*
     * parameters
     */
    std::array<double, 2> x;

    double radius;
    std::vector<std::array <double, 2>> protumorneighbors;
    std::vector<std::array <double, 2>> vesselneighbors;
    std::array<double, 2> target;
    double CCNdistance;
    double sproutSpeed;
    double sproutBias;
    bool canExtravasate;
    int vesselState;
    int id;
    //double sproutProb;
    bool canSprout;


private:
    std::mt19937 mt;
};

#endif //VESSEL_H
