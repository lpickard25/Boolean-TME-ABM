

#include <iostream>
#include "Environment.h"
#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

/*
 * CHANGE summing influences to multiplying (1 - inf)
 */

/*
 * MAIN THINGS TO LOOK INTO RIGHT NOW
 * ----------------------------------
 * - impact of macrophages on CD8 infiltration: https://www.pnas.org/doi/10.1073/pnas.1720948115
 * - how M1 and M2 should impact CD8: M2 is easy to find (Petty and Yang, Tumor-associated macrophages: implications in cancer immunotherapy, 2017)
 *      - M2
 *          - reduced CD8 proliferation via cytokine secretion - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6715643/ (media from m2 macrophages)
 *
 *      - M1
 *          - increased T cell infiltration and CD8 activation - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6715643/ (but is this from M1, or simply the removal of M2?)
 *          - conversion to M1 "regulated T-cell response by relieving immunosuppression" - https://pubs.acs.org/doi/full/10.1021/acsami.1c07626
 *          - direct tumor killing - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7751482/
 *   but M1 is a bit more difficult
 *   - changes in CD8 killing effect
 *   - changes in CD8 cytokine secretion
 *   - changes in CD8 proliferation
 *   - changes in CD8 recruitment
 *   - changes in tumor proliferation
 * - M1 cytotoxic effect
 * - CD8 proliferation
 * - cell recruitment
 * - inclusion of Th2
 * - functions of CD4
 *      - Th1
 *          - IFN-y secretion (https://www.nature.com/articles/s41417-020-0183-x)
 *          - IL-2 secretion: drive CD8 effector function, differentiation, and proliferation (https://www.nature.com/articles/s41417-020-0183-x)
 *          - increase number of CD8 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1403291/)
 *          - can replace exogenous IL-2, due to IL-2 production (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1403291/)
 *          - removing Treg is not enough to remove tumor, T help must be provided (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1403291/)
 *          - specific impacts on CD8 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2970736/)
 *              - increased CD8 recruitment to tumor site (IFN-y)
 *              - increased CD8 survival (IL-2)
 *              - increased CD8 proliferation at tumor site (IL-2)
 *              - increased CD8 granzyme B (IL-2)
 *      - Treg
 *          - suppress effects of exogenous IL-2 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1403291/)
 *
 * - suppressed CD8 IFN-y
 *
 * MODEL FRAMEWORK THINGS
 * ----------------------
 * - what is going on with migration bias?
 */

/*
 * MODEL OF TUMOR-IMMUNE INTERACTIONS
 * ----------------------------------
 * - PUT MODEL DESCRIPTION HERE
 */

/*
 * THERAPY (future)
 * ----------------
 * - in whole tumor mode, assume sufficient vascularization outside of tumor, then decrease effect going into the tumor using exponential decay
 *
 * POTENTIAL ANALYSES
 * ------------------
 * - don't let cells die of age. at end of simulation, output distributions of cancer cell phenotype/distance to nearest immune cell that induces that phenotype
 */

int main(int argc, char **argv) {
    /*
     * takes all the program arguments, runs the genParams.py, creates the environment object
     * initializes, begins and times the simulation
     */
    //std::cout << "argument 0: " << argv[0] << std::endl;
    std::string folder = argv[1];  //string: SimulationSet
    std::string set = argv[2];     //int: 1
    std::string pST = argv[3];     //int: 1
    std::string dp_fac = argv[4];  //float: 1
    std::string kp_fac = argv[5];  //float: 1
    std::string numClusters = argv[6];  //int:2 (stod() is used on it)
    std::string sizeClusters = argv[7]; //char: 3.7 (means one cluster is 3 and the other is 7 cell in radius
    std::string distance = argv[8];    //int: 1000 (stod() is used on it)
    std::string rec_rate_fac = argv[9]; //float: 1
    std::string mig_fac = argv[10];  //float: 5
    std::string mig_cells = argv[11];  //int: either 3 for all immune cells or 1 for cd8 t cells only (stod())

    std::string str = "python genParams.py "+folder+" "+set + " "+ pST + " " + dp_fac + " " + kp_fac
      + " " + rec_rate_fac + " " + mig_fac;
    const char* command = str.c_str();
    std::system(command);

    std::cout << "params done" << std::endl;


    Environment model(folder, set, "Model\\phenotype_out\\"); //can replace with a directory representing any other phenotype state
    model.visualize = true;
    std::cout << "initialize begun " << std::endl;
    model.initialize(numClusters, sizeClusters, distance);
    std::cout << "initialize finish " << std::endl;
    double start = omp_get_wtime();
    model.simulate(1, mig_cells);
    double stop = omp_get_wtime();
    std::cout << "Duration: " << (stop-start)/(60*60) << std::endl;

    return 0;
}
