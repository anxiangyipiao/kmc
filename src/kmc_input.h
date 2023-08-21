#ifndef _KMC_INPUT_H_
#define _KMC_INPUT_H_
#include <random>
#include <string>
#include "yaml-cpp/yaml.h"

extern std::mt19937 e;
extern std::uniform_real_distribution<float> u;


extern int   v1nbr_bcc[2][3][8];
extern int   v2nbr_bcc[3][6];
extern int   v1nbr_fcc[3][4][4][4];

/*<summary>
par_ltc: a string that specifies the lattice structure of the system
nx, ny, nz: integers that specify the number of unit cells in the x, y, and z directions, respectively
parallelx, parallely, parallelz : integers that specify the number of processors to use in each direction for parallel computing
par_compB, par_compV, par_compC, par_compD : doubles that specify the compression factors for different types of interactions in the system
Mutiple, MutipleSize, Paralle : boolean variables that specify whether to use multiple replicas of the system, whether to change the size of the replicas, and whether to use parallel computing, respectively
par_time : a double that specifies the total simulation time in units of picoseconds
time_conf : a double that specifies the time interval for recording the system configuration
par_step : an unsigned long long integer that specifies the number of time steps in the simulation
step_log : an unsigned long long integer that specifies the time interval for recording log data
outposcar : an unsigned long long integer that specifies the time interval for outputting the system configuration in VASP POSCAR format
read_file : a boolean variable that specifies whether to read in a configuration file
filepath : a string that specifies the path of the configuration file to read in
filepackage : a string that specifies the package to use for reading in the configuration file
par_temp : a double that specifies the temperature of the system in units of Kelvin
par_beta : a double that specifies the inverse temperature of the system
par_dis_rec : a double that specifies the distance between atoms that should trigger the energy cutoff
par_muvA : a double that specifies the viscosity of the system
par_radius_start : an integer that specifies the starting number of atoms
par_eSPA... par_eSPV2V : doubles that specify the energy parameters for different types of interactions in the system.
*/
class SimulationParameters {

public:
    std::string par_ltc;
    int nx, ny, nz, parallelx, parallely, parallelz;
    double par_compB, par_compV, par_compC, par_compD;
    double par_time, time_conf;
    unsigned long long par_step, step_log, outposcar;
    bool read_file;
    bool Mutiple, MutipleSize,Paralle;
    std::string filepath, filepackage;
    double par_temp, par_beta, par_dis_rec, par_muvA;
    int par_radius_start;
    double par_eSPA, par_eSPB, par_eSPC, par_eSPD;
    double par_eSPA1A, par_eSPA2A, par_eSPA1B, par_eSPA2B, par_eSPA1V, par_eSPA2V, par_eSPA1C, par_eSPA2C, par_eSPA1D, par_eSPA2D;
    double par_eSPB1B, par_eSPB2B, par_eSPB1V, par_eSPB2V, par_eSPB1C, par_eSPB2C, par_eSPB1D, par_eSPB2D;
    double par_eSPC1C, par_eSPC2C, par_eSPC1D, par_eSPC2D, par_eSPC1V, par_eSPC2V;
    double par_eSPD1D, par_eSPD2D, par_eSPD1V, par_eSPD2V;
    double par_eSPV1V, par_eSPV2V;
    SimulationParameters(std::string filename) {
        YAML::Node config = YAML::LoadFile(filename);

        // System parameters
        par_ltc = config["par_ltc"].as<std::string>();
        nx = config["nx"].as<int>();
        ny = config["ny"].as<int>();
        nz = config["nz"].as<int>();
        parallelx = config["parallelx"].as<int>();
        parallely = config["parallely"].as<int>();
        parallelz = config["parallelz"].as<int>();

        par_compB = config["par_compB"].as<double>();
        par_compV = config["par_compV"].as<double>();
        par_compC = config["par_compC"].as<double>();
        par_compD = config["par_compD"].as<double>();
        Mutiple = config["Mutiple"].as<bool>();
        MutipleSize = config["MutipleSize"].as<bool>();
        Paralle = config["paralle"].as<bool>();
        // Simulation time parameters
        par_time = config["par_time"].as<double>();
        time_conf = config["time_conf"].as<double>();
        par_step = config["par_step"].as<unsigned long long>();
        step_log = config["step_log"].as<unsigned long long>();
        outposcar = config["Outposcar"].as<unsigned long long>();
        // Read File
        read_file = config["read_file"].as<bool>();
        filepath = config["filepath"].as<std::string>();
        filepackage = config["filepackage"].as<std::string>();
        // Kinetic parameters
        par_temp = config["par_temp"].as<double>();
        par_beta = config["par_beta"].as<double>();
        par_dis_rec = config["par_dis_rec"].as<double>();
        par_muvA = config["par_muvA"].as<double>();

        // Sets the starting number of atoms
        par_radius_start = config["par_radius_start"].as<int>();

        // Energy parameters
        par_eSPA = config["par_eSPA"].as<double>();
        par_eSPB = config["par_eSPB"].as<double>();
        par_eSPC = config["par_eSPC"].as<double>();
        par_eSPD = config["par_eSPD"].as<double>();

        par_eSPA1A = config["par_eSPA1A"].as<double>();
        par_eSPA2A = config["par_eSPA2A"].as<double>();
        par_eSPA1B = config["par_eSPA1B"].as<double>();
        par_eSPA2B = config["par_eSPA2B"].as<double>();
        par_eSPA1V = config["par_eSPA1V"].as<double>();
        par_eSPA2V = config["par_eSPA2V"].as<double>();
        par_eSPA1C = config["par_eSPA1C"].as<double>();
        par_eSPA2C = config["par_eSPA2C"].as<double>();
        par_eSPA1D = config["par_eSPA1D"].as<double>();
        par_eSPA2D = config["par_eSPA2D"].as<double>();

        par_eSPB1B = config["par_eSPB1B"].as<double>();
        par_eSPB2B = config["par_eSPB2B"].as<double>();
        par_eSPB1V = config["par_eSPB1V"].as<double>();
        par_eSPB2V = config["par_eSPB2V"].as<double>();
        par_eSPB1C = config["par_eSPB1C"].as<double>();
        par_eSPB2C = config["par_eSPB2C"].as<double>();
        par_eSPB1D = config["par_eSPB1D"].as<double>();
        par_eSPB2D = config["par_eSPB2D"].as<double>();

        par_eSPC1C = config["par_eSPC1C"].as<double>();
        par_eSPC2C = config["par_eSPC2C"].as<double>();
        par_eSPC1D = config["par_eSPC1D"].as<double>();
        par_eSPC2D = config["par_eSPC2D"].as<double>();
        par_eSPC1V = config["par_eSPC1V"].as<double>();
        par_eSPC2V = config["par_eSPC2V"].as<double>();

        par_eSPD1D = config["par_eSPD1D"].as<double>();
        par_eSPD2D = config["par_eSPD2D"].as<double>();
        par_eSPD1V = config["par_eSPD1V"].as<double>();
        par_eSPD2V = config["par_eSPD2V"].as<double>();

        par_eSPV1V = config["par_eSPV1V"].as<double>();
        par_eSPV2V = config["par_eSPV2V"].as<double>();
      
    };


};












#endif