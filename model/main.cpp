#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <ctime>
//#include <random>
//#include <chrono>
#include <stdlib.h>
#include <json/json.h>

#include "matrix.h"
#include "player.h"
#include "config.h"

using namespace matrix;
using namespace player;
using namespace config;

// command line arguments <number of batch intervals 'per side'> <batch x index, starting at 0> <batch y index, starting at 0>
// eg. a single batch (equivalent to no argument default) would take the arguments '1 0 0'

const double pi(4 * atan(1.0));

double myRandom();
void randomGame(Matrix &payoff, Matrix &oppPayoff, const double gamma, const double ep_x, const double ep_y, const double trace);

// definitions for normal form game payoff matrices
// TODO? add these to their own source file
matrix::Matrix ipd_payoff(Json::Value config)
{
  double t = config["games"]["ipd"]["t"]
  double r = config["games"]["ipd"]["r"]
  double p = config["games"]["ipd"]["p"]
  double s = config["games"]["ipd"]["s"]
  double m = config["games"]["ipd"]["m"]
  double c = config["games"]["ipd"]["c"]

          //AllC                    //ALLD                //TFT
/*ALLC*/   payoff(1,1) = r;         payoff(1,2) = s;      payoff(1,3) = r;
/*ALLD*/   payoff(2,1) = t;         payoff(2,2) = p;      payoff(2,3) = (t + p * (m - 1)) / m;
/*TFT*/    payoff(3,1) = r  - c/m;  payoff(3,2) = (s + p *(m - 1) - c)/m;  payoff(3,3) = r  - c/m;

  return payoff;
}

matrix::Matrix rps_payoff()
{
/*rock*/        rps(1,1) = 0;    rps(1,2) = -1;  rps(1,3) = 1;
/*paper*/       rps(2,1) = 1;    rps(2,2) = 0;   rps(2,3) = -1;
/*scissors*/    rps(3,1) = -1;   rps(3,2) = 1;   rps(3,3) = 0;

  return payoff;
}

int main(int argc, char* argv[])
{
    try
    {
        // initialise random number generator
        srand (time(NULL))

        //parse configuration json file
        std::ifstream file_input("config.json");
        Json::Reader reader;
        Json::Value config;
        reader.parse(file_input, config);

        // TODO phase out the homemade config file parser in favour of a json-based system
        // TODO remove these value declarations and declare at initialisation instead
        int N, dataPoints, xIntervals, yIntervals, noSamples, batch_x(0), batch_y(0), batchSize(1);
        double  x0, y0, z0, x20, y20, z20, e_x, e_y, gamma;
        double gMin, gMax, bMin, bMax, lMin, lMax, beta;

        std::string gameType;

        configSys.set(gameType, "gameType");

        configSys.set(gamma, "gamma");
        configSys.set(beta, "beta");

        configSys.set(gMin, "gMin");
        configSys.set(gMax, "gMax");
        configSys.set(bMin, "bMin");
        configSys.set(bMax, "bMax");
        configSys.set(lMin, "lMin");
        configSys.set(lMax, "lMax");
        configSys.set(noSamples, "noSamples");
        configSys.set(xIntervals, "xIntervals");
        configSys.set(yIntervals, "yIntervals");

        if (argc == 4) // if command line arguments are given for batch processing, adjust limits accordingly
        {
            batchSize = atoi(argv[1]);
            batch_x = atoi(argv[2]);
            batch_y = atoi(argv[3]);
        }

        configSys.set(dataPoints, "dataPoints");
        configSys.set(N, "N");

        configA.set(e_x, "ep");
        configB.set(e_y, "ep");

        configA.set(x0, "x0");
        configA.set(y0, "y0");
        configA.set(z0, "z0");

        configB.set(x20, "x0");
        configB.set(y20, "y0");
        configB.set(z20, "z0");

        // Initialise payoff matrices
        Matrix ipd(3,3), rps(3,3), I3(3,3), rps_x(3,3), rps_y(3,3), largeA(N, N), largeB(N, N), large(N, N);
        I3(1,1) = 1; I3(2,2) = 1; I3(3,3) = 1;

        // RPS BROKEN SYMMETRY (Sato et Al. 2001) //

        rps_x = rps + I3 * e_x;
        rps_y = rps + I3 * e_y;

        // RANDOM LARGE GAMES //

        std::ofstream matrix("matrix.dat");

        Matrix startpointA(N), startpointB(N);

        startpointA[1] = log(10*x0);
        startpointA[2] = log(10*y0);
        startpointA[3] = log(10*z0);

        startpointB[1] = log(10*x20);
        startpointB[2] = log(10*y20);
        startpointB[3] = log(10*z20);

        // instantiate players

        std::cout << "Batch (" << batch_x+1 << "/" << batch_y+1 << ") out of (" << batchSize << "/" << batchSize << ")" << std::endl;
        std::cout << "\nGametype = " << gameType << ", N = " << N;

        Matrix gameA, gameB;

        if (gameType == "rps") {
            gameA = rps;
            gameB = rps;
            }

        else if (gameType == "rps_b") {
            gameA = rps_x;
            gameB = rps_y;
            }

        else if (gameType == "large") {
            randomGame(largeA, largeB, gamma, 0, 0, 1);
            gameA = largeA;
            gameB = largeB;
            }

        else if (gameType == "brokenSymmetry") {
            randomGame(largeA, largeB, -1, e_x, e_y, 1);
            gameA = largeA;
            gameB = largeB;
            }

        else if (gameType == "brokenSymmetryTL") {
            randomGame(largeA, largeB, -1, e_x, e_y, 0);
            gameA = largeA;
            gameB = largeB;
            }

        Player alice(gameA, configA), bob(gameB, configB);
        std::cout << ", k(" << alice.kLevel << ", " << bob.kLevel << ")" << std::endl;

        alice.mapDivergence(bob, configSys, batch_x, batch_y, batchSize);


//        alice.plotDivergence(bob, batch_x, pertMax, arcLength, noPerturbations, dataPoints);
//        std::cout << std::endl;
//            alice.kLevel = (alice.kLevel % 2) + 1;
//            bob.kLevel = (bob.kLevel % 2) + 1;
//            alice.plotDivergence(bob, batch_x, pertMax, arcLength, noPerturbations, dataPoints);

        //(P &opponent, i index, d pertMax, i arcLength, i noPerturbations, i noPoints)

        //alice.mapDivergenceBL(bob, batch_x, batch_y, pertMax, arcLength, noPerturbations, xIntervals, yIntervals, bMin, bMax, lMin, lMax, dataPoints);

//        P &opponent, i index, i indey, d pertMax, i arcLength, i noPerturbations, i xIntervals, i yIntervals, i noSamples, d xMin, d xMax, d yMin, d yMax, i T);

//        alice.plotTrajectories(bob, 0, 3, dataPoints);
//        alice.plotVectorField(xIntervals);

        //alice.mapVarianceBL(bob, startpointA, 0, noSamples, noIntervals, bMin, bMax, lMin, lMax, dataPoints);
        //(P &opponent, M startpoint, i index, i noSamples, i noIntervals, d betaMin, d betaMax, d lambdaMin, d lambdaMax, i T)

        //alice.mapVarianceGL(bob, startpointA, 0, noSamples, noIntervals, GMin, GMax, lMin, lMax, dataPoints);
        //(P &opponent, M startpoint, i index, i noSamples, i noIntervals, d gammaMin, d gammaMax, d lambdaMin, d lambdaMax, i T)

        return 0;
    }

    catch(const std::string &err)
    {
        std::cout << "\nException: " << err << std::endl;
    }

    catch(std::exception &ex)
    {
        std::cout << "\nException: " << ex.what() << std::endl;
    }
}
