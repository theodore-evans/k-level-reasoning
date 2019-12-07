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

int main(int argc, char* argv[])
{
    try
    {
        srand (time(NULL));
//        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
//        normal_distribution<double> normal(0, 1);
//        ranlux64_base_01 generator(seed);
        // initialise random number generators

        int N, dataPoints, xIntervals, yIntervals, noSamples, batch_x(0), batch_y(0), batchSize(1);

        Config configSys("system.ini"), configA("alice.ini"), configB("bob.ini");

        double  R, S, T, P, c, m, x0, y0, z0, x20, y20, z20, e_x, e_y, gamma;
        double gMin, gMax, bMin, bMax, lMin, lMax, beta;

        string gameType;

        configSys.set(gameType, "gameType");

        configSys.set(R, "R");
        configSys.set(S, "S");
        configSys.set(T, "T");
        configSys.set(P, "P");
        configSys.set(m, "m");
        configSys.set(c, "c");
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

        Matrix prisoners(3,3), rps(3,3), I3(3,3), rps_x(3,3), rps_y(3,3), largeA(N, N), largeB(N, N), large(N, N);
        I3(1,1) = 1; I3(2,2) = 1; I3(3,3) = 1;

        // IMHOF IPD //
                    //AllC              //ALLD              //TFT
        /*ALLC*/    prisoners(1,1) = R;    prisoners(1,2) = S;    prisoners(1,3) = R;
        /*ALLD*/    prisoners(2,1) = T;    prisoners(2,2) = P;    prisoners(2,3) = (T + P * (m - 1)) / m;
        /*TFT*/     prisoners(3,1) = R  - c/m; prisoners(3,2) = (S + P *(m - 1) - c)/m; prisoners(3,3) = R  - c/m;

        // RPS //

        /*rock*/        rps(1,1) = 0;    rps(1,2) = -1;    rps(1,3) = 1;
        /*paper*/       rps(2,1) = 1;    rps(2,2) = 0;    rps(2,3) = -1;
        /*scissors*/    rps(3,1) = -1;  rps(3,2) = 1;   rps(3,3) = 0;

        // RPS BROKEN SYMMETRY (Sato et Al. 2001) //

        rps_x = rps + I3 * e_x;
        rps_y = rps + I3 * e_y;

        // RANDOM LARGE GAMES //

        ofstream matrix("matrix.dat");

        Matrix startpointA(N), startpointB(N);

        startpointA[1] = log(10*x0);
        startpointA[2] = log(10*y0);
        startpointA[3] = log(10*z0);

        startpointB[1] = log(10*x20);
        startpointB[2] = log(10*y20);
        startpointB[3] = log(10*z20);

        // instantiate players

        std::cout << "Batch (" << batch_x+1 << "/" << batch_y+1 << ") out of (" << batchSize << "/" << batchSize << ")" << endl;
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
        std::cout << ", k(" << alice.kLevel << ", " << bob.kLevel << ")" << endl;

        alice.mapDivergence(bob, configSys, batch_x, batch_y, batchSize);


//        alice.plotDivergence(bob, batch_x, pertMax, arcLength, noPerturbations, dataPoints);
//        std::cout << endl;
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

    catch(const string &err)
    {
        std::cout << "\nException: " << err << endl;
    }

    catch(exception &ex)
    {
        std::cout << "\nException: " << ex.what() << endl;
    }
}
