#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
//#include <chrono>
//#include <random>
#include <cmath>
#include <time.h>
#include <stdlib.h>
#include "matrix.h"
#include "player.h"
#include "config.h"

using namespace player;
using namespace std;
using namespace matrix;
using namespace config;

const double sqrt3over2(sqrt(3.0)/2);
const double pi(4 * atan(1.0));

int factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

double myRandom()
{
    return ((double)rand()/(double)RAND_MAX);
}

void randomGame(Matrix &payoff, Matrix &oppPayoff, const double gamma, const double ep_x, const double ep_y, const double trace)
{
    double u, v, U, V;
    for (int m(1); m <= payoff.columns(); ++m)
    {
        for (int n(1); n <= payoff.rows(); ++n)
        {
            u = myRandom();
            v = myRandom();

            U = sqrt(-2* log(u)) * cos(2 * pi * v);
            V = sqrt(-2* log(u)) * sin(2 * pi * v);

            payoff(n, m) = U;
            oppPayoff(m, n) = (gamma * U + sqrt(1 - gamma*gamma) * V);

            // generate two matrices with correlation determined by G (gamma)
        }
    }
}

Player::Player(const Matrix gamePayoff, Config config) :
    payoff(gamePayoff), N(payoff.rows()), strategy(N, 1), newStrategy(N, 1), perturb(N, 1)
{
    config.set(beta, "beta");
    config.set(lambda, "lambda");
    config.set(tau, "tau");
    config.set(kappa, "kappa");
    config.set(kLevel, "kLevel");
    config.set(label, "label");

    //poissonD.resize(5);
    //kStrategies.resize(5);
    // following Camerer & Ho Poisson CH model, poissonD is a poisson distn. with mean tau.

//    for (int p(0); p < (int)poissonD.size(); ++p)
//    {
//        poissonD[p] = pow(tau, (double)p) * exp(-tau) / (double)factorial(p);
//    }
//
//    Matrix randomStrategy(N, 1);
//    for (int q(1); q <= N; ++q) {
//        randomStrategy[q] = 1 / N;
//    }
//
//    kStrategies[0] = randomStrategy;

    initialise();
}

Player::Player(const Player& parent)
    :   payoff(parent.payoff), N(parent.N), kLevel(parent.kLevel), beta(parent.beta), lambda(parent.lambda), kappa(parent.kappa), label(parent.label),
        strategy(parent.strategy), newStrategy(N, 1), perturb(N, 1)
{}

void Player::initialise()
{
    for (int i(1); i <= N; ++i) {
        strategy[i] = 1. / (double)N;
    }
}

void Player::initialiseRandom()
{
    double normalisation(0);
    for (int i(1); i <= N; ++i)
    {
        strategy[i] = myRandom();
        normalisation += strategy[i];
    }

    strategy = strategy * (1 / normalisation);
}

void Player::initialise(Matrix startPoint)
{
    strategy = startPoint;
}

void Player::match(const Player &opponent)
{
    strategy = opponent.strategy;
}

void Player::duplicate(const Player &target)
{
    strategy = target.strategy;
    payoff = target.payoff;
    lambda = target.lambda;
    beta = target.beta;
}

void Player::interact(Player &opponent) // closed map
{
    opponent.strategyBuffer = opponent.mapStrategy(opponent, *this, opponent.kLevel);
    strategy = mapStrategy(*this, opponent, kLevel);
    opponent.strategy = opponent.strategyBuffer;
}

matrix::Matrix Player::mapStrategy(Player &me, Player &opponent, const int K)
{
    if (K == 1)
    {
        kStrategyBuffer = payoff * opponent.strategy;

        for (int i(1); i <= N; ++i) {
            newStrategy[i] = pow(me.strategy[i], (1 - lambda)) * exp(beta * kStrategyBuffer[i]);
        }
    }

    else
    {
        kBuffer0 = mapStrategy(opponent, me, K - 1);

        if (kappa == 1) // 2nd order kLR map
        {
            kBuffer1 = opponent.payoff * me.strategy;
            kBuffer2 = kBuffer0.transpose() * (opponent.payoff * me.strategy);

            for (int l(1); l <= N; ++l)
            {
                accum = 0;
                for (int m(1); m <= N; ++m) {
                    accum += me.payoff(l, m) * kBuffer0[m] * (1 + beta * (kBuffer1[m] - kBuffer2[1]));
                }

                kStrategyBuffer[l] = accum;
            }
        }

        else kStrategyBuffer = me.payoff * kBuffer0; // 1st order kLR map

        for (int j(1); j <= N; ++j) {
            newStrategy[j] = pow(me.strategy[j], (1 - lambda)) * exp(beta * kStrategyBuffer[j]);
        }
    }

    denom = 0;
    for (int i(1); i <= N; ++i)
    {
        denom += newStrategy[i];
    }

    return newStrategy * (1 / denom);

}

//double Player::componentVariance(Player &opponent, const int component, const int T)
//{
//    int T0((int)(T * 0.666));
//    double x(0), S(0), S2(0), recipT(0), recipT2(0);
//    vector<double> xComponents(N, 0), x2Components(N, 0);
//
//    recipT =  1/( (double)T - (double)T0 );
//    recipT2 = recipT * recipT;
//
//    for (int q(0); q < T0; ++q) {
//        interact(opponent);         // equilibrate
//    }
//
//    for (int q(T0); q < T; ++q)
//    {
//        interact(opponent);
//
//        for (int p(0); p < N; ++p)
//        {
//            x = strategy[p+1];
//            x2Components[p] += x*x;
//            xComponents[p] += x;
//        }
//    }
//
//    S2 = x2Components[component - 1];
//    S = xComponents[component - 1] * xComponents[component - 1];
//
//    return recipT * S2 - recipT2 * S;
//}
//
//double Player::variance(Player &opponent, const int T)
//{
//    int T0(T * 0.8);
//    double x(0), S(0), S2(0), recipT, recipT2, recipN;
//    vector<double> xComponents(N, 0), x2Components(N, 0);
//
//    recipT =  1/( (double)T - (double)T0 );
//    recipT2 = recipT * recipT;
//    recipN = 1 / (double)N;
//
//    for (int q(0); q < T0; ++q) {
//        interact(opponent);         // equilibrate
//    }
//
//    for (int q(T0); q < T; ++q)
//    {
//        interact(opponent);
//
//        for (int p(0); p < N; ++p)
//        {
//            x = strategy[p+1];
//            x2Components[p] += x*x;
//            xComponents[p] += x;
//        }
//    }
//
//    for (int p(0); p < N; ++p)
//    {
//        S2 += x2Components[p];
//        S += xComponents[p] * xComponents[p];
//    }
//
////    cout << " " << S2 << " " << S << endl
//
//    return recipN * (recipT * S2 - recipT2 * S);
//}

void Player::normalise(Matrix &vector, const double &length)
{
    vector = vector * ( length / vector.mod() );
}

double Player::divergence(Player &opponent, const int index, const double pertMax, const int arcLength, const int noPerturbations)
{
    double accumulator, phi(1.), lya, normalisation, normalisation2;

    Player dummy(*this), dummyOpponent(opponent);

    accumulator = 0;
    for (int d(0); d < noPerturbations; ++d)
    {
        // generate perturbation vector
        for (int i(1); i <= N; ++i) {
            perturb[i] = myRandom();
        }
        normalise(perturb, pertMax);

        // if the perturbation takes the strategy outside of the simplex, perturb toward the interior
        for (int i(1); i <= N; ++i)
        {
            if ( (strategy + perturb)[i] <= 0 || (strategy + perturb)[i] >= 1)
            {
                normalisation = 0;
                for (int i(1); i <= N; ++i)
                {
                    perturb[i] = myRandom(); // random point in interior
                    normalisation += perturb[i];
                }
                perturb = perturb * (1 / normalisation);
                perturb = perturb - strategy;
                normalise(perturb, pertMax);

                break;
            }
        }

        // perturb strategy
        dummy.strategy = strategy + perturb;

        // likewise for the opponent
            for (int i(1); i <= N; ++i) {
                perturb[i] = myRandom();
            }
            normalise(perturb, pertMax);

            for (int i(1); i <= N; ++i)
            {
                if ( (opponent.strategy + perturb)[i] <= 0 || (opponent.strategy + perturb)[i] >= 1)
                {
                    normalisation = 0;
                    for (int i(1); i <= N; ++i)
                    {
                        perturb[i] = myRandom(); // random point in interior
                        normalisation += perturb[i];
                    }
                    perturb = perturb * (1 / normalisation);
                    perturb = perturb - opponent.strategy;
                    normalise(perturb, pertMax);

                    break;
                }
            }

            dummyOpponent.strategy = opponent.strategy + perturb;

        // normalise the perturbed strategy
            normalisation = 0;
            normalisation2 = 0;

            for (int i(1); i <= N; ++i)
            {
                normalisation += dummy.strategy[i];
                normalisation2 += dummyOpponent.strategy[i];
            }

            dummy.strategy = dummy.strategy * (1 / normalisation);
            dummyOpponent.strategy = dummyOpponent.strategy  * (1 / normalisation);
        //

        // run perturbed trajectory
        for (int q(0); q < arcLength; ++q)
        {
            interact(opponent);
            dummy.interact(dummyOpponent);

            // divergence vectors
            perturb = dummy.strategy - strategy;
            perturbOpp = dummyOpponent.strategy - opponent.strategy;

            lya = perturb.mod() / pertMax;

            if (lya == lya && q > 0) phi *= lya; // excludes nan and drops first value

            normalise(perturb, pertMax);
            normalise(perturbOpp, pertMax);

            dummy.strategy = strategy + perturb;
            dummyOpponent.strategy = opponent.strategy + perturbOpp;
        }

        accumulator += log(phi); // divergence
    }

    return accumulator / ((double)noPerturbations * (double)arcLength);
}

//void Player::plotDivergence(Player &opponent, const int index, const double pertMax, const int arcLength, const int noPerturbations, const int noPoints)
//{
//    double accumulator(0), normalisation(0), normalisation2(0), phi(1.), lya;
//
//    stringstream ssFilename;
//    string sFilename;
//
//    sFilename = "";
//    ssFilename.str("");
//    ssFilename << label << kLevel << opponent.kLevel << "_divA" << index << ".dat";
//    sFilename = ssFilename.str();
//
//    ofstream foutPlayer(sFilename.c_str());
//
//    sFilename = "";
//    ssFilename.str("");
//    ssFilename << label << kLevel << opponent.kLevel << "_divB" << index << ".dat";
//    sFilename = ssFilename.str();
//
//    ofstream foutDummy(sFilename.c_str());
//
//    Player dummy(*this), dummyOpponent(opponent);
//
//    initialise();
//    opponent.match(*this);
//
//    cout << "\nRunning k(" << kLevel << ", " << opponent.kLevel << ")" << endl;
//    cout << "\nb_1 = " << beta << ", l_2 = " << lambda << endl << "b_2 = " << opponent.beta << ", l_2 = " << opponent.lambda << endl;
//
//    //int T0 = (int)((double)noPoints * 0.66);
//    for (int q(0); q < noPoints; ++q) // equilibrate
//    {
//        interact(opponent);
//        foutPlayer << *this << endl;
//    }
//
//    accumulator = 0;
//
//    for (int d(0); d < noPerturbations; ++d)
//    {
//        foutDummy << endl;
//
//        // generate perturbation vector
//        for (int i(1); i <= N; ++i) {
//            perturb[i] = myRandom();
//        }
//        normalise(perturb, pertMax);
//
//        // if the perturbation takes the strategy outside of the simplex, perturb toward the interior
//        for (int i(1); i <= N; ++i)
//        {
//            if ( (strategy + perturb)[i] <= 0 || (strategy + perturb)[i] >= 1)
//            {
//                normalisation = 0;
//                for (int i(1); i <= N; ++i)
//                {
//                    perturb[i] = myRandom(); // random point in interior
//                    normalisation += perturb[i];
//                }
//                perturb = perturb * (1 / normalisation);
//                perturb = perturb - strategy;
//                normalise(perturb, pertMax);
//
//                break;
//            }
//        }
//
//        dummy.strategy = strategy + perturb;
//
//        for (int i(1); i <= N; ++i) {
//            perturb[i] = myRandom();
//        }
//        normalise(perturb, pertMax);
//
//        for (int i(1); i <= N; ++i)
//        { // likewise for the opponent
//            if ( (opponent.strategy + perturb)[i] <= 0 || (opponent.strategy + perturb)[i] >= 1)
//            {
//                normalisation = 0;
//                for (int i(1); i <= N; ++i)
//                {
//                    perturb[i] = myRandom(); // random point in interior
//                    normalisation += perturb[i];
//                }
//                perturb = perturb * (1 / normalisation);
//                perturb = perturb - opponent.strategy;
//                normalise(perturb, pertMax);
//
//                break;
//            }
//        }
//
//        dummyOpponent.strategy = opponent.strategy + perturb;
//
//        // try to normalise the perturbed strategy as best you can
//        normalisation = 0;
//        normalisation2 = 0;
//        for (int i(1); i <= N; ++i)
//        {
//            normalisation += dummy.strategy[i];
//            normalisation2 += dummyOpponent.strategy[i];
//        }
//        dummy.strategy = dummy.strategy * (1 / normalisation);
//        dummyOpponent.strategy = dummyOpponent.strategy  * (1 / normalisation);
//        //
//
//        interact(opponent);
//        dummy.interact(dummyOpponent);
//
//        foutPlayer << *this << endl;
//        foutDummy << "1" << endl << endl;
//
//        perturb = dummy.strategy - strategy;
//        perturbOpp = dummyOpponent.strategy - opponent.strategy;
//
//        normalise(perturb, pertMax);
//        normalise(perturbOpp, pertMax);
//
//        dummy.strategy = strategy + perturb;
//        dummyOpponent.strategy = opponent.strategy + perturbOpp;
//
//        // run perturbed trajectory
//        for (int q(1); q < arcLength; ++q)
//        {
//            interact(opponent);
//            dummy.interact(dummyOpponent);
//
//            foutPlayer << *this << endl;
//            foutDummy << dummy << endl;
//
//            perturb = dummy.strategy - strategy;    // divergence vectors
//            perturbOpp = dummyOpponent.strategy - opponent.strategy;
//
//            lya = perturb.mod() / pertMax;
//
//            if (lya == lya) phi *= lya;
//
//            normalise(perturb, pertMax);
//            normalise(perturbOpp, pertMax);
//
//            dummy.strategy = strategy + perturb;
//            dummyOpponent.strategy = opponent.strategy + perturbOpp;
//        }
//
////        perturb = dummy.strategy - strategy;    // divergence vectors
////        perturbOpp = dummyOpponent.strategy - opponent.strategy;
////
////        lya = perturb.mod() / pertMax;
//
//        accumulator += log(phi); // divergence
//    }
//
//    cout << "\nLargest Lyapunov exponent = " << accumulator / ((double)noPerturbations * (double)arcLength); // average Lyapunov over phase space evolution
//}


void Player::mapDivergence(Player &opponent, Config &configSys, const int batch_x, const int batch_y, const int batchSize)
{

    int counter(0);
    double accumulator, accumulator1, accumulator2, gamma, pertMax, xMin, xMax, yMin, yMax, ep_x(0), ep_y(0), trace(1), dummy;
    double temp0, temp1, temp2;
    double *x1, *x2, *y1, *y2;
    int dataPoints, noPerturbations, arcLength, noSamples, xIntervals, yIntervals;
    string xParameter, yParameter, gameType;

    Player control1(*this), control1Opponent(opponent), control2(*this), control2Opponent(opponent);
    //control1.kLevel = 1;
    control1Opponent.kLevel = 1;
    control2Opponent.kLevel = kLevel - 1;

    configSys.set(dataPoints, "dataPoints");
    configSys.set(noPerturbations, "noPerturbations");
    configSys.set(arcLength, "arcLength");
    configSys.set(noSamples, "noSamples");
    configSys.set(pertMax, "pertMax");
    configSys.set(gameType, "gameType");
    configSys.set(xParameter, "xParameter");
    configSys.set(yParameter, "yParameter");
    configSys.set(xIntervals, "xIntervals");
    configSys.set(yIntervals, "yIntervals");
    configSys.set(gamma, "gamma");
    configSys.set(trace, "trace");

    stringstream ssFilename;
    string sFilename;

    sFilename = "";
    ssFilename.str("");

    ssFilename << gameType << "_";

    if (xParameter == "lambda")
    {
        configSys.set(xMin, "lMin");
        configSys.set(xMax, "lMax");
        x1 = &lambda;
        x2 = &opponent.lambda;

        ssFilename << "L";
    }

    else if (xParameter == "ep_x")
    {
        configSys.set(xMin, "eMin");
        configSys.set(xMax, "eMax");
        x1 = &ep_x;
        x2 = &dummy;

        ssFilename << "eX";
    }

    else
    {
        x1 = &dummy;
        x2 = &dummy;
        throw ("Invalid plotting parameter for x");
    }

    if (yParameter == "gamma")
    {
        configSys.set(yMin, "gMin");
        configSys.set(yMax, "gMax");
        y1 = &gamma;
        y2 = &dummy;

        ssFilename << "G";
    }

    else if (yParameter == "ep_y")
    {
        configSys.set(yMin, "eMin");
        configSys.set(yMax, "eMax");
        y1 = &ep_y;
        y2 = &dummy;

        ssFilename << "eY";
    }

    else
    {
        y1 = &dummy;
        y2 = &dummy;
        throw ("Invalid plotting parameter for y");
    }

    ssFilename << "_lyap_k" << kLevel << opponent.kLevel << "_";
        if (batch_x < 10) ssFilename << "0"; // for proper filename formatting
    ssFilename << batch_x << "_";
        if (batch_y < 10) ssFilename << "0"; // for proper filename formatting
    ssFilename << batch_y;
    sFilename = ssFilename.str() + ".dat";

    cout << "\nWriting to " << sFilename << endl;

    ofstream foutPlayer(sFilename.c_str());

    ssFilename << "all";
    sFilename = ssFilename.str() + ".dat";

    ofstream foutPlayerAll(sFilename.c_str());

    // batch processing modifications
    double dx( (xMax - xMin)/(double)xIntervals ), dy( (yMax - yMin)/(double)yIntervals );

    xIntervals /= batchSize;
    yIntervals /= batchSize;

    xMin += dx * xIntervals * batch_x;
    xMax = xMin + dx * xIntervals;

    yMin += dy * yIntervals * batch_y;
    yMax = yMin + dy * yIntervals;

    cout << "Calculating LLE over (" << noPerturbations << "|" << arcLength << ") perturbations after " << (int)(dataPoints * 0.66) << " rounds." << endl;

    *y1 = yMin;
    *y2 = yMin;

    for (int y(0); y < yIntervals; ++y)
    {
        *y1 += dy;
        *y2 += dy;

        *x1 = xMin;
        *x2 = xMin;

        for (int x(0); x < yIntervals; ++x)
        {
            *x1 += dx;
            *x2 += dx;
            counter++;

            accumulator = 0;
            accumulator1 = 0;
            accumulator2 = 0;

            for (int p(0); p < noSamples; ++p)
            {
                cout << "[" << p << "/" << counter << "/" << xIntervals*yIntervals << "] " << yParameter << " = " << *y1 << ", " << xParameter << " = " << *x1 << endl;

                // generate random game if required
                if (gameType == "large") randomGame(payoff, opponent.payoff, gamma, ep_x, ep_y, trace);

                // add diagonal term
                for (int t(1); t <= N; ++t)
                {
                    payoff(t, t) = payoff(t, t) * trace + ep_x;
                    opponent.payoff(t, t) = opponent.payoff(t, t) * trace + ep_y;
                }

                // don't initialise at the fixed point for RPS, otherwise start in the centre.
                if (gameType == "rps") initialiseRandom();
                else initialise();

                opponent.match(*this);

                control1.duplicate(*this);
                control2.duplicate(*this);

                control1Opponent.duplicate(opponent);
                control2Opponent.duplicate(opponent);

                // equilibrate
                for (int q(0); q < dataPoints * 0.666; ++q)
                {
                    interact(opponent);
                    control1.interact(control1Opponent);
                    control2.interact(control2Opponent);
                }

                // calculate trajectory divergence

                temp0 = divergence(opponent, 0, pertMax, arcLength, noPerturbations); // A vs B
                accumulator += temp0;

                temp1 = control1.divergence(control1Opponent, 0, pertMax, arcLength, noPerturbations);
                accumulator1 += temp1;

                temp2 = control2.divergence(control2Opponent, 0, pertMax, arcLength, noPerturbations);
                accumulator2 += temp2;

                foutPlayerAll << *x1 << " " << *y1 << " " << temp0 << " " << temp2 << " " << temp1 << endl;

            }

            // output averaged values
            foutPlayer << *x1 << " " << *y1 << " " << accumulator / (double)noSamples << " " << accumulator2 / (double)noSamples << " " << accumulator1 / (double)noSamples << endl;
        }
    }
}

//// varying beta and lambda
//void Player::mapVarianceBL(Player &opponent, Matrix startpoint, const int index, const int noSamples, const int xIntervals, const int yIntervals, const double xMin, const double xMax, const double yMin, const double yMax, const int T)
//{
//    stringstream ssFilename;
//    string sFilename;
//
//    sFilename = "";
//    ssFilename.str("");
//    ssFilename << label << "_variance" << index << ".dat";
//    sFilename = ssFilename.str();
//
//    ofstream foutPlayer(sFilename.c_str());
//
//    double xInterval( (xMax - xMin)/(double)xIntervals ), yInterval( (yMax - yMin)/(double)yIntervals );
//
//    cout << "Calculating variance between rounds " << T * 0.8 << " - " << T << endl;
//
//    for (double x(xMin + xInterval); x <= xMax; x += xInterval)
//    {
//        for (double y(yMin + yInterval); y <= yMax; y += yInterval)
//        {
//            beta = x;
//            opponent.beta = x;
//
//            lambda = y;
//
//            opponent.lambda = y;
//
//            cout << "beta = " << beta << ", lambda = " << lambda;
//
//            initialise(startpoint);
//            opponent.match(*this);
//
//            foutPlayer << beta << " " << lambda << " " << variance(opponent, T);
//
//            for (int p(1); p <= N; ++p) {
//                foutPlayer << " " << componentVariance(opponent, p, T); // individual component variance
//            }
//
//            for (int s(1); s < noSamples; ++s)
//            {
//                initialiseRandom(); // random initial conditions
//                opponent.match(*this);
//
//                foutPlayer << " " << variance(opponent, T); // average variance
//
//                for (int p(1); p <= N; ++p) {
//                    foutPlayer << " " << componentVariance(opponent, p, T); // individual component variance
//                }
//            }
//
//            foutPlayer << endl;
//            cout << " ... done" << endl;
//        }
//    }
//}
//
//// varying lambda and gamma
//void Player::mapVarianceGL(Player &opponent, Matrix startpoint, const int index, const int noSamples, const int xIntervals, const int yIntervals, const double xMin, const double xMax, const double yMin, const double yMax, const int T)
//{
//    stringstream ssFilename;
//    string sFilename;
//
//    sFilename = "";
//    ssFilename.str("");
//    ssFilename << label << "_varianceB" << index << ".dat";
//    sFilename = ssFilename.str();
//
//    ofstream foutPlayer(sFilename.c_str());
//
//    double xInterval( (xMax - xMin)/(double)xIntervals ), yInterval( (yMax - yMin)/(double)yIntervals );
//    double gamma, accumulator, counter(0);
//
//    cout << "Calculating variance between rounds " << T * 0.666 << " - " << T << endl;
//
//    for (double x(xMin + xInterval); x < xMax; x += xInterval)
//    {
//        for (double y(yMin + yInterval); y < yMax; y += yInterval)
//        {
//            gamma = x;
//
//            lambda = y;
//            alpha = (1 - y);
//
//            opponent.lambda = y;
//            opponent.alpha = (1 - y);
//
//            foutPlayer << lambda << " " << gamma;
//
//            accumulator = 0;
//
//            for (int s(0); s < noSamples; ++s)
//            {
//                counter++;
//
//                initialise(); //  initial conditions
//                opponent.match(*this);
//
//                cout << fixed << "gamma = " << gamma << ", lambda = " << lambda << " " << "[" << counter << "/" << noSamples*xIntervals*yIntervals << "]" << endl;
//
//                randomGame(payoff, opponent.payoff, gamma, 0, 0, 1);
//
//                accumulator += variance(opponent, T);
//            }
//
//            foutPlayer << " " << accumulator / (double)noSamples;
//            foutPlayer << endl;
//
//        }
//    }
//}

void Player::plotTrajectories(Player &opponent, Matrix startpointA, Matrix startpointB, const int index, const int noTrajectories, const int noPoints)
{
    stringstream ssFilename;
    string sFilename;

    sFilename = "";
    ssFilename.str("");
    ssFilename << label << "_traj" << index << ".dat";
    sFilename = ssFilename.str();

    ofstream foutPlayer(sFilename.c_str());

    sFilename = "";
    ssFilename.str("");
    ssFilename << opponent.label << "_traj" << index << ".dat";
    sFilename = ssFilename.str();

    ofstream foutOpponent(sFilename.c_str());

    for (int p(0); p < noTrajectories; ++p) // plot m trajectories
    {
        cout << "\nPlotting trajectory " << p << "; lambda=" << lambda <<  " beta=" << beta << endl;

        initialise(startpointA);
        opponent.initialise(startpointB);

        for (int q(0); q < noPoints; ++q) // plot q data points per trajectory
        {
            interact(opponent);

            foutPlayer << fixed << setprecision(8) << *this << endl;
            foutOpponent << fixed << setprecision(8) << opponent << endl;
        }
    }
}

void Player::plotTrajectories(Player &opponent, const int index, const int noTrajectories, const int noPoints)
{
    stringstream ssFilename;
    string sFilename;

    sFilename = "";
    ssFilename.str("");
    ssFilename << label << "_traj" << index << ".dat";
    sFilename = ssFilename.str();

    ofstream foutPlayer(sFilename.c_str());

    sFilename = "";
    ssFilename.str("");
    ssFilename << opponent.label << "_traj" << index << ".dat";
    sFilename = ssFilename.str();

    ofstream foutOpponent(sFilename.c_str());

    for (int p(0); p < noTrajectories; ++p) // plot m trajectories
    {
        cout << "\nPlotting trajectory " << p << "; lambda=" << lambda <<  " beta=" << beta << endl;

        initialiseRandom();
        opponent.match(*this);

        for (int q(0); q < noPoints; ++q) // plot q data points per trajectory
        {
            interact(opponent);

            foutPlayer << fixed << setprecision(8) << *this << endl;
            foutOpponent << fixed << setprecision(8) << opponent << endl;
        }
    }
}

void Player::plotVectorField(const int intervals)
{
    stringstream ssFilename;
    string sFilename;

    ssFilename.str("");
    ssFilename << label << "_vector.dat";
    sFilename = ssFilename.str();

    ofstream fout(sFilename.c_str());

    double  x(0), y(0), z(0);
    pair<double,double> ternary;
    Matrix vector(N);

    for (int n(0); n < intervals; ++n)
    {
        for (int m(0); m < intervals ; ++m)
        {
            x = (double)n / (double)intervals;
            y = (double)m / (double)intervals;
            z = 1 - x - y;

            if (z >= 0)
            {
                strategy[1] = x;
                strategy[2] = y;
                strategy[3] = z;

                ternary = tern(z, x);

                fout << fixed << setprecision(8) << ternary.first << " " << ternary.second;

                fout << "\t";
                vector = mapStrategy(*this, *this, kLevel) - strategy;
                ternary = tern(vector[3], vector[1]);

                fout << ternary.first << " " << ternary.second << endl;
            }
        }
    }
}

pair<double, double> Player::tern(double a, double b) const
{
	pair <double, double> ternary;

	ternary.first = 0.5*(2*a+b);
	ternary.second = sqrt3over2 * b;

	return ternary;
}

ostream &player::operator<<(ostream &os, const Player &player)
{
    for (int i(1); i <=player.N; ++i) {
        os << player.strategy[i] << " ";
    }

    pair<double,double> ternary = player.tern(player.strategy[3], player.strategy[1]);
    os << "\t" << ternary.first << " " << ternary.second << "\t";

    return os;
}
