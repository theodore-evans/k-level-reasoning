#ifndef PLAYER_H_INCLUDED
#define PLAYER_H_INCLUDED

#include "matrix.h"
#include "config.h"
#include <vector>

namespace player
{
    class Player
    {
        friend std::ostream & operator<<(std::ostream &os, const Player &player);

        protected:

        public:
            matrix::Matrix payoff;

            int N, kLevel;
            double beta, lambda, tau, kappa, denom, accum;
            std::string label;
            Player(const matrix::Matrix payoff, config::Config config);
            Player(const Player& parent);
            ~Player() {}

			matrix::Matrix strategy, newStrategy, kBuffer0, kBuffer1, kBuffer2, kStrategyBuffer, strategyBuffer, perturb, perturbOpp;
//			std::vector<double> poissonD;
//			std::vector<matrix::Matrix> kStrategies;

            void initialise();
            void initialiseRandom();
            void initialise(matrix::Matrix startPoint);
            void match(const Player &opponent);
            void duplicate(const Player &target);

            matrix::Matrix mapStrategy(player::Player &me, player::Player &opponent, const int K);

            void interact(Player &opponent);

//            double variance(Player &opponent, const int T);
//            double componentVariance(Player &opponent, const int component, const int T);

            void normalise(matrix::Matrix &vector, const double &length);
            double divergence(Player &opponent, const int index, const double pertMax, const int arcLength, const int noPerturbations);

//            void plotDivergence(Player &opponent, const int index, const double pertMax, const int arcLength, const int noPerturbations, const int noPoints);
//            void mapDivergenceBL(Player &opponent, const int index, const int indey, const double pertMax, const int arcLength, const int noPerturbations, const int xIntervals, const int yIntervals, const double xMin, const double xMax, const double yMin, const double yMax, const int T);

            void mapDivergence(Player &opponent, config::Config &configSys, const int batch_x, const int batch_y, const int batchSize);

//            void mapVarianceBL(Player &opponent, matrix::Matrix startpoint, const int index, const int noSamples, const int xIntervals, const int yIntervals, const double xMin, const double xMax, const double yMin, const double yMax, const int T);
//            void mapVarianceGL(Player &opponent, matrix::Matrix startpoint, const int index, const int noSamples, const int xIntervals, const int yIntervals, const double xMin, const double xMax, const double yMin, const double yMax, const int T);

            void plotVectorField(const int intervals);

            void plotTrajectories(Player &opponent, matrix::Matrix startpointA, matrix::Matrix startpointB, const int index, const int noTrajectories, const int noPoints);
            void plotTrajectories(Player &opponent, const int index, const int noTrajectories, const int noPoints);

            std::pair<double, double> tern(double a, double b) const;
    };

    std::ostream & operator<<(std::ostream &os, const Player &player);
}

#endif

