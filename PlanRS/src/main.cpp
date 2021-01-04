#include "ReedsSheppPlanner.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ompl/base/spaces/ReedsSheppStateSpace.h>
#include <ompl/base/ScopedState.h>

int main(int argc, char *argv[])
{
    float turningRad = 1.0; //default value

    float startX = 2.5;
    float startY = 0.01;
    float startTh = 0.5 * boost::math::constants::pi<double>();

    float goalX = 7.5;
    float goalY = 0.01;
    float goalTh = -0.5 * boost::math::constants::pi<double>();

    float xHigh = 20.0;
    float yHigh = 10.0;
    float robotWidth = 0.0;
    float robotLength = 0.0;

    int easyCollision = 1;
    float solveTime = 30.0;
    float stateValidityCheckingRes = 0.005;
    int interpolateNum = 1000;

    if (argc > 1) //read in turning radius
    {
        turningRad = atof(argv[1]);
    }
    else if (argc > 2)
    {
        interpolateNum = atof(argv[2]);
    }
    
    ReedsShepp::RSPlanner RSPlanner(turningRad, startX, startY, startTh, goalX, goalY, goalTh, xHigh, yHigh, robotWidth, robotLength);

    RSPlanner.planRS(easyCollision, solveTime, stateValidityCheckingRes, interpolateNum);

    return 0;
}