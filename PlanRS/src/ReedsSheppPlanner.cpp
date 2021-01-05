#include "ReedsSheppPlanner.h"
#include <boost/program_options.hpp>
#include <ompl/base/objectives/PathLengthOptimizationObjective.h>
#include <ompl/geometric/planners/rrt/InformedRRTstar.h>
#include <ompl/geometric/planners/rrt/RRTstar.h>
#include <ompl/geometric/planners/rrt/RRTConnect.h>
#include <ompl/base/PlannerData.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>

namespace ob = ompl::base;
namespace og = ompl::geometric;
namespace po = boost::program_options;

namespace ReedsShepp
{
    RSPlanner::RSPlanner(float arg_turningRadius, float arg_startX, float arg_startY, float arg_startTh, float arg_goalX, float arg_goalY, float arg_goalTh, float arg_xHigh, float arg_yHigh, float arg_robotWidth, float arg_robotLength) : m_turningRadius(arg_turningRadius), m_startX(arg_startX) , m_startY(arg_startY), m_startTh(arg_startTh), m_goalX(arg_goalX), m_goalY(arg_goalY), m_goalTh(arg_goalTh), m_xHigh(arg_xHigh), m_yHigh(arg_yHigh), m_robotWidth(arg_robotWidth), m_robotLength(arg_robotLength) {}

    bool isStateValidEasy(ob::SpaceInformationPtr si, const ob::State *state, RSPlanner thisRSPlanner)
    {
        const auto *s = state->as<ob::SE2StateSpace::StateType>();
        double x = s->getX(), y = s->getY();

        float minDist = fmin(fmin(x, thisRSPlanner.m_xHigh - x), fmin(y, thisRSPlanner.m_yHigh - y));
        float maxLen = fmax(thisRSPlanner.m_robotLength, thisRSPlanner.m_robotWidth / 2.0);

        bool clearance = minDist > maxLen;
        return si->satisfiesBounds(s); // && clearance;
    }

    bool isStateValidHard(ob::SpaceInformationPtr si, const ob::State *state, RSPlanner thisRSPlanner)
    {
        return 1;
    }

    void RSPlanner::planRS(int easy, float solveTime, float stateValidityCheckingRes, int interpolateNum)
    {
        ompl::base::StateSpacePtr space(std::make_shared<ompl::base::ReedsSheppStateSpace>(this->m_turningRadius));

        //set Goal States
        ob::ScopedState<> start(space), goal(space);
        start[0] = this->m_startX;
        start[1] = this->m_startY;
        start[2] = this->m_startTh;

        goal[0] = this->m_goalX;
        goal[1] = this->m_goalY;
        goal[2] = this->m_goalTh;

        ob::RealVectorBounds bounds(2);
        bounds.setLow(0);
        bounds.high[0] = this->m_xHigh;
        bounds.high[1] = this->m_yHigh;
        space->as<ob::SE2StateSpace>()->setBounds(bounds);

        ob::SpaceInformationPtr si(new ob::SpaceInformation(space));

        auto isStateValid = easy ? isStateValidEasy : isStateValidHard;

        si->setStateValidityChecker([isStateValid, si, this](const ob::State *state) {
            return isStateValid(si, state, *this);
        });

        si->setStateValidityCheckingResolution(stateValidityCheckingRes);

        //og::SimpleSetup ss(space);
        //si->setStateValidityChecker

        si->setup();

        ob::PlannerPtr optimizingPlanner = std::make_shared<og::RRTstar>(si);

        auto pdef(std::make_shared<ob::ProblemDefinition>(si));

        pdef->setOptimizationObjective(std::make_shared<ob::PathLengthOptimizationObjective>(si));

        pdef->setStartAndGoalStates(start, goal);
        //ss.setStartAndGoalStates(start, goal);

        optimizingPlanner->setProblemDefinition(pdef);

        optimizingPlanner->setup();

        //ss.setup();
        //ss.print();

        // attempt to solve the problem within solveTime seconds of planning time
        //ob::PlannerStatus solved = ss.solve(solveTime);

        pdef->print();

        ob::PlannerStatus solved = optimizingPlanner->solve(solveTime);

        if (solved)
        {
            // Output the length of the path found
            std::vector<double> reals;

            std::cout << "Found solution." << std::endl;
            //ss.simplifySolution();
            ob::PathPtr path = pdef->getSolutionPath();            
            //path.subdivide();
            //std::cout << path.getStateCount() << std::endl;
            //path.printAsMatrix(std::cout);

            // If a filename was specified, output the path as a matrix to
            // that file for visualization

            //ob::PlannerData plannerData(si);
            //ss.getPlannerData(plannerData);
            //plannerData.getVertex(1);
            //std::cout << ss.haveExactSolutionPath() << std::endl;

            //std::cout << ss.getPlanner() << std::endl;

            std::static_pointer_cast<og::PathGeometric>(pdef->getSolutionPath())->interpolate(interpolateNum);
            std::filebuf fb;
            fb.open("./../rspath.txt", std::ios::out);
            std::ostream os(&fb);
            std::static_pointer_cast<og::PathGeometric>(pdef->getSolutionPath())->printAsMatrix(os);
            fb.close();
        }
        else
        {
            std::cout << "No solution found." << std::endl;
        }
    }

} // namespace ReedsShepp
