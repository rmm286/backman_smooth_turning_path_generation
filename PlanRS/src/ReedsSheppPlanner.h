#include <ompl/base/spaces/ReedsSheppStateSpace.h>
#include <ompl/base/ScopedState.h>
#include <ompl/geometric/SimpleSetup.h>

namespace ReedsShepp
{
    class RSPlanner
    { //smooth path planner type
    public:
        RSPlanner(float arg_turningRadius, float arg_startX, float arg_startY, float arg_startTh, float arg_goalX, float arg_goalY, float arg_goalTh, float arg_xHigh, float arg_yHigh, float arg_robotWidth, float arg_robotLength);

        float m_turningRadius;
        float m_startX;
        float m_startY;
        float m_startTh;
        float m_goalX;
        float m_goalY;
        float m_goalTh;
        float m_xHigh;
        float m_yHigh;
        float m_robotWidth;
        float m_robotLength;

        void planRS(int easy, float solveTime, float stateValidityCheckingRes, int interpolateNum);
    };

    bool isStateValidEasy(ompl::base::SpaceInformationPtr si, const ompl::base::State *state, RSPlanner thisRSPlanner);

    bool isStateValidHard(ompl::base::SpaceInformationPtr si, const ompl::base::State *state, RSPlanner thisRSPlanner);
} // namespace ReedsShepp