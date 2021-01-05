import numpy as np
from ODESolver import RungeKutta4
from StateSpaces import SmoothPathState


class SmoothPathPlanner:
    """ class for implementation of backman algorithm"""

    def __init__(self):

        self.path = []

    def setConstraints(kMax, kMin, kDotMax, kDotMin, vDotMax, vDotMin):
        self.kMax = kMax
        self.kMin = kMin
        self.kDotMax = kDotMax
        self.kDotMin = kDotMin
        self.vDotMax = vDotMax
        self.vDotMin = vDotMin

    def setStartAndGoal(initalState, finalState):
        if not (isinstance(initalState, StateSpaces.SE2State) and isinstance(initalState, StateSpaces.SE2State)):
            raise TypeError('Start and Goal not SE2States')
        else:
            self.initialState = initalState
            self.finalState = finalState
            self.headlandSpeed = finalState.v
            self.k_C0 = initalState.k

    def setNominalCurvatures(kStart, kCenter, kEnd, reverse):

        self.k_C1 = kStart
        self.k_C2 = kCenter
        self.k_C3 = kEnd
        self.reverse = reverse

    def generateCurvatureTrajectory(self, kInit, kFinal, dT, tInit):
        """
        Eqn. 4 from paper

        TODO: remove reference to tInit
        """
        kTolerance = 0.05

        k = kInit
        t = tInit

        kTrajectory = []

        while np.abs(k - kFinal) < kTolerance:
            if k < kFinal:
                k = k + dT*self.kDotMax
            else:
                k = k + dT*self.kDotMin
            t = t + dT

            kTrajectory.append([k, t])

        return kTrajectory

    def generateSpeedTrajectory(self, vInit, vFinal, dT, tInit):
        """
        Eqn. 7 from paper

        TODO: remove reference to tInit
        """
        vTolerance = 0.05

        v = vInit
        t = tInit

        vTrajectory = []

        while np.abs(v - vFinal) < vTolerance:
            if v < vFinal:
                v = v + dT*self.vDotMax
            else:
                v = v + dT*self.vDotMin
            t = t + dT

            vTrajectory.append([v, t])

        return vTrajectory

    def f(u, t):
        v, k = u[0], u[1]

        thetaDot =
        y =
        theta = u[2]

    def integrateTrajectory(kTrajectory, vTrajectory, Xo):

        X = []
        X.append(Xo)

        if (not (len(kTrajectory) == len(vTrajectory))):
            print("curvature and speed trajectories not same length")
            return 0
        else:

            return X

    def plan(self):

        path_is_not_feasible = True

        while path_is_not_feasible:
            K_S1 = self.generateCurvatureTrajectory(
                self.k_C0, self.k_C1, 0.005, 0)
            v_S1 = self.generateSpeedTrajectory(
                self.initialState.v, self.headlandSpeed, 0.005, 0)

            # make trajectories equal length
            v_S1 = v_S1[0:len(K_S1)]

            #S1 = integrateTrajectory(K_S1,v_S1)

        return self.path


def main():

    # x pos., ypos., orientation, speed, curvature
    initialState = SmoothPathState(0, 0, 0.5*np.pi, 1, 0)
    finalState = SmoothPathState(2, 0, -0.5*np.pi, 1, 0)
    turningRadius = 1  # m
    kMax = 1/turningRadius
    kMin = -kMax
    kDotMax = 1  # max derivative of curvature
    kDotMin = -1  # min derivative of curvature
    vDotMax = 1
    vDotMin = -1
    headlandSpeed = 1

    RSRPathDetail = [kMin, 0, kMin, False]
    #LSLPathDetail = [KMax, 0, KMax, False]
    #LRLPathDetail = [KMax, KMin, KMax, False]
    #RLRPathDetail = [KMin, KMax, KMin, False]
    #LSRPathDetail = [KMax, 0, KMin, False]
    #RSLPathDetail = [KMin, 0, KMax, False]
    #R1L1RPAthDetail = [KMin, KMax, KMin, True]
    #L1R1LPAthDetail = [KMax, KMin, KMax, True]

    planSmoothInst = SmoothPathPlanner(kMax, kMin, initialState, finalState, kDotMax, kDotMin, vDotMax, vDotMin,
                                       RSRPathDetail[0],  RSRPathDetail[1],  RSRPathDetail[3],  RSRPathDetail[4], headlandSpeed)

    path = planSmoothInst.plan()

    print(path)


main()
