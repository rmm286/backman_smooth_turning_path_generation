import numpy as np
from ODESolver import RungeKutta4
from StateSpaces import SmoothPathState


class SmoothPathPlanner:
    """ class for implementation of backman algorithm"""

    def __init__(self, kMax, kDotMax, kDotMin, vDotMax, vDotMin, initalState, finalState, kStart, kCenter, kEnd, reverse, headlandSpeed):
        self.kMax = kMax
        self.kDotMax = kDotMax
        self.kDotMin = kDotMin
        self.vDotMax = vDotMax
        self.vDotMin = vDotMin
        self.headlandSpeed = headlandSpeed
        self.path = []

        self.initialState = initalState
        self.finalState = finalState
        self.k_C0 = initalState.k
        self.k_C1 = kStart
        self.k_C2 = kCenter
        self.k_C3 = kEnd

    def generateCurvatureTrajectory(self, kInit, kFinal, dT, tInit):
        """
        Eqn. 4 from paper
        """
        kTolerance = 0.05

        k = kInit
        t = tInit

        kTrajectory = []

        while np.abs(k - kFinal) < kTolerance:
            if k < kFinal:
                k = k + dT*self.kDotMax
            else:
                k = k - dT*self.kDotMin
            t = t + dT

            kTrajectory.append([k, t])

        return kTrajectory

    def generateSpeedTrajectory(self, vInit, vFinal, dT, tInit):
        """
        Eqn. 7 from paper
        """
        vTolerance = 0.05

        v = vInit
        t = tInit

        vTrajectory = []

        while np.abs(v - vFinal) < vTolerance:
            if v < vFinal:
                v = v + dT*self.vDotMax
            else:
                v = v - dT*self.vDotMin
            t = t + dT

            vTrajectory.append([v, t])

        return vTrajectory

    def f(u,t):
        x = u[0]
        y = u[1]
        theta = u[2]

        

    def integrateTrajectory(kTrajectory,vTrajectory, Xo):
        
        X = []
        X.append(Xo)

        f = 
        
        if not (len(kTrajectory) == len(vTrajectory)):
            print("curvature and speed trajectories not same length")
            return 0
        else:
            
            return X

    def plan(self):

        path_is_not_feasible = True

        while path_is_not_feasible:
            K_S1 = self.generateCurvatureTrajectory(
                self.k_C0, self.k_C1, 0.005, 0)
            v_S1 = self.generateSpeedTrajectory(self.initialState.v, self.headlandSpeed, 0.005, 0)

            #make trajectories equal length
            v_S1 = v_S1[0:len(K_S1)]

            S1 = integrateTrajectory(K_S1,v_S1)

        return self.path


def main():

    # x pos., ypos., orientation, speed, curvature
    initialState = SmoothPathState(0, 0, 0.5*np.pi, 1, 0)
    finalState = SmoothPathState(2, 0, -0.5*np.pi, 1, 0)
    turningRadius = 1  # m
    kMax = 1/turningRadius
    kDotMax = 1  # max derivative of curvature
    kDotMin = -1  # min derivative of curvature
    vDotMax = 1
    vDotMin = -1
    headlandSpeed = 1

    RSRPathDetail = [KMax, 0, KMax, False]
    #LSLPathDetail = [KMax, 0, KMax, False]
    #LRLPathDetail = [KMax, KMax, KMax, False]
    #RLRPathDetail = [KMax, KMax, KMax, False]
    #LSRPathDetail = [KMax, 0, KMax, False]
    #RSLPathDetail = [KMax, 0, KMax, False]
    #R1L1RPAthDetail = [KMax, KMax, KMax, True]
    #L1R1LPAthDetail = [KMax, KMax, KMax, True]

    planSmoothInst = SmoothPathPlanner(kMax, initialState, finalState, kDotMax, kDotMin, vDotMax, vDotMin,
                                       RSRPathDetail[0],  RSRPathDetail[1],  RSRPathDetail[3],  RSRPathDetail[4], headlandSpeed)

    path = planSmoothInst.plan()

    print(path)


main()
