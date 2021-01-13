import numpy as np
from numpy import sin, cos, tan
from ODESolver import RungeKutta4
from StateSpaces import SmoothPathState
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


class SmoothPathPlanner:
    """ class for implementation of backman algorithm"""

    def __init__(self):

        self.path = []

    def setConstraints(self, kConstraints, vConstraints, headlandSpeed, headlandSpeedReverse):

        self.kMax = kConstraints[0]
        self.kMin = kConstraints[1]
        self.kDotMax = kConstraints[2]
        self.kDotMin = kConstraints[3]
        self.kDDotMax = kConstraints[4]
        self.kDDotMin = kConstraints[5]

        self.vMax = kConstraints[0]
        self.vMin = kConstraints[1]
        self.vDotMax = kConstraints[2]
        self.vDotMin = kConstraints[3]
        self.vDDotMax = kConstraints[4]
        self.vDDotMin = kConstraints[5]

        self.headlandSpeed = headlandSpeed
        self.headlandSpeedReverse = headlandSpeedReverse

    def setStartAndGoal(self, initalState, finalState):
        if not (isinstance(initalState, SmoothPathState) and isinstance(initalState, SmoothPathState)):
            raise TypeError('Start and Goal not SE2States')
        else:
            self.initialState = initalState
            self.finalState = finalState
            self.k_C0 = initalState.k
            self.k_C4 = finalState.k

    def setNominalCurvatures(self, kStart, kCenter, kEnd, reverse):

        self.k_C1 = kStart
        self.k_C2 = kCenter
        self.k_C3 = kEnd
        self.reverse = reverse

    def generateCurvatureTrajectory(self, kInit, kFinal, tInit):
        """
        Eqn. 4 from paper

        TODO: optimize
        """
        kTolerance = 0.05

        k = kInit
        t = tInit

        kTrajectory = np.array([[k, t]])

        while np.abs(k - kFinal) > kTolerance:

            if k < kFinal:
                k = k + self.dT*self.kDotMax
            else:
                k = k + self.dT*self.kDotMin
            t = t + self.dT

            kTrajectory = np.append(kTrajectory, np.array([[k, t]]), axis=0)

        return kTrajectory

    def generateSpeedTrajectory(self, vInit, vFinal, tInit):
        """
        Eqn. 7 from paper

        TODO: optimize
        """
        vTolerance = 0.05

        v = vInit
        t = tInit

        vTrajectory = np.array([[v, t]])

        while np.abs(v - vFinal) > vTolerance:
            if v < vFinal:
                v = v + self.dT*self.vDotMax
            else:
                v = v + self.dT*self.vDotMin
            t = t + self.dT

            vTrajectory = np.append(vTrajectory, np.array([[v, t]]), axis=0)

        return vTrajectory

    def vehicleModel(self, x, t, u):
        v = u[0]
        k = u[1]
        theta = x[2]

        return [v*cos(theta), v*sin(theta), k]

    def integrateTrajectory(self, kTrajectory, vTrajectory, xo):

        if (not (len(kTrajectory) == len(vTrajectory))):
            raise ValueError(
                "curvature and speed trajectories not same length")
            return 0
        else:

            t = [i[1] for i in kTrajectory]

            u = np.empty_like(kTrajectory)

            for i in range(len(kTrajectory)):

                u[i][0] = vTrajectory[i][0]
                u[i][1] = kTrajectory[i][0]

            x = np.empty([len(kTrajectory), 3])
            x[0] = xo

            for i in range(1, len(t)):

                # integrate using LSODA
                resultVal = odeint(
                    self.vehicleModel, x[i-1], [t[i-1], t[i]], args=([u[i][0], u[i][1]],))

                x[i] = resultVal[1]

            return x

    def makeTrajectoriesEqualLength(self, kTraj, vTraj, fromStart=False):
        """
        TODO: increment time
        """

        if (len(kTraj) < len(vTraj)) and not fromStart:  # cut from end of vTraj
            vTraj = vTraj[0:len(kTraj)]
        elif (len(kTraj) < len(vTraj)) and fromStart:  # cut from start of vTraj
            vTraj = vTraj[len(vTraj) - len(kTraj):len(vTraj)]
        elif (len(kTraj) > len(vTraj)) and not fromStart:
            vTraj = np.append(vTraj, np.array(
                [vTraj[-1] for i in range((len(kTraj) - len(vTraj)))]), axis=0)
        elif (len(kTraj) > len(vTraj)) and fromStart:
            vTraj = np.append(
                np.array([vTraj[0] for i in range(len(kTraj) - len(vTraj))]), vTraj, axis=0)

        return {'kTraj': kTraj, 'vTraj': vTraj}

    def makeTrajectoriesEqualLengthAndMatchZeros(self, kTraj, vTraj):
        """
        From section 2.2 curvature and velocity must be zero at point where velocity changes
        From predefined turning types the curvature trajectory will always have a zero point when manuever involved reversing direction
        """

        # find zero point in kTraj
        kValArray = [np.abs(i[0]) for i in kTraj]
        vValArray = [np.abs(i[0]) for i in vTraj]
        kZeroIndex = np.argmin(kValArray)
        vZeroIndex = np.argmin(vValArray)

        kZeroTime = kTraj[kZeroIndex][1]

        if kZeroIndex < vZeroIndex:
            diff = vZeroIndex - kZeroIndex
            kTrajExtension = np.array(
                [[0, kTraj[kZeroIndex] + (i+1)*dT] for i in range(diff)])
            kTraj = np.insert(kTraj, kZeroIndex, kTrajExtension, axis=0)

            return self.makeTrajectoriesEqualLength(kTraj, vTraj)
        else:
            diff = kZeroIndex - vZeroIndex
            vTrajExtension = np.array(
                [[0, vTraj[kZeroIndex] + (i+1)*dT] for i in range(diff)])
            vTraj = np.insert(vTraj, vZeroIndex, vTrajExtension, axis=0)

            return self.makeTrajectoriesEqualLength(kTraj, vTraj)

    def relocatePath(self, path, point, relocateStartPoint=True):
        """Use Homogenous Transform to relocate a path to a point
            TODO: optimize function and use Homogenous TF efficiently
            TODO: get orienttion values
        """

        if not (len(path[1]) == len(point)):
            raise ValueError("Path and Point are not of suitable dimension")

        index = 0 if relocateStartPoint else -1
        pathPoint = path[index]
        tf = point - pathPoint
        H = np.array([[cos(tf[2]), -sin(tf[2]), 0],
                      [sin(tf[2]), cos(tf[2]), 0], [0, 0, 1]])
        pathHomogenous = np.array(
            [[i[0]-pathPoint[0], i[1]-pathPoint[1], 1] for i in path]).T
        pathTF = np.matmul(H, pathHomogenous).T * [1, 1, 0] + np.array([[0, 0, i[2]] for i in path])

        pathTF = (pathTF) + [point[0] - pathTF[-1]
                             [0], point[1] - pathTF[-1][1],  tf[2]]

        return pathTF

    def calculateCenterArc(self):
        S1 = self.S1
        S2 = self.S2
        S3 = self.S3
        S4 = self.S4
        k_C1 = self.k_C1
        k_C2 = self.k_C2
        k_C3 = self.k_C3

        omega_S2_tC2 = np.array(
            [S2[-1][0] - (k_C2**-1)*sin(S2[-1][2]), S2[-1][1] + (k_C2**-1)*cos(S2[-1][2])])
        omega_S2_tS2 = np.array(
            [S2[0][0] - (k_C1**-1)*sin(S2[0][2]), S2[0][1] + (k_C1**-1)*cos(S2[0][2])])

        omega_S3_tC3 = np.array(
            [S3[-1][0] - (k_C3**-1)*sin(S3[-1][2]), S3[-1][1] + (k_C3**-1)*cos(S3[-1][2])])
        omega_S3_tS3 = np.array(
            [S3[0][0] - (k_C2**-1)*sin(S3[0][2]), S3[0][1] + (k_C2**-1)*cos(S3[0][2])])

        self.omega_k = np.array(
            [S1[-1][0] - (k_C1**-1)*sin(S1[-1][2]), S1[-1][1] + (k_C1**-1)*cos(S1[-1][2])])
        self.omega_kplus2 = np.array(
            [S4[0][0] - (k_C3**-1)*sin(S4[0][2]), S4[0][1] + (k_C3**-1)*cos(S4[0][2])])

        d1 = np.linalg.norm(self.omega_kplus2 - self.omega_k)
        d2 = np.linalg.norm(omega_S2_tC2 - omega_S2_tS2)
        d3 = np.linalg.norm(omega_S3_tC3 - omega_S3_tS3)

        l1 = (d2**2 - d3**2 + d1**2)/(2*d1)
        if(l1 > d2):
            return False
        l2 = np.sqrt(d2**2 - l1**2)

        signVal = 1 if (k_C2 < 0) else -1
        self.omega_kplus1 = self.omega_k + l1*(self.omega_kplus2 - self.omega_k)/d1 + signVal*l2*np.matmul(
            np.array([[0, -1], [1, 0]]), (self.omega_kplus2 - self.omega_k))/d1

        return self.omega_kplus1

    def searchC1C2(self):
        S1 = self.S1
        S2 = self.S2
        S3 = self.S3
        S4 = self.S4
        k_C1 = self.k_C1
        k_C2 = self.k_C2
        k_C3 = self.k_C3
        omega_1 = self.omega_k
        omega_2 = self.omega_kplus1
        omega_3 = self.omega_kplus2

        searchResolution = 200
        serachTolerance = 0.5
        thetaSearch = np.linspace(0, 2*np.pi, searchResolution)

        thetaAB = S2[-1][2] - S2[0][2]
        rAB = np.array([S2[-1][0] - S2[0][0], S2[-1][1] - S2[0][1]])
        
        bestNorm = 100000

        for theta in thetaSearch:
            thetaT = theta + np.pi/2.0
            R = np.array([[cos(thetaT), -sin(thetaT)], [sin(thetaT), cos(thetaT)]])
            rAB = np.matmul(R,rAB)
            pointA = np.array(
                [(k_C1**-1)*cos(theta) + omega_1[0], (k_C1**-1)*sin(theta) + omega_1[1]])
            thetaTangentA = np.pi/2.0 + theta
            pointB = np.array([pointA[0] + rAB[0], pointA[0] + rAB[1]])
            thetaTangentB = np.pi/2.0 + thetaAB

            theta2 = -1*np.pi/2.0 - thetaTangentB

            pointOnCircle2 = np.array(
                [-1*(k_C2**-1)*cos(-1*theta2) + omega_2[0], (k_C2**-1)*sin(-1*theta2) + omega_2[1]])

            currentNorm = np.linalg.norm(pointB - pointOnCircle2)
            if (currentNorm < bestNorm):
                bestNorm = currentNorm
                bestPoint = pointA

        if (bestNorm < serachTolerance):
            return np.array([bestPoint[0], bestPoint[1], -1*thetaTangentB])
        else:
            raise ValueError("Search of C1 C2 failed.")
            print("Best Norm", bestNorm)

    def plan(self, dT):

        self.path_is_not_feasible = True
        self.dT = dT

        while self.path_is_not_feasible:

            ################ generate first connecting spiral ################
            K_S1 = self.generateCurvatureTrajectory(
                self.k_C0, self.k_C1, 0)
            v_S1 = self.generateSpeedTrajectory(
                self.initialState.v, self.headlandSpeed, 0)

            trajectories = self.makeTrajectoriesEqualLength(K_S1, v_S1, False)
            if len(v_S1) > len(trajectories['vTraj']):
                cut_v_S1 = True
                cut_v_S1_index = len(trajectories['vTraj'])
            else:
                cut_v_S1 = False
            xo = [self.initialState.x, self.initialState.y,
                  self.initialState.theta]
            self.S1 = self.integrateTrajectory(
                trajectories['kTraj'], trajectories['vTraj'], xo)
            v_C1 = trajectories['vTraj'][-1][0]

            ################ generate last connecting spiral ################
            K_S4 = self.generateCurvatureTrajectory(self.k_C3, self.k_C4, 0)
            v_S4 = self.generateSpeedTrajectory(
                self.headlandSpeed, self.finalState.v, 0)
            trajectories = self.makeTrajectoriesEqualLength(K_S4, v_S4, True)

            if len(v_S4) > len(trajectories['vTraj']):
                cut_v_S4 = True
                cut_v_S4_index = len(v_S4) - len(trajectories['vTraj'])
            else:
                cut_v_S4 = False

            v_C3 = trajectories['vTraj'][0][0]  # first value of v_s4

            xo = [0, 0, 0]
            self.S4 = self.integrateTrajectory(
                trajectories['kTraj'], trajectories['vTraj'], xo)
            relocationPoint = np.array(
                [self.finalState.x, self.finalState.y, self.finalState.theta])
            self.S4 = self.relocatePath(self.S4, relocationPoint, False)

            ################ generate second conneting spiral ################
            K_S2 = self.generateCurvatureTrajectory(self.k_C1, self.k_C2, 0)
            xo = [0, 0, 0]

            if self.reverse:
                v_S2 = generateSpeedTrajectory(self.headlandSpeed, 0, 0)
                v_S2 = np.append(v_S2, generateSpeedTrajectory(
                    0, self.headlandSpeedReverse, v_S2[-1][1]), axis=0)
                trajectories = self.makeTrajectoriesEqualLengthAndMatchZeros(
                    K_S2, v_S2)
                self.S2 = self.integrateTrajectory(
                    trajectories['kTraj'], trajectories['vTraj'], xo)
            else:
                if cut_v_S1:
                    v_S2 = v_S1[cut_v_S1_index:-1]
                    len_v_S2 = len(v_S2)
                    extend_v_S2 = np.array([[v_S1[-1][0], 0]
                                            for i in range(len(K_S2) - len_v_S2)])
                    v_S2 = np.append(v_S2, extend_v_S2, axis=0)
                    # print(len(v_S2))
                    # print(len(K_S2))
                else:
                    v_S2 = np.array([[v_C1, 0] for i in range(len(K_S2))])

                self.S2 = self.integrateTrajectory(K_S2, v_S2, xo)

            ################ generate thrid connecting spiral ################
            K_S3 = self.generateCurvatureTrajectory(self.k_C2, self.k_C3, 0)

            if self.reverse:
                v_S3 = generateSpeedTrajectory(self.headlandSpeedReverse, 0, 0)
                v_S3 = np.append(v_S3, generateSpeedTrajectory(
                    0, self.headlandSpeed, v_S3[-1][1]), axis=0)
                trajectories = self.makeTrajectoriesEqualLengthAndMatchZeros(
                    K_S3, v_S3)
                self.S3 = self.integrateTrajectory(
                    trajectories['kTraj'], trajectories['vTraj'], xo)
            else:
                if cut_v_S4:
                    v_S3 = v_S4[0:cut_v_S4_index]
                    len_v_S3 = len(v_S3)
                    extend_v_S3 = np.array([[v_S4[-1][0], 0]
                                            for i in range(len(K_S3) - len_v_S3)])
                    v_S3 = np.append(v_S3, extend_v_S3, axis=0)
                    # print(len(v_S2))
                    # print(len(K_S2))
                else:
                    v_S3 = np.array([[v_C3, 0] for i in range(len(K_S3))])

                self.S3 = self.integrateTrajectory(K_S3, v_S3, xo)

            ################ generate center arc ################
            if np.abs(self.k_C2) > 0.05:  # center is an arc
                omega_C2 = self.calculateCenterArc()

                if omega_C2.all() == False:
                    return False

            relocatePoint = self.searchC1C2()
            self.S2 = self.relocatePath(self.S2, relocatePoint)

            self.path_is_not_feasible = False

        # plotting stuff
        plt.plot([i[0] for i in self.S1], [i[1] for i in self.S1])
        plt.plot([i[0] for i in self.S2], [i[1] for i in self.S2])
        #plt.plot([i[0] for i in self.S3], [i[1] for i in self.S3])
        plt.plot([i[0] for i in self.S4], [i[1] for i in self.S4])

        plt.plot(self.omega_k[0], self.omega_k[1], 'b^')
        plt.plot(self.omega_kplus1[0], self.omega_kplus1[1], 'b^')
        plt.plot(self.omega_kplus2[0], self.omega_kplus2[1], 'b^')

        plt.plot([(self.k_C1**-1)*cos(theta) + self.omega_k[0] for theta in np.linspace(0, 2*np.pi, 25)],
                 [(self.k_C1**-1)*sin(theta) + self.omega_k[1] for theta in np.linspace(0, 2*np.pi, 25)], 'r--')
        plt.plot([(self.k_C2**-1)*cos(theta) + self.omega_kplus1[0] for theta in np.linspace(0, 2*np.pi, 25)],
                 [(self.k_C2**-1)*sin(theta) + self.omega_kplus1[1] for theta in np.linspace(0, 2*np.pi, 25)], 'r--')
        plt.plot([(self.k_C3**-1)*cos(theta) + self.omega_kplus2[0] for theta in np.linspace(0, 2*np.pi, 25)],
                 [(self.k_C3**-1)*sin(theta) + self.omega_kplus2[1] for theta in np.linspace(0, 2*np.pi, 25)], 'r--')

        plt.xlim([-1, 4])
        plt.ylim([-1, 4])
        plt.savefig("trajectory.png")

        return self.path


def main():

    # x pos., ypos., orientation, speed, curvature
    initialState = SmoothPathState(0, 0, 0.5*np.pi, 1, 0)
    finalState = SmoothPathState(1.5, 0, -0.5*np.pi, 1, 0)
    turningRadius = 1.0  # m
    dT = 0.005

    kMax = 1/float(turningRadius)
    kMin = -kMax
    kDotMax = 5.0  # max derivative of curvature
    kDotMin = -kDotMax  # min derivative of curvature
    kDDotMax = 1.0
    kDDotMin = -kDDotMax

    kConstraints = [kMax, kMin, kDotMax, kDotMin, kDDotMax, kDDotMin]

    vMax = 1.0
    vMin = -vMax
    vDotMax = 1.0
    vDotMin = -.75
    vDDotMax = 1.0
    vDDotMin = -1.0
    headlandSpeed = vMax + 0.5
    headlandSpeedReverse = vMin - 0.25

    vConstraints = [vMax, vMin, vDotMax, vDotMin, vDDotMax, vDDotMin]

    RSR = [kMin, 0.0, kMin, False]
    # LSL = [KMax, 0, KMax, False]
    LRL = [kMax, kMin, kMax, False]
    # RLR = [KMin, KMax, KMin, False]
    # LSR = [KMax, 0, KMin, False]
    # RSL = [KMin, 0, KMax, False]
    # R1L1R = [KMin, KMax, KMin, True]
    # L1R1L = [KMax, KMin, KMax, True]

    planSmoothInst = SmoothPathPlanner()

    pathType = LRL
    planSmoothInst.setConstraints(
        kConstraints, vConstraints, headlandSpeed, headlandSpeedReverse)
    planSmoothInst.setNominalCurvatures(
        pathType[0], pathType[1], pathType[2], pathType[3])
    planSmoothInst.setStartAndGoal(initialState, finalState)

    path = planSmoothInst.plan(dT)

    # print(path)


main()
