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

    def setStartAndGoal(self,initalState, finalState):
        if not (isinstance(initalState, SmoothPathState) and isinstance(initalState, SmoothPathState)):
            raise TypeError('Start and Goal not SE2States')
        else:
            self.initialState = initalState
            self.finalState = finalState
            self.k_C0 = initalState.k
            self.k_C4 = finalState.k

    def setNominalCurvatures(self,kStart, kCenter, kEnd, reverse):

        self.k_C1 = kStart
        self.k_C2 = kCenter
        self.k_C3 = kEnd
        self.reverse = reverse

    def generateCurvatureTrajectory(self, kInit, kFinal, tInit):
        """
        Eqn. 4 from paper

        TODO: optimize
        """
        kTolerance = 0.001

        k = kInit
        t = tInit

        kTrajectory = np.array([[k,t]])
        
        while np.abs(k - kFinal) > kTolerance:
            
            if k < kFinal:
                k = k + self.dT*self.kDotMax
            else:
                k = k + self.dT*self.kDotMin
            t = t + self.dT

            kTrajectory = np.append(kTrajectory, np.array([[k,t]]), axis = 0)

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

            vTrajectory = np.append(vTrajectory, np.array([[v,t]]), axis = 0)

        return vTrajectory

    def vehicleModel(self,x,t,u):
        v = u[0]
        k = u[1]
        theta = x[2]

        return [v*cos(theta), v*sin(theta), k]

    def integrateTrajectory(self, kTrajectory, vTrajectory, xo):

        if (not (len(kTrajectory) == len(vTrajectory))):
            raise ValueError("curvature and speed trajectories not same length")
            return 0
        else:

            t = [i[1] for i in kTrajectory]

            u = np.empty_like(kTrajectory)

            for i in range(len(kTrajectory)):
                
                u[i][0] = vTrajectory[i][0]
                u[i][1] = kTrajectory[i][0]

            x = np.empty([len(kTrajectory),3])
            x[0] = xo

            for i in range(1,len(t)):

                resultVal = odeint(self.vehicleModel, x[i-1], [t[i-1],t[i]], args=([u[i][0], u[i][1]],)) #integrate using LSODA
                
                x[i] = resultVal[1]

            return x

    def makeTrajectoriesEqualLength(self, kTraj, vTraj, fromStart = False):
        """
        TODO: increment time
        """

        if (len(kTraj) < len(vTraj)) and not fromStart: #cut from end of vTraj
            vTraj = vTraj[0:len(kTraj)]
        elif (len(kTraj) < len(vTraj)) and fromStart: #cut from start of vTraj
            vTraj = vTraj[len(vTraj)- len(kTraj):len(vTraj)]
        elif (len(kTraj) > len(vTraj)) and not fromStart:
            vTraj = np.append(vTraj, np.array([vTraj[-1] for i in range((len(kTraj)- len(vTraj)))]), axis = 0)
        elif (len(kTraj) > len(vTraj)) and fromStart:
            vTraj = np.append(np.array([vTraj[0] for i in range(len(kTraj)- len(vTraj))]),vTraj, axis = 0)

        return {'kTraj': kTraj, 'vTraj':vTraj}
    
    def makeTrajectoriesEqualLengthAndMatchZeros(self, kTraj, vTraj):
        """ 
        From section 2.2 curvature and velocity must be zero at point where velocity changes 
        From predefined turning types the curvature trajectory will always have a zero point when manuever involved reversing direction
        """

        #find zero point in kTraj
        kValArray = [np.abs(i[0]) for i in kTraj]
        vValArray = [np.abs(i[0]) for i in vTraj]
        kZeroIndex = np.argmin(kValArray)
        vZeroIndex = np.argmin(vValArray)

        kZeroTime = kTraj[kZeroIndex][1]

        if kZeroIndex < vZeroIndex:
            diff = vZeroIndex - kZeroIndex
            kTrajExtension = np.array([[0, kTraj[kZeroIndex] + (i+1)*dT] for i in range(diff)])
            kTraj = np.insert(kTraj,kZeroIndex, kTrajExtension, axis = 0)
            
            return self.makeTrajectoriesEqualLength(kTraj, vTraj)
        else:
            diff = kZeroIndex - vZeroIndex
            vTrajExtension = np.array([[0, vTraj[kZeroIndex] + (i+1)*dT] for i in range(diff)])
            vTraj = np.insert(vTraj,vZeroIndex, vTrajExtension, axis = 0)

            return self.makeTrajectoriesEqualLength(kTraj, vTraj)


    def relocatePath(self, path, point, relocateStartPoint = True):
        """Use Homogenous Transform to relocate a path to a point
            TODO: optimize function and use Homogenous TF efficiently
        """

        if not (len(path[1]) == len(point)):
            raise ValueError("Path and Point are not of suitable dimension")
        
        index = 0 if relocateStartPoint else -1
        pathPoint = path[index]
        tf = point - pathPoint
        H = np.array([ [cos(tf[2]), -sin(tf[2]), 0], [sin(tf[2]), cos(tf[2]), 0], [0, 0, 1] ])
        pathHomogenous = np.array([[i[0]-pathPoint[0], i[1]-pathPoint[1], 1] for i in path]).T 
        pathTF = np.matmul(H, pathHomogenous).T * [1, 1, 0]

        pathTF = (pathTF) + [point[0] - pathTF[-1][0], point[1] - pathTF[-1][1],  0]
        
        return pathTF


    def plan(self, dT):

        self.path_is_not_feasible = True
        self.dT = dT

        while self.path_is_not_feasible:
            #generate first connecting spiral
            K_S1 = self.generateCurvatureTrajectory(
                self.k_C0, self.k_C1, 0)
            v_S1 = self.generateSpeedTrajectory(
                self.initialState.v, self.headlandSpeed, 0)

            trajectories = self.makeTrajectoriesEqualLength(K_S1, v_S1, False)
            
            if len(v_S1) < len(trajectories['vTraj']):
                cut_v_S1 = True
                cut_v_S1_index = len(trajectories['vTraj'])
                
            xo = [self.initialState.x, self.initialState.y, self.initialState.theta]
            S1 = self.integrateTrajectory(trajectories['kTraj'], trajectories['vTraj'], xo)

            v_C1 = trajectories['vTraj'][-1]

            #generate last connecting spiral
            K_S4 = self.generateCurvatureTrajectory(self.k_C3, self.k_C4, 0)
            v_S4 = self.generateSpeedTrajectory(self.headlandSpeed, self.finalState.v, 0)
            trajectories = self.makeTrajectoriesEqualLength(K_S4, v_S4, True)
            xo = [0, 0, 0]
            S4 = self.integrateTrajectory(trajectories['kTraj'], trajectories['vTraj'], xo)
            relocationPoint = np.array([self.finalState.x, self.finalState.y, self.finalState.theta])
            S4 = self.relocatePath(S4, relocationPoint, False)

            #generate second conneting spiral
            K_S2 = self.generateCurvatureTrajectory(self.k_C1, self.k_C2, 0)
            xo = [0, 0, 0]

            if self.reverse:
                v_S2 = generateSpeedTrajectory(self.headlandSpeed, 0, 0)
                v_S2 = np.append(v_S2, generateSpeedTrajectory(0, self.headlandSpeedReverse, v_S2[-1][1]), axis = 0)
                trajectories = self.makeTrajectoriesEqualLengthAndMatchZeros(K_S2, v_S2)
                S2 = self.integrateTrajectory(trajectories['kTraj'], trajectories['vTraj'], xo)
            else:
                if cut_v_S1:
                    v_S2 = v_S1[cut_v_S1_index:-1]
                    #v_S2 = np.append(v_S2,[[v_S1[-1][0], 0] for i in range(len(K_S2) - len(v_S2))])
                    print(len(v_S2))
                    print(len(K_S2))
                    
                else:
                    v_S2 = [[v_C1, 0] for i in len(K_S2)]
                
                S2 = self.integrateTrajectory(K_S2, v_S2, xo)

            self.path_is_not_feasible = False

            
        ##plotting stuff
        plt.plot([i[0] for i in S1], [i[1] for i in S1])
        plt.xlim([-0.1,1])
        plt.ylim([-0.1,1])
        plt.savefig("trajectory.png")

        print("theta final:", S1[-1][2])
        print("k final:", (S1[-2][2]-S1[-1][2])/self.dT)

        plt.plot([i[0] for i in S4], [i[1] for i in S4])
        #plt.plot([i*0.01 for i in range(0,200)], [0 for i in range(0,200)])
        plt.xlim([-2.5,2.5])
        plt.ylim([-2.5,2.5])
        plt.savefig("trajectory3.png")

        print("theta final:", S4[-1][2])
        print("k final:", (S4[-2][2]-S4[-1][2])/self.dT)

        return self.path

def main():

    # x pos., ypos., orientation, speed, curvature
    initialState = SmoothPathState(0, 0, 0.5*np.pi, 1, 0)
    finalState = SmoothPathState(2, 0, -0.5*np.pi, 1, 0)
    turningRadius = 1  # m
    dT = 0.05
    
    kMax = 1/turningRadius
    kMin = -kMax
    kDotMax = 1  # max derivative of curvature
    kDotMin = -1  # min derivative of curvature
    kDDotMax = 1
    kDDotMin = -1

    kConstraints = [kMax, kMin, kDotMax, kDotMin, kDDotMax, kDDotMin]
    
    vMax = 1
    vMin = -1
    vDotMax = 1
    vDotMin = -1
    vDDotMax = 1
    vDDotMin = -1
    headlandSpeed = vMax
    headlandSpeedReverse = vMin

    vConstraints = [vMax, vMin, vDotMax, vDotMin, vDDotMax, vDDotMin]

    RSRPathDetail = [kMin, 0, kMin, False]
    #LSLPathDetail = [KMax, 0, KMax, False]
    #LRLPathDetail = [KMax, KMin, KMax, False]
    #RLRPathDetail = [KMin, KMax, KMin, False]
    #LSRPathDetail = [KMax, 0, KMin, False]
    #RSLPathDetail = [KMin, 0, KMax, False]
    #R1L1RPAthDetail = [KMin, KMax, KMin, True]
    #L1R1LPAthDetail = [KMax, KMin, KMax, True]

    planSmoothInst = SmoothPathPlanner()

    planSmoothInst.setConstraints(kConstraints,vConstraints,headlandSpeed,headlandSpeedReverse  )
    planSmoothInst.setNominalCurvatures(RSRPathDetail[0], RSRPathDetail[1], RSRPathDetail[2], RSRPathDetail[3])
    planSmoothInst.setStartAndGoal(initialState, finalState)

    path = planSmoothInst.plan(dT)

    #print(path)


main()
