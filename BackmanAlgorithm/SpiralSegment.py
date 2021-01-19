import numpy as np
from numpy import sin, cos, tan
import matplotlib.pyplot as plt
import scipy
from scipy import integrate
from scipy.integrate import odeint

class SpiralSegment:

    def __init__(self, kTrajectory, vTrajectory, xo):
        self.poses = self.integrateTrajectory(kTrajectory, vTrajectory, xo)

    def placePath(self, x, y, theta, wrtStart=True):
        """ 
        Places the path at a point 

        Input: 
            theta: orientation of start or end pose(rad)
            x: x position of start or end pose (m)
            y: y position of start or end pose (m)
            wrtStart: perform operation wrt to the start point (bool)

         """

        index = 0 if wrtStart else -1

        thetaRotate = theta - self.poses[index][2]
        a = self.poses[index][0]
        b = self.poses[index][1]

        h =  x - a
        k = y - b

        homogenousCoords = self.poses.T * np.array([[1], [1], [0]]) + np.array([[0], [0], [1]])

        H = np.array([[cos(thetaRotate), -1*sin(thetaRotate), -1*a*cos(thetaRotate) + b*sin(thetaRotate) + a],
                      [sin(thetaRotate), cos(thetaRotate), -1*a*sin(thetaRotate) - b*cos(thetaRotate) + b], 
                      [0, 0, 1]])
        T = np.array([[1, 0, h],[ 0, 1, k], [0,0,1]])

        TxH = np.matmul(T,H)

        rotAndTranslatedPoints = np.matmul(TxH,homogenousCoords)
        finalPath = rotAndTranslatedPoints + (self.poses.T * np.array([[0], [0], [1]])) + np.array([[0], [0], [thetaRotate]])

        self.poses = finalPath.T

    def integrateTrajectory(self, kTrajectory, vTrajectory, xo):
        """ 
        Performs integration on dyanmic model using LSODA from odeint. kTrajetory and vTrajectory must be of equal length.

        Input:
            kTrajectory: Array of [curvature, time] values for each timestep
            vTrajectory: Array of [velocity, time] values for each timestep
            xo: Inital pose of model, x[0] = x, x[1] = y, x[2] = theta (orientation)

        Output:
            x: A Nx3 numpy array of [x,y,theta] poses for each timestep.
        """

        if (not (len(kTrajectory) == len(vTrajectory))):
            raise ValueError(
                "curvature and speed trajectories not same length")

        t = [i[1] for i in kTrajectory]
        u = np.empty_like(kTrajectory)
        for i in range(len(kTrajectory)):
            u[i][0] = vTrajectory[i][0]
            u[i][1] = kTrajectory[i][0]
        x = np.empty([len(kTrajectory), 3])
        x[0] = xo
        for i in range(1, len(t)):
            resultVal = odeint(
                self.vehicleModel, x[i-1], [t[i-1], t[i]], args=([u[i][0], u[i][1]],))
            x[i] = resultVal[1]

        return x
    
    def vehicleModel(self, x, t, u):
        """ 
        Dynamic model for advancing integrator. Based on common bicycle model.

        Input:
            x: current state of model, x[0] = x, x[1] = y, x[2] = theta (orientation)
            t: time point, not used in model
            u: Control vector, u[0] = velocity control value u[1] = curvature control value

        Output:
            x: next state of model based on control and previous state, x[0] = x, x[1] = y, x[2] = theta (orientation)
        """

        v = u[0]
        k = u[1]
        theta = x[2]

        return [v*cos(theta), v*sin(theta), k]

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