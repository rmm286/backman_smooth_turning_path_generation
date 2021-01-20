import numpy as np
from numpy import sin, cos, tan
import matplotlib.pyplot as plt
import scipy
from scipy import integrate
from scipy.integrate import odeint

class PathSegment: 

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

        #next line does a bunch of broadcasting to strip the ones from homogenous coords and add back the orientation data
        finalPath = rotAndTranslatedPoints * np.array([[1], [1], [0]]) + (self.poses.T * np.array([[0], [0], [1]])) + np.array([[0], [0], [thetaRotate]])

        self.poses = finalPath.T
    
    def rotateAboutPoint(self, a, b, theta):
        """ 
        Rotates the path about a point [a,b] by theta radians.

        Input: 
            
            a: x position of center of rotation (m)
            b: y position of center of rotation (m)
            theta: orientation of start or end pose (rad)
         """

        homogenousCoords = self.poses.T * np.array([[1], [1], [0]]) + np.array([[0], [0], [1]])

        H = np.array([[cos(theta), -1*sin(theta), -1*a*cos(theta) + b*sin(theta) + a],
                      [sin(theta), cos(theta), -1*a*sin(theta) - b*cos(theta) + b], 
                      [0, 0, 1]])

        rotatedPoints = np.matmul(H,homogenousCoords)
        finalPath = rotatedPoints* np.array([[1], [1], [0]]) + (self.poses.T * np.array([[0], [0], [1]])) + np.array([[0], [0], [theta]])

        self.poses = finalPath.T

class SpiralSegment(PathSegment):

    def __init__(self, kTrajectory, vTrajectory, xo):
        self.poses = self.integrateTrajectory(kTrajectory, vTrajectory, xo)
    
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

class CCSegment(PathSegment):

    def __init__(self, curvature, start, end):
        self.poses = np.array([[0],[0],[0]])