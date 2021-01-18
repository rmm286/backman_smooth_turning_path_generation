import numpy as np
from numpy import sin, cos, tan
import matplotlib.pyplot as plt

class PathSegment:

    def __init__(self, kTrajectory, vTrajectory, xo):
        self.poses = integrateTrajectory(kTrajectory, vTrajectory, xo)


    def placePath(self, x, y, theta, wrtStart = True):
        """ 
        Places the path at a point 

        Input: 
            theta: orientation of start or end pose(rad)
            x: x position of start or end pose (m)
            y: y position of start or end pose (m)
            wrtStart: perform operation wrt to the start point (bool)

         """
         self.rotatePath()
         self.translatePath()

         return 0



    def rotatePath(self, rot, a, b, wrtStart = True):
        """ 
        Rotates the path about point (a,b) 

        Input: 
            rot: rotation about point [a,b] (rad)
            a: x position of center of rotation (m)
            b: y position of center of rotation (m)
            wrtStart: perform operation wrt to the start point (bool)
        """
        return 0



    def translatePath(self, xTrans, yTrans, wrtStart = True):
        """
        Translates the path to point 
        Input: 
            xTrans: displacement in x (m)
            yTrans: displacement in y (m)
            wrtStart: perform operation wrt to the start point (bool)
        """
        return 0

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