import numpy as np
from numpy import sin, cos, tan
import matplotlib.pyplot as plt
from scipy.integrate import odeint

class PathSegment: 

    def placePath(self, x, y, theta, wrtStart=True):
        """ 
        Places a path so that either the first or last element of the path is equal to the given [x,y,th].

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

        #next line does a bunch of broadcasting to strip the ones from homogenous coords and add back the orientation values
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
    
    def pathIntersectsWith(self, otherPath):
        """
        Function returns true if the two paths intersect. This is implemented as a simple bounding box, since paths shouldn't really get very close to each other.

        Output:
            intersects: True if paths intersect (bool)
        
        """
        leftSideBox1 = np.amin(np.array([i[0] for i in self.poses]))
        rightSideBox1 = np.amax(np.array([i[0] for i in self.poses]))
        leftSideBox2 = np.amin(np.array([i[0] for i in otherPath.poses]))
        rightSideBox2 = np.amax(np.array([i[0] for i in otherPath.poses]))
        bottomSideBox1 = np.amin(np.array([i[1] for i in self.poses]))
        topSideBox1 = np.amax(np.array([i[1] for i in self.poses]))
        bottomSideBox2 = np.amin(np.array([i[1] for i in otherPath.poses]))
        topSideBox2 = np.amax(np.array([i[1] for i in otherPath.poses]))

        return not (rightSideBox1 < leftSideBox2 or leftSideBox1 > rightSideBox2 or topSideBox1 < bottomSideBox2 or bottomSideBox1 > topSideBox2)

            

class SpiralSegment(PathSegment):

    def __init__(self, kTrajectory, vTrajectory, xo):
        self.poses = self.integrateTrajectory(kTrajectory, vTrajectory, xo)
        self.controls = np.array([[kTrajectory[i][0], vTrajectory[i][0]] for i in range(len(kTrajectory))]) #TODO: optimize with array operations
    
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

        return [v*cos(theta), v*sin(theta), k*v]

class CCSegment(PathSegment):

    def __init__(self, curvature, vTraj, start, end, center, dT):
        """
        Generates a constant curvature segment from start to end coordinates, with given curvature, center of rotation and speed profile. 

        Input: 
            curvature: curvature of constant arc (float)
            vTraj: velocity profile of constant arc segment (Nx2 array of floats)
            start: SE2 position of start point (1x2 array of floats)
            end: SE2 position of end point (1x2 array of floats)
            center: R2 position of center of turning (1x2 of floats)
            dT: timestep, must be timestep of overall planning profile (float)
        """
        
        ang1 = np.arctan2(start[1]-center[1],start[0]-center[0])
        ang2 = np.arctan2(end[1]-center[1], end[0]-center[0])
        r = np.linalg.norm(start[0:2]-center)

        if(curvature > 0 and ang2 < ang1):
            ang2 = ang2 + 2*np.pi
        
        if(curvature < 0 and ang2 > ang1):
            ang2 = ang2 - 2*np.pi

        arcLen = r*np.abs(ang2-ang1)
        maxTimeToTraverse = arcLen/np.amin(vTraj.T[0])
        maxStepsToTraverse = int(maxTimeToTraverse/dT)

        self.poses = np.zeros([maxStepsToTraverse+1,3])
        self.controls = np.zeros([maxStepsToTraverse+1,2])

        ang = ang1
        i = 0 #counter
        v_index = i

        while np.sign(curvature)*(ang - ang2) <= 0:
            w = abs(vTraj[v_index][0])/r
            ang = ang + np.sign(curvature)*w*dT
            
            self.poses[i] = np.array([center[0] + r*cos(ang), center[1] + r*sin(ang), ang + np.sign(curvature)*np.pi/2.0])
            self.controls[i] = np.array([vTraj[v_index][0], curvature])
            
            i = i + 1
            if v_index < len(vTraj) - 1:
                v_index = v_index + 1

        self.poses = self.poses[0:i] #trim excess

class LineSegment(PathSegment):

    def __init__(self, speed, start, end, phi, dT):
        """
        Generates a line segment from start to end coordinates, with given direction, and speed profile. 

        Input: 
            curvaure: curvature of constant arc (float)
            vTraj: velocity profile of constant arc segment (Nx2 array of floats)
            start: SE2 position of start point (1x2 array of floats)
            end: SE2 position of end point (1x2 array of floats)
            phi: direction angle wrt to horizontal axis
            dT: timestep, must be timestep of overall planning profile (float)
        """
        
