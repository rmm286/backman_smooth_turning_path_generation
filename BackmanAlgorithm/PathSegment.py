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

    def __init__(self, kTrajectory, vTrajectory, xo, dT):
        self.poses = self.integrateTrajectory(kTrajectory, vTrajectory, xo, dT)
        self.controls = np.vstack([kTrajectory,vTrajectory])
    
    def integrateTrajectory(self, kTrajectory, vTrajectory, xo, dT):
        """ 
        Performs integration on dyanmic model using LSODA from odeint. kTrajetory and vTrajectory must be of equal length.

        Input:
            kTrajectory: Array of curvature values for each timestep
            vTrajectory: Array of velocity values for each timestep
            xo: Inital pose of model, x[0] = x, x[1] = y, x[2] = theta (orientation)

        Output:
            x: A Nx3 numpy array of [x,y,theta] poses for each timestep.
        """

        if (not (len(kTrajectory) == len(vTrajectory))):
            raise ValueError(
                "curvature and speed trajectories not same length")

        t = [i*dT for i in range(len(kTrajectory))]
        u = np.empty([len(kTrajectory), 2])
        for i in range(len(kTrajectory)):
            u[i][0] = vTrajectory[i]
            u[i][1] = kTrajectory[i]
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
            vTraj: velocity profile of constant arc segment (Nx1 array of floats)
            start: SE2 position of start point (1x3 array of floats)
            end: SE2 position of end point (1x3 array of floats)
            center: R2 position of center of turning (float)
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
        maxTimeToTraverse = arcLen/np.amin(vTraj)
        maxStepsToTraverse = int(maxTimeToTraverse/dT) +1

        self.poses = np.zeros([maxStepsToTraverse+1,3])
        self.controls = np.zeros([maxStepsToTraverse,2])

        ang = ang1
        i = 0 #counter
        v_index = i

        while np.sign(curvature)*(ang - ang2) <= 0:
            w = abs(vTraj[v_index])/r
            ang = ang + np.sign(curvature)*w*dT
            
            self.poses[i] = np.array([center[0] + r*cos(ang), center[1] + r*sin(ang), ang + np.sign(curvature)*np.pi/2.0])
            self.controls[i] = np.array([vTraj[v_index], curvature])
            
            i = i + 1
            if v_index < len(vTraj) - 1:
                v_index = v_index + 1

        self.poses = self.poses[0:i] #trim excess
        self.controls = self.controls[0:i] #trim excess


class C2ArcSegment(PathSegment):

    def __init__(self, curvature, v1, v2, start, end, center, reverse, dT):
        """
        Generates the center C2 segment from start to end coordinates, with given curvature, center of rotation and speed profile. 

        Input: 
            curvature: curvature of constant arc (float)
            v1: velocity profile of V_C21 arc segment (array of floats)
            v2: velocity profile of V_C22 arc segment (array of floats)
            start: SE2 position of start point (1x3 array of floats)
            end: SE2 position of end point (1x3 array of floats)
            center: R2 position of center of turning (1x2 of floats)
            dT: timestep, must be timestep of overall planning profile (float)
        """
        
        ang1 = np.arctan2(start[1]-center[1],start[0]-center[0])
        ang2 = np.arctan2(end[1]-center[1], end[0]-center[0])

        if reverse:
            ang1, ang2 = ang2, ang1
            
        r = np.linalg.norm(start[0:2]-center)

        if(curvature > 0 and ang2 < ang1):
            ang2 = ang2 + 2*np.pi
        
        if(curvature < 0 and ang2 > ang1):
            ang2 = ang2 - 2*np.pi

        arcLen = r*np.abs(ang2-ang1)
        maxTimeToTraverse = arcLen/np.amin(np.abs(np.append(v1,v2))) #TODO: change this to average
        maxStepsToTraverse = int(maxTimeToTraverse/dT) + 1

        self.poses = np.zeros([maxStepsToTraverse,3])
        self.controls = np.zeros([maxStepsToTraverse,2])

        i = 0 #counter
        i1 = 0
        i2 = maxStepsToTraverse - 1
        v_1_index = 0
        v_2_index = len(v2)-1

        while np.sign(curvature)*(ang1 - ang2) <= 0:
            
            if i%2 == 0: #forward step
                w = abs(v1[v_1_index])/r
                ang1 = ang1 + np.sign(curvature)*w*dT
                
                self.poses[i1] = np.array([center[0] + r*cos(ang1), center[1] + r*sin(ang1), ang1 + np.sign(curvature)*np.pi/2.0])
                self.controls[i1] = np.array([v1[v_1_index], curvature])
                
                i = i + 1
                i1 = i1 + 1
                if v_1_index < len(v1) - 1:
                    v_1_index = v_1_index + 1

            else: #back step
                w = abs(v2[v_2_index])/r
                ang2 = ang2 + -1*np.sign(curvature)*w*dT
                
                self.poses[i2] = np.array([center[0] + r*cos(ang2), center[1] + r*sin(ang2), ang2 + np.sign(curvature)*np.pi/2.0])
                self.controls[i2] = np.array([v2[v_2_index], curvature])
                
                i = i + 1
                i2 = i2 - 1
                if v_2_index > 0:
                    v_2_index = v_2_index - 1

        self.poses = np.append(self.poses[0:i1], self.poses[i2+1:len(self.poses)], axis = 0) #trim excess
        self.controls = np.append(self.controls[0:i1], self.controls[i2+1:len(self.controls)], axis = 0) #trim excess       

class C2LineSegment(PathSegment):

    def __init__(self, v1, v2, start, end, dT):
        """
        Generates the center C2 segment from start to end coordinates, with given curvature, center of rotation and speed profile. 

        Input: 
            curvature: curvature of constant arc (float)
            v1: velocity profile of V_C21 arc segment (array of floats)
            v2: velocity profile of V_C22 arc segment (array of floats)
            start: SE2 position of start point (1x3 array of floats)
            end: SE2 position of end point (1x3 array of floats)
            center: R2 position of center of turning (1x2 of floats)
            dT: timestep, must be timestep of overall planning profile (float)
        """
        segmentLength = np.linalg.norm(end[0:2]-start[0:2])
        maxTimeToTraverse = segmentLength/np.amin(np.append(v1,v2))
        maxStepsToTraverse = int(maxTimeToTraverse/dT) + 1

        self.poses = np.zeros([maxStepsToTraverse,3])
        self.controls = np.zeros([maxStepsToTraverse,2])

        pose1 = start
        pose2 = end
        v_1_index = 0
        v_2_index = len(v2) - 1
        i = 0
        i1 = 0
        i2 = maxStepsToTraverse - 1

        dirVec = (end[0:2]-start[0:2])/segmentLength
        tolerance = (dT*np.amax(np.append(v1,v2)))

        while np.linalg.norm(pose2[0:2] - pose1[0:2]) > tolerance:

            if i%2 == 0: #fwd
                pose1 = np.append(pose1[0:2] + np.array(dirVec*v1[v_1_index]*dT), start[2])
                self.poses[i1] = pose1
                self.controls[i1] = np.array([v1[v_1_index], 0.0])

                i = i + 1
                i1 = i1 + 1
                if v_1_index < len(v1) - 1:
                    v_1_index = v_1_index + 1

            else: #back step
                pose2 = np.append(pose2[0:2] - np.array(dirVec*v2[v_2_index]*dT), start[2])
                self.poses[i2] = pose2
                self.controls[i2] = np.array([v2[v_2_index], 0.0])
                
                i = i + 1
                i2 = i2 - 1
                if v_2_index > 0:
                    v_2_index = v_2_index - 1
        
        self.poses = np.append(self.poses[0:i1], self.poses[i2+1:len(self.poses)], axis = 0) #trim excess
        self.controls = np.append(self.controls[0:i1], self.controls[i2+1:len(self.controls)], axis = 0) #trim excess

class LineSegment(PathSegment):

    def __init__(self, vTraj, start, end, dT):
        """
        Generates a line segment from start to end coordinates, with given direction, and speed profile. 

        Input: 
            vTraj: velocity profile of constant arc segment (Nx2 array of floats)
            start: SE2 position of start point (1x2 array of floats)
            end: SE2 position of end point (1x2 array of floats)
            dT: timestep, must be timestep of overall planning profile (float)
        """
        
        segmentLength = np.linalg.norm(end[0:2]-start[0:2])
        maxTimeToTraverse = segmentLength/np.amin(vTraj)
        maxStepsToTraverse = int(maxTimeToTraverse/dT) + 1

        self.poses = np.zeros([maxStepsToTraverse,3])
        self.controls = np.zeros([maxStepsToTraverse,2])

        pose = start
        v_index = 0
        i = 0

        dirVec = (end[0:2]-start[0:2])/segmentLength
        tolerance = (dT*np.amax(vTraj))

        while np.linalg.norm(end[0:2] - pose[0:2]) > tolerance:
            
            pose = np.append(pose[0:2] + np.array(dirVec*vTraj[v_index]*dT), start[2]) 
    
            self.poses[i] = pose
            self.controls[i] = np.array([vTraj[v_index], 0.0])

            i = i + 1
            if v_index < len(vTraj) - 1:
                v_index = v_index + 1

        self.poses = self.poses[0:i] #trim excess
        self.controls = self.controls[0:i] #trim excess