#!/usr/bin/env python3
import numpy as np
from numpy import sin, cos, tan
import matplotlib.pyplot as plt
from PathSegment import SpiralSegment, CCSegment, LineSegment, C2ArcSegment, C2LineSegment, FullPath

class SmoothPathPlanner:
    """
    Class for implementation of Smooth curvature path generating algorithm as described in Juha Backman's 2015 paper "Smooth Turning PathGeneration for Agricultural Vehicles in Headlands".   
    
    Note that variable names are chosen either to conform to standard conventions or in some cases to reflect the naming conventions in the paper.
    """

    def __init__(self, dT):
        self.dT = dT

    def setConstraints(self, kConstraints, vConstraints, headlandSpeed, headlandSpeedReverse):
        """
        Set constraints on K, Kdot, Kddot, V, Vdot, Vddot, speed in headland and reverse speed in headland. 

        Input: 
            kConstraints: constraints on curvature and derivatives [kMax, kMin, kDotMax, kDotMin, kDDotMax, kDDotMin]
            vConstraints: constraints on velocity and derivatives [vMax, vMin, vDotMax, vDotMin, vDDotMax, vDDotMin]
            headlandSpeed: Desired speed in headland
            headlandSpeedReverse: Desired reverse speed in headland
        """

        if headlandSpeed > vConstraints[0]:
            raise ValueError("Headland Speed should not be larger than V Max")

        if headlandSpeedReverse > vConstraints[1]:
            raise ValueError("Headland Speed in reverse should not be smaller than V Min")

        if not (len(kConstraints) == 6 and len(vConstraints) == 6):
            raise TypeError('kConstraintsand/or vConstraints not of correct dimension.')
        else:
            self.kMax = kConstraints[0]
            self.kMin = kConstraints[1]
            self.kDotMax = kConstraints[2]
            self.kDotMin = kConstraints[3]
            self.kDDotMax = kConstraints[4]
            self.kDDotMin = kConstraints[5]

            self.vMax = vConstraints[0]
            self.vMin = vConstraints[1]
            self.vDotMax = vConstraints[2]
            self.vDotMin = vConstraints[3]
            self.vDDotMax = vConstraints[4]
            self.vDDotMin = vConstraints[5]

        self.headlandSpeed = headlandSpeed
        self.headlandSpeedReverse = headlandSpeedReverse

    def setStartAndGoal(self, initialState, finalState):
        """
        Set Start and Goal States.

        Input:
            initialState: start state of vehicle [x, y, theta, v, k]
            finalState: end or goal state of vehicle [x, y, theta, v, k]
        """

        if not (len(initialState) == 5 and len(finalState) == 5):
            raise TypeError('Start and/or Goal not of correct dimension: [x, y, th, v, k]')
        else:
            self.initialState = initialState
            self.finalState = finalState
            self.k_C0 = initialState[4]
            self.k_C4 = finalState[4]

    def setNominalCurvatures(self, kStart, kCenter, kEnd, reverse):
        """ 
        Set the curvature for each constant curvature segment.
            
        Input:
            kStart: Curvature of first constant curvature segment
            kCenter: Curvature of middle constant curvature segment
            kEnd: Curvature of final constant curvature segment
            reverse: Whether the middle segment is driven with reverse speed
        """

        self.k_C1 = kStart
        self.k_C2 = kCenter
        self.k_C3 = kEnd
        self.reverse = reverse

    def generateCurvatureTrajectory(self, kInit, kFinal):
        """
        Generates a curvature trajectory which has starting curvature equal to kInit and final curvature equal to kFinal. Time values start form tInit and increment by self.dT
        Uses Eqn. 4 from paper.

        Input:  
            kInit: starting curvature of trajectory
            kFinal: final curvature of trajectory
        
        Output:
            kTrajectory: np array of curvature values for ever timestep in trajectory 
        """

        kTolerance = self.dT*max(np.abs(self.kDotMax),np.abs(self.kDotMin))
        k = kInit
        kTrajectory = np.array([k])
        while np.abs(k - kFinal) > kTolerance:
            if k < kFinal:
                k = k + self.dT*self.kDotMax
            else:
                k = k + self.dT*self.kDotMin
            kTrajectory = np.append(kTrajectory, np.array([k]), axis=0)

        return kTrajectory

    def generateSpeedTrajectory(self, vInit, vFinal):
        """
        Generates a velocity trajectory which has starting velocity equal to vInit and final velocity equal to vFinal. Time values start form tInit and increment by self.dT
        Eqn. 7 from paper

        Input:  
            vInit: starting velocity of trajectory
            vFinal: final velocity of trajectory
        
        Output:
            vTrajectory: np array of velocity values for ever timestep in trajectory
        """

        vTolerance = self.dT*max(np.abs(self.vDotMax), np.abs(self.vDotMin))
        v = vInit
        vTrajectory = np.array([v])

        while np.abs(v - vFinal) > vTolerance:
            if v < vFinal:
                v = v + self.dT*self.vDotMax
            else:
                v = v + self.dT*self.vDotMin
            vTrajectory = np.append(vTrajectory, np.array([v]), axis=0)

        return vTrajectory

    def generateOptimalTrajectory(self, x0, xFinal, xDotMax, xDotMin, xDDotMax, xDDotMin):
        """
        Analytically solves the optimal trajectory problem to find the x trajectory which moves from x0 to xFinal in minimum time subject to arbitrary boundary constraints on the first and second derivative of x. Returns optimal trajectory as a 1xN numpy array of floats.
        
        Input:
            x0: Initial state
            xFinal: Final state
            xDotMax: Upper bound on first derivative of x
            xDotMin: Lower bound on first derivative of x
            xDDotMax: Upper bound on second derivative of x
            xDDotMin: Lower bound on second derivative of x (This value must be negative)
        
        Output:
            trajectory: Nx1 np array of x values for every timestep.
        """

        dT = self.dT

        if x0 < xFinal: #increase x to xf
            riseTime = (xDotMax)/(xDDotMax)
            fallTime = (-1*xDotMax)/(xDDotMin)
            xRT = x0 + (xDDotMax/2.0)*riseTime**2
            xFT = xRT + (xDDotMax*riseTime)*fallTime + (xDDotMin/2.0)*fallTime**2

            if xFT < xFinal:
                diff = xFinal - xFT
                tTop = diff/xDotMax
                t = np.linspace(0,riseTime,int(riseTime/dT))
                x0toRT = x0 + (xDDotMax/2.0)*t**2
                t = np.linspace(dT, tTop, int(tTop/dT))
                xRTtoDiff = x0toRT[-1] + xDDotMax*riseTime*t
                t = np.linspace(dT,fallTime,int(fallTime/dT))
                xDifftoFT = xRTtoDiff[-1] + xDDotMax*riseTime*t + 0.5*xDDotMin*t**2
                xTrajectory = np.append(x0toRT,np.append(xRTtoDiff,xDifftoFT,axis=0),axis=0)
            else:
                t1 = np.sqrt((xFinal - x0)/((xDDotMax/2.0)*(1-xDDotMax/xDDotMin)))
                t2 = -1*(xDDotMax/xDDotMin)*t1
                t = np.linspace(0,t1,int(t1/dT))
                x0tot1 = x0 + xDDotMax/2.0*t**2
                t = np.linspace(dT,t2,int(t2/dT))
                xt1tot2 = x0tot1[-1] + (xDDotMax*t1)*t + (xDDotMin/2.0)*t**2
                xTrajectory = np.append(x0tot1,xt1tot2,axis=0)
        elif x0 > xFinal: #decrease x to xf
            fallTime = (xDotMin)/(xDDotMin)
            riseTime = (-1*xDotMin)/(xDDotMax)
            xFT = x0 + (xDDotMin/2.0)*fallTime**2
            xRT = xFT + (xDDotMin*fallTime)*riseTime + (xDDotMax/2.0)*riseTime**2

            if xRT > xFinal:
                diff = xFinal - xRT
                tBottom = diff/xDotMin
                t = np.linspace(0,fallTime,int(fallTime/dT))
                x0toFT = x0 + (xDDotMin/2.0)*t**2
                t = np.linspace(dT, tBottom, int(tBottom/dT))
                xFTtoDiff = x0toFT[-1] + xDDotMin*fallTime*t
                t = np.linspace(dT,riseTime,int(riseTime/dT))
                xDifftoRT = xFTtoDiff[-1] + xDDotMin*fallTime*t + 0.5*xDDotMax*t**2
                xTrajectory = np.append(x0toFT,np.append(xFTtoDiff,xDifftoRT,axis=0),axis=0)
            else:
                t1 = np.sqrt((xFinal - x0)/((xDDotMin/2.0)*(1-xDDotMin/xDDotMax)))
                t2 = -1*(xDDotMin/xDDotMax)*t1
                t = np.linspace(0,t1,int(t1/dT))
                x0tot1 = x0 + xDDotMin/2.0*t**2
                t = np.linspace(dT,t2,int(t2/dT))
                xt1tot2 = x0tot1[-1] + (xDDotMin*t1)*t + (xDDotMax/2.0)*t**2
                xTrajectory = np.append(x0tot1,xt1tot2,axis=0)
        else: #x0 = xFinal
            xTrajectory = np.array([x0, xFinal])
        
        return xTrajectory

    def generateCurvatureTrajectoryDDot(self, k0, kFinal):
        """
        Helper function to call generateOptimalTrajectory and return a trajectory with curvature equal to kInit and final curvature equal to kFinal.
        This function differs from generateCurvatureTrajectoryDDot in that it produces a transition which is both continuous and differentiable, and respects kDDot constraints.
        
        Input:
            k0: Initial curvature
            kFinal: Final curvature
        
        Output:
            kTrajectory: 1xN numpy array of curvature values for each timestep
        """
        xDotMax = self.kDotMax
        xDotMin = self.kDotMin
        xDDotMax = self.kDDotMax
        xDDotMin = self.kDDotMin
        
        return self.generateOptimalTrajectory(k0, kFinal, xDotMax, xDotMin, xDDotMax, xDDotMin)

    def generateSpeedTrajectoryDDot(self, v0, vFinal):
        """
        Helper function to call generateOptimalTrajectory and return a speed trajectory which has starting velocity equal to vInit and final speed equal to vFinal.
        
        Input:
            v0: Initial speed
            vFinal: Final speed

        Output:
            vTrajectory: 1xN numpy array of velocity values
        """
        xDotMax = self.vDotMax
        xDotMin = self.vDotMin
        xDDotMax = self.vDDotMax
        xDDotMin = self.vDDotMin

        return self.generateOptimalTrajectory(v0, vFinal, xDotMax, xDotMin, xDDotMax, xDDotMin)

    def makeTrajectoriesEqualLength(self, kTrajIn, vTrajIn, fromStart=False):
        """
        Takes curvature and velocity trajectories and makes them the same length. Either cutting vTraj from the start or end, or by adding values to the start or end of vTraj.
        
        Input:
            kTrajIn: Input curvature trajectory
            vTrajIn: input velocity trajectory
            fromStart: set to true if values will be cut or added to the front of the array, False if values are cut or added to end of the array
        
        Output: 
            Output Dict:
                kTraj: return curvature trajectory
                vTraj: return velocity trajectory
                cutV: bool value, true if the velocity trajectory was shortened by the operation
                leftover: array of velocity values that was cut from input array, if cutV is False then this is a 2x1 array of either the first or last value of the input velocity array
        """

        cutV = False
        if (len(kTrajIn) < len(vTrajIn)) and not fromStart:  # cut from end of vTraj
            vTraj = vTrajIn[0:len(kTrajIn)]
            leftover = vTrajIn[len(kTrajIn):len(vTrajIn)]
            cutV = True
        elif (len(kTrajIn) < len(vTrajIn)) and fromStart:  # cut from start of vTraj
            vTraj = vTrajIn[len(vTrajIn) - len(kTrajIn):len(vTrajIn)]
            leftover = vTrajIn[0:len(vTrajIn) - len(kTrajIn)]
            cutV = True
        elif (len(kTrajIn) > len(vTrajIn)) and not fromStart: # add to end of vTraj
            extension = vTrajIn[-1]*np.ones(len(kTrajIn) - len(vTrajIn))
            vTraj = np.append(vTrajIn, extension, axis=0)
            leftover = vTrajIn[-1]*np.ones(2)
        elif (len(kTrajIn) > len(vTrajIn)) and fromStart: # add to start of vTraj
            extension = vTrajIn[0]*np.ones(len(kTrajIn) - len(vTrajIn))
            vTraj = np.append(extension, vTrajIn, axis = 0)
            leftover = vTrajIn[0]*np.ones(2)
        elif (len(kTrajIn) == len(vTrajIn)) and not fromStart:
            leftover = vTrajIn[-1]*np.ones(2)
            vTraj = vTrajIn
        elif (len(kTrajIn) == len(vTrajIn)) and fromStart:
            leftover = vTrajIn[0]*np.ones(2)
            vTraj = vTrajIn

        return {'kTraj': kTrajIn, 'vTraj': vTraj, 'cutV': cutV, 'leftover': leftover }

    def calculateCenterLine(self):
        """
        Uses equations (17) through (21) of paper to calculate orientation of the center line.

        Output:
            Phi: The direction of the center line represented as a radian value
        """
        S1 = self.S1.poses
        S2 = self.S2.poses
        S3 = self.S3.poses
        S4 = self.S4.poses
        k_C1 = self.k_C1
        k_C3 = self.k_C3

        omega_S2_tS2 = np.array(
            [S2[0][0] - (k_C1**-1)*sin(S2[0][2]), S2[0][1] + (k_C1**-1)*cos(S2[0][2])]) #instantaneous center of turning at beginning of spiral segment S2.
        omega_S3_tC3 = np.array(
            [S3[-1][0] - (k_C3**-1)*sin(S3[-1][2]), S3[-1][1] + (k_C3**-1)*cos(S3[-1][2])]) #instantaneous center of turning at end of spiral segment S3.
        #Note here that the spiral segments S2 and S3 have not been placed yet, so these centers can tell us relative displacement, they cannot give us position of the center (C2) CC segment. 

        self.omega_k = np.array(
            [S1[-1][0] - (k_C1**-1)*sin(S1[-1][2]), S1[-1][1] + (k_C1**-1)*cos(S1[-1][2])]) #center of turning of C1
        self.omega_kplus2 = np.array(
            [S4[0][0] - (k_C3**-1)*sin(S4[0][2]), S4[0][1] + (k_C3**-1)*cos(S4[0][2])]) #center of turning of C3
        rk = cos(S2[-1][2]) * (omega_S2_tS2[1] - S2[-1][1]) - sin(S2[-1][2]) * (omega_S2_tS2[0] - S2[-1][0])
        rkplus1 = cos(S3[0][2]) * (omega_S3_tC3[1] - S3[0][1]) - sin(S3[0][2]) * (omega_S3_tC3[0] - S3[0][0])
        d = np.linalg.norm(self.omega_kplus2 - self.omega_k)
        
        if np.abs(rkplus1 - rk) > d:
            raise Exception

        w = np.arcsin((rkplus1 - rk)/d)
        vec = self.omega_kplus2 - self.omega_k

        return np.arctan2(vec[1], vec[0]) - w #phi

    def calculateConstantArcs(self):
        """
        Calculates the positions of the center of turning for each constant curvature segments. Uses equations (10) through (16) from paper.

        Follow the following conventions from paper: 
            Omega_Ck is the center of turning for the Kth constant curvature segment
            Omega_Sk_tSk is the instantaneous center of turning for the Kth spiral segment at the beginning of the segment Sk.
            Omega_Sk_tCk+1 is the  instantaneous center of turning for the Kth spiral segment at the end of the segment Sk. This is also equal to the center of turning for the Ck+1 CC segment.

        Output: 
            omega_kplus1: center of turning for the center arc.
        """
        S1 = self.S1.poses
        S2 = self.S2.poses
        S3 = self.S3.poses
        S4 = self.S4.poses
        k_C1 = self.k_C1
        k_C2 = self.k_C2
        k_C3 = self.k_C3

        omega_S2_tC2 = np.array(
            [S2[-1][0] - (k_C2**-1)*sin(S2[-1][2]), S2[-1][1] + (k_C2**-1)*cos(S2[-1][2])]) #instantaneous center of turning at end of spiral segment S2.
        omega_S2_tS2 = np.array(
            [S2[0][0] - (k_C1**-1)*sin(S2[0][2]), S2[0][1] + (k_C1**-1)*cos(S2[0][2])]) #instantaneous center of turning at beginning of spiral segment S2.

        omega_S3_tC3 = np.array(
            [S3[-1][0] - (k_C3**-1)*sin(S3[-1][2]), S3[-1][1] + (k_C3**-1)*cos(S3[-1][2])]) #instantaneous center of turning at end of spiral segment S3.
        omega_S3_tS3 = np.array(
            [S3[0][0] - (k_C2**-1)*sin(S3[0][2]), S3[0][1] + (k_C2**-1)*cos(S3[0][2])]) #instantaneous center of turning at start of spiral segment S3.

        #Note here that the spiral segments S2 and S3 have not been placed yet, so these centers can tell us relative displacement, they cannot give us position of the center (C2) CC segment. 

        self.omega_k = np.array([S1[-1][0] - (k_C1**-1)*sin(S1[-1][2]), S1[-1][1] + (k_C1**-1)*cos(S1[-1][2])]) #center of turning of C1
        self.omega_kplus2 = np.array([S4[0][0] - (k_C3**-1)*sin(S4[0][2]), S4[0][1] + (k_C3**-1)*cos(S4[0][2])]) #center of turning of C3
        d1 = np.linalg.norm(self.omega_kplus2 - self.omega_k)
        d2 = np.linalg.norm(omega_S2_tC2 - omega_S2_tS2)
        d3 = np.linalg.norm(omega_S3_tC3 - omega_S3_tS3)
        l1 = (d2**2 - d3**2 + d1**2)/(2*d1)
        if(l1 > d2):
            raise Exception
        l2 = np.sqrt(d2**2 - l1**2)
        signVal = 1 if (k_C2 < 0) else -1
        signVal = -signVal if (self.omega_kplus2[0] > self.omega_k[0] and k_C2 > 0) else signVal
        #signVal = -signVal if self.crossover else signVal
        self.omega_kplus1 = self.omega_k + l1*(self.omega_kplus2 - self.omega_k)/d1 + signVal*l2*np.matmul(
            np.array([[0, -1], [1, 0]]), (self.omega_kplus2 - self.omega_k))/d1 #center of turning of C2

        return self.omega_kplus1

    def plotPaths(self, plotCircles = False, plotArrows = False):
        """
        Plots all segments of the paths. Can optionally plot the turning circles of constant curvature segments and orientation data.

        Input:
            plotCircles (optional): True if constant curvature circles are to be plotted
            plotArrows (optional): True if orientation data are to be plotted
        
        Output:
            Saves plot to file "trajectory.png"
        """

        plt.figure(0)
        plt.clf()
        plt.title("Path Segments")
        plt.xlabel("x (m)")
        plt.ylabel("y (m)")
        
        if hasattr(self, 'S1'):
            plt.plot([i[0] for i in self.S1.poses], [i[1] for i in self.S1.poses])
        if hasattr(self, 'S2'):    
            plt.plot([i[0] for i in self.S2.poses], [i[1] for i in self.S2.poses])
        if hasattr(self, 'S3'):
            plt.plot([i[0] for i in self.S3.poses], [i[1] for i in self.S3.poses])
        if hasattr(self, 'S4'):
            plt.plot([i[0] for i in self.S4.poses], [i[1] for i in self.S4.poses])
        if hasattr(self, 'C1'):
            plt.plot([i[0] for i in self.C1.poses], [i[1] for i in self.C1.poses])
        if hasattr(self, 'C2'):
            plt.plot([i[0] for i in self.C2.poses], [i[1] for i in self.C2.poses])
        if hasattr(self, 'C3'):
            plt.plot([i[0] for i in self.C3.poses], [i[1] for i in self.C3.poses])

        if plotCircles: 
            plt.plot(self.omega_k[0], self.omega_k[1], 'b^')
            if hasattr(self,'omega_kplus1'):
                plt.plot(self.omega_kplus1[0], self.omega_kplus1[1], 'b^')

            plt.plot(self.omega_kplus2[0], self.omega_kplus2[1], 'b^')

            plt.plot([(self.k_C1**-1)*cos(theta) + self.omega_k[0] for theta in np.linspace(0, 2*np.pi, 25)],[(self.k_C1**-1)*sin(theta) + self.omega_k[1] for theta in np.linspace(0, 2*np.pi, 25)], 'r--')
            if hasattr(self,'omega_kplus1'):
                plt.plot([(self.k_C2**-1)*cos(theta) + self.omega_kplus1[0] for theta in np.linspace(0, 2*np.pi, 25)],[(self.k_C2**-1)*sin(theta) + self.omega_kplus1[1] for theta in np.linspace(0, 2*np.pi, 25)], 'r--')
                        
            plt.plot([(self.k_C3**-1)*cos(theta) + self.omega_kplus2[0] for theta in np.linspace(0, 2*np.pi, 25)], [(self.k_C3**-1)*sin(theta) + self.omega_kplus2[1] for theta in np.linspace(0, 2*np.pi, 25)], 'r--')
            plt.arrow(self.S1.poses[-1][0], self.S1.poses[-1][1], 0.1*cos(self.S1.poses[-1][2]), 0.1*sin(self.S1.poses[-1][2]), length_includes_head = True, width = 0.02, head_width = 0.03, color = 'r', alpha = 0.5)
        
        if plotArrows:
            for i in range(0, len(self.S2.poses), int(len(self.S2.poses)/10)):
                plt.arrow(self.S2.poses[i][0], self.S2.poses[i][1], 0.1*cos(self.S2.poses[i][2]), 0.1*sin(self.S2.poses[i][2]), length_includes_head = True, width = 0.01, head_width = 0.03, color = 'r', alpha = 0.5)
            
            for i in range(0, len(self.S4.poses), int(len(self.S2.poses)/10)):
                plt.arrow(self.S4.poses[i][0], self.S4.poses[i][1], 0.1*cos(self.S4.poses[i][2]), 0.1*sin(self.S4.poses[i][2]), length_includes_head = True, width = 0.01, head_width = 0.03, color = 'r', alpha = 0.5)

        plt.savefig("trajectory.png")

    def plotControls(self):
        """
        Plots all control inputs for each segment. Can only be called once full path has been created.

        Output:
            Saves to two files "velProfile.png" and "kProfile.png"
        """
        plt.figure(1)
        plt.clf()
        plt.title("Speed Profile")
        plt.xlabel("time (s)")
        plt.ylabel("speed (m/s)")

        timeseries_S1 = np.array([i*self.dT for i in range(len(self.S1.controls))])
        timeseries_C1 = np.array([i*self.dT for i in range(len(self.C1.controls))]) + timeseries_S1[-1]
        timeseries_S2 = np.array([i*self.dT for i in range(len(self.S2.controls))]) + timeseries_C1[-1]
        timeseries_C2 = np.array([i*self.dT for i in range(len(self.C2.controls))]) + timeseries_S2[-1]
        timeseries_S3 = np.array([i*self.dT for i in range(len(self.S3.controls))]) + timeseries_C2[-1]
        timeseries_C3 = np.array([i*self.dT for i in range(len(self.C3.controls))]) + timeseries_S3[-1]
        timeseries_S4 = np.array([i*self.dT for i in range(len(self.S4.controls))]) + timeseries_C3[-1]

        plt.plot(timeseries_S1, [i[0] for i in self.S1.controls])
        plt.plot(timeseries_C1, [i[0] for i in self.C1.controls])
        plt.plot(timeseries_S2, [i[0] for i in self.S2.controls])
        plt.plot(timeseries_C2, [i[0] for i in self.C2.controls])
        plt.plot(timeseries_S3, [i[0] for i in self.S3.controls])
        plt.plot(timeseries_C3, [i[0] for i in self.C3.controls])
        plt.plot(timeseries_S4, [i[0] for i in self.S4.controls])

        plt.savefig("velProfile.png")

        plt.figure(2)
        plt.clf()
        plt.title("Curvature Profile")
        plt.xlabel("time (s)")
        plt.ylabel("curvature (m^-1)")

        plt.plot(timeseries_S1, [i[1] for i in self.S1.controls])
        plt.plot(timeseries_C1, [i[1] for i in self.C1.controls])
        plt.plot(timeseries_S2, [i[1] for i in self.S2.controls])
        plt.plot(timeseries_C2, [i[1] for i in self.C2.controls])
        plt.plot(timeseries_S3, [i[1] for i in self.S3.controls])
        plt.plot(timeseries_C3, [i[1] for i in self.C3.controls])
        plt.plot(timeseries_S4, [i[1] for i in self.S4.controls])

        plt.savefig("kProfile.png")

    def plan(self):
        """
        Central planner function for the object. Returns a path from start to end state that conforms to all constraints
        
        Output:
            FullPath: Full Path object
        """

        if not hasattr(self, 'headlandSpeed'):
            raise ValueError('setConstraints has not been called successfully on this instance.')

        if not hasattr(self, 'initialState'):
            raise ValueError('setStartAndGoal has not been called successfully on this instance.')

        if not hasattr(self, 'k_C1'):
            raise ValueError('setNominalCurvatures has not been called successfully on this instance.')

        self.path_is_not_feasible = True

        while self.path_is_not_feasible:

            ################ generate first connecting spiral ################
            K_S1 = self.generateCurvatureTrajectoryDDot(
                self.k_C0, self.k_C1)
            v_S1 = self.generateSpeedTrajectoryDDot(
                self.initialState[3], self.headlandSpeed)
            trajectories = self.makeTrajectoriesEqualLength(K_S1, v_S1, False)
            self.cut_v_S1 = trajectories['cutV']
            v_C1 = trajectories['leftover']
            xo = [self.initialState[0], self.initialState[1],self.initialState[2]]
            self.S1 = SpiralSegment(trajectories['kTraj'], trajectories['vTraj'], xo, self.dT)
            
            ################ generate last connecting spiral ################
            K_S4 = self.generateCurvatureTrajectoryDDot(self.k_C3, self.k_C4)
            v_S4 = self.generateSpeedTrajectoryDDot(self.headlandSpeed, self.finalState[3])
            trajectories = self.makeTrajectoriesEqualLength(K_S4, v_S4, True)
            self.cut_v_S4 = trajectories['cutV']
            v_C3 = trajectories['leftover']
            xo = [0, 0, 0]
            self.S4 = SpiralSegment(
                trajectories['kTraj'], trajectories['vTraj'], xo, self.dT)

            self.S4.placePath(self.finalState[0], self.finalState[1], self.finalState[2], False)

            ################ generate second connecting spiral ################
            
            xo = [0, 0, 0]
            if self.cut_v_S1:
                if not(hasattr(self,'vC1_index')): 
                    #initial guess: C1 is long enough that vC1 is used completely
                    self.vC1_index = len(v_C1) - 1
            else: #know VC1 and VC2
                self.vC1_index = len(v_C1)  - 1
            
            if self.vC1_index < 0 or self.vC1_index > len(v_C1) - 1:
                    raise ValueError("vC1_index is out of range")

            if self.reverse:
                
                v_S21 = self.generateSpeedTrajectoryDDot(v_C1[self.vC1_index], 0) #max function is in case vC1_index = 0
                K_S21 = self.generateCurvatureTrajectoryDDot(self.k_C1, 0)

                if len(v_S21) > len(K_S21):
                    diff = len(v_S21) - len(K_S21)
                    K_S21 = np.append(K_S21, np.zeros(diff))
                elif len(v_S21) < len(K_S21):
                    diff = len(K_S21) - len(v_S21)
                    v_S21 = np.append(v_S21, np.zeros(diff))

                v_S22 = self.generateSpeedTrajectoryDDot(0, self.headlandSpeedReverse)
                K_S22 = self.generateCurvatureTrajectoryDDot(0, self.k_C2)
                trajectories = self.makeTrajectoriesEqualLength(K_S22,v_S22, False)
                self.cut_v_S2 = trajectories['cutV']
                v_C21 = trajectories['leftover']
                v_S2 = np.append(v_S21, trajectories['vTraj'], axis = 0)
                K_S2 = np.append(K_S21, trajectories['kTraj'], axis = 0)
                self.S2 = SpiralSegment(K_S2, v_S2, xo, self.dT)
            else:
                K_S2 = self.generateCurvatureTrajectoryDDot(self.k_C1, self.k_C2)
                v_S2 = v_C1[self.vC1_index:len(v_C1)]
                trajectories = self.makeTrajectoriesEqualLength(K_S2, v_S2, False)
                self.cut_v_S2 = trajectories['cutV']
                v_C21 = trajectories['leftover']
                self.S2 = SpiralSegment(trajectories['kTraj'], trajectories['vTraj'], xo, self.dT)

            ################ generate third connecting spiral ################
            
            if self.cut_v_S4:
                if not(hasattr(self,'vC3_index')):
                    #initial guess: C3 is long enough that vC3 is used completely
                    self.vC3_index = 0
            else: #know value of v_C4
                self.vC3_index = 0

            if self.vC3_index < 0 or self.vC3_index > len(v_C3)-1:
                raise ValueError("vC3_index is out of range.")
            
            if self.reverse:
                v_S31 = self.generateSpeedTrajectoryDDot(self.headlandSpeedReverse, 0)
                K_S31 = self.generateCurvatureTrajectoryDDot(self.k_C2, 0)
                trajectories = self.makeTrajectoriesEqualLength(K_S31, v_S31, True)
                self.cut_v_S3 = trajectories['cutV']
                v_C22 = trajectories['leftover']
                v_S32 = self.generateSpeedTrajectoryDDot(0, v_C3[self.vC3_index])
                K_S32 = self.generateCurvatureTrajectoryDDot(0, self.k_C3)
                
                if len(v_S32) > len(K_S32):
                    diff = len(v_S32) - len(K_S32)
                    K_S32 = np.append(np.zeros(diff), K_S32, axis = 0)
                elif len(v_S32) < len(K_S32):
                    diff = len(K_S32) - len(v_S32)
                    v_S32 = np.append(np.zeros(diff),v_S32, axis = 0 )

                v_S3 = np.append(trajectories['vTraj'], v_S32, axis = 0)
                K_S3 = np.append(trajectories['kTraj'], K_S32, axis = 0)

                self.S3 = SpiralSegment(K_S3, v_S3, xo, self.dT)
            else:

                K_S3 = self.generateCurvatureTrajectoryDDot(self.k_C2, self.k_C3)
                v_S3 = self.generateSpeedTrajectoryDDot(self.headlandSpeed, v_C3[self.vC3_index])
                trajectories = self.makeTrajectoriesEqualLength(K_S3, v_S3, True)
                self.cut_v_S3 = trajectories['cutV']
                v_C22 = trajectories['leftover']
                self.S3 = SpiralSegment(trajectories['kTraj'], trajectories['vTraj'], xo, self.dT)

            ################ generate center CC segment ################
            if np.abs(self.k_C2) > 0.005:  # center is an arc
                
                try:
                    self.calculateConstantArcs()
                except: #can't use this turning type
                    return 0 
                
                ################ place S2 and S3 segments ################
                
                self.S2.placePath(self.S1.poses[-1][0], self.S1.poses[-1][1], self.S1.poses[-1][2])
                omega_S2_tC2 = np.array([self.S2.poses[-1][0] - self.k_C2**-1 *sin(self.S2.poses[-1][2]), self.S2.poses[-1][1] + self.k_C2**-1*cos(self.S2.poses[-1][2])])
                r = np.linalg.norm(omega_S2_tC2 - self.omega_k)
                rotAngle = np.arccos((self.omega_kplus1[0] - self.omega_k[0])/r) - np.arccos((omega_S2_tC2[0] - self.omega_k[0]) /r)
                self.S2.rotateAboutPoint(self.omega_k[0], self.omega_k[1], rotAngle)
                self.S3.placePath(self.S4.poses[0][0], self.S4.poses[0][1], self.S4.poses[0][2], False)
                omega_S3_tS3 = np.array([self.S3.poses[0][0] - self.k_C2**-1 *sin(self.S3.poses[0][2]), self.S3.poses[0][1] + self.k_C2**-1*cos(self.S3.poses[0][2])])
                r = np.linalg.norm(omega_S3_tS3 - self.omega_kplus2)
                rotAngle =  np.arccos((self.omega_kplus1[0] - self.omega_kplus2[0])/r) - np.arccos((omega_S3_tS3[0]- self.omega_kplus2[0])/r)
                self.S3.rotateAboutPoint(self.omega_kplus2[0], self.omega_kplus2[1], rotAngle)

                ################ make C2 segment ################
                self.C2 = C2ArcSegment(self.k_C2, v_C21, v_C22, self.S2.poses[-1], self.S3.poses[0], self.omega_kplus1, self.reverse, self.dT)

            else: #center is a line

                ################ place S2 and S3 segments ################
                try:
                    self.phi_C2 = self.calculateCenterLine()
                except: #can't use this turning type
                    return 0
                
                self.S2.placePath(self.S1.poses[-1][0], self.S1.poses[-1][1], self.S1.poses[-1][2])
                rotAngle = self.phi_C2 - self.S2.poses[-1][2]
                self.S2.rotateAboutPoint(self.omega_k[0], self.omega_k[1], rotAngle)
                self.S3.placePath(self.S4.poses[0][0], self.S4.poses[0][1], self.S4.poses[0][2], False)
                rotAngle = self.phi_C2 - self.S3.poses[0][2]
                self.S3.rotateAboutPoint(self.omega_kplus2[0], self.omega_kplus2[1], rotAngle)

                if self.S2.pathIntersectsWith(self.S3):
                    return 0
                    #TODO: implement fix for this

                ################ make C2 segment ################
                self.C2 = C2LineSegment(v_C21, v_C22, self.S2.poses[-1], self.S3.poses[0], self.dT)

            ################ generate C1 and C3 ################
            self.C1 = CCSegment(self.k_C1, v_C1, self.S1.poses[-1], self.S2.poses[0], self.omega_k, self.dT)
            
            self.C3 = CCSegment(self.k_C3, v_C3[min(self.vC3_index,len(v_C3)-1):len(v_C3)], self.S3.poses[-1], self.S4.poses[0], self.omega_kplus2, self.dT)

            ################ Check Feasibility ################
            self.path_is_not_feasible = False

            if np.abs(self.C1.angle) > 1.5*np.pi:
                self.k_C1 = self.k_C1 - np.sign(self.k_C1)*0.1*self.kMax
                if np.abs(self.k_C1) < 10**-10:
                    return 0
                
                self.path_is_not_feasible = True
                continue

            if np.abs(self.C3.angle) > 1.5*np.pi:
                self.k_C3 = self.k_C3 - np.sign(self.k_C3)*0.1*self.kMax
                if np.abs(self.k_C3) <= 10**-10:
                    return 0 

                self.path_is_not_feasible = True
                continue
            
            if not(hasattr(self,'vC1_index_new')):
                self.vC1_index_new = self.C1.v_index
                if (self.vC1_index != self.vC1_index_new):
                    self.vC1_index = self.vC1_index_new
                    self.path_is_not_feasible = True

            if not(hasattr(self,'vC3_index_new')):
                self.vC3_index_new = (len(v_C3) - 1) - self.C3.v_index
                if (self.vC3_index != self.vC3_index_new):
                    self.vC3_index = self.vC3_index_new
                    self.path_is_not_feasible = True
            
            if (self.C1.arcLen < self.dT*self.vMax) and (np.linalg.norm(self.S2.poses[0][0:2] - self.S1.poses[-1][0:2]) > self.dT*self.vMax):
                return 0

            if (self.C3.arcLen < self.dT*self.vMax) and (np.linalg.norm(self.S4.poses[0][0:2] - self.S3.poses[-1][0:2]) > self.dT*self.vMax):
                return 0

        return FullPath(self.dT, self.S1, self.C1, self.S2, self.C2, self.S3, self.C3, self.S4)

    def planShortest(self):
        """
        Generic path planner which returns shortest path from start to goal.

        Output:
            FullPath: Full Path object
        """
        RSR = [self.kMin, 0, self.kMin, False]
        LSL = [self.kMax, 0, self.kMax, False]
        LRL = [self.kMax, self.kMin, self.kMax, False]
        RLR = [self.kMin, self.kMax, self.kMin, False]
        LSR = [self.kMax, 0, self.kMin, False]
        RSL = [self.kMin, 0, self.kMax, False]
        R1L1R = [self.kMin, self.kMax, self.kMin, True]
        L1R1L = [self.kMax, self.kMin, self.kMax, True]

        turningTypes = [RSR,LSL,LRL,RLR,LSR,RSL,R1L1R,L1R1L]

        shortestPathFinalTime = np.inf

        for i in range(len(turningTypes)):
            #turningTypeNames = ['RSR','LSL','LRL','RLR','LSR','RSL','R1L1R','L1R1L']
            pathType = turningTypes[i]
            self.setNominalCurvatures(pathType[0], pathType[1], pathType[2], pathType[3])
            path = self.plan()

            if path == 0:
                continue
            elif path.finalTime < shortestPathFinalTime:
                shortestPath = path
                shortestPathFinalTime = path.finalTime

        return shortestPath

# make variable names consistent
# make spacing consistent
# make sure tolerances make sense, maybe pass tolerances as parameters