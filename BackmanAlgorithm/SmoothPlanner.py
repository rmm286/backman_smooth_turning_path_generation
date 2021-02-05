import numpy as np
from numpy import sin, cos, tan
from ODESolver import RungeKutta4
from StateSpaces import SmoothPathState
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from PathSegment import SpiralSegment, CCSegment, LineSegment, C2ArcSegment, C2LineSegment

class SmoothPathPlanner:
    """ class for implementation of backman algorithm"""

    def __init__(self):

        self.path = []

    def setConstraints(self, kConstraints, vConstraints, headlandSpeed, headlandSpeedReverse):
        """ Set constraints on K, Kdot, Kddot, V, Vdot, Vddot, speed in headland and reverse speed in headland. """

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

    def setStartAndGoal(self, initalState, finalState):
        """ Set Start and Goal States, both states must be SE2States """

        if not (isinstance(initalState, SmoothPathState) and isinstance(initalState, SmoothPathState)):
            raise TypeError('Start and Goal not SE2States')
        else:
            self.initialState = initalState
            self.finalState = finalState
            self.k_C0 = initalState.k
            self.k_C4 = finalState.k

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

        Eqn. 4 from paper.

        TODO: implement kDDot constraints
        """

        kTolerance = 0.01
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
        Generates a velocity trajectory which has starting velcoity equal to vInit and final velocity equal to vFinal. Time values start form tInit and increment by self.dT

        Eqn. 7 from paper.

        TODO: implement vDDot constraints
        TODO: optimize with matrix operations
        """

        vTolerance = 0.05
        v = vInit
        vTrajectory = np.array([v])

        while np.abs(v - vFinal) > vTolerance:
            if v < vFinal:
                v = v + self.dT*self.vDotMax
            else:
                v = v + self.dT*self.vDotMin
            vTrajectory = np.append(vTrajectory, np.array([v]), axis=0)

        return vTrajectory

    def makeTrajectoriesEqualLength(self, kTrajIn, vTrajIn, fromStart=False):
        """
        Takes curvature and velocity trajectories and makes them the same length. Either cutting vTraj from the start or end, or by adding values to the start or end of vTraj.
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

        return {'kTraj': kTrajIn, 'vTraj': vTraj, 'cutV': cutV, 'leftover': leftover }

    def calculateCenterLine(self):
        """
        Uses equations (17) through (21) of paper to create center line.
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

        self.omega_k = np.array(
            [S1[-1][0] - (k_C1**-1)*sin(S1[-1][2]), S1[-1][1] + (k_C1**-1)*cos(S1[-1][2])]) #center of turning of C1
        self.omega_kplus2 = np.array(
            [S4[0][0] - (k_C3**-1)*sin(S4[0][2]), S4[0][1] + (k_C3**-1)*cos(S4[0][2])]) #center of turning of C3
        d1 = np.linalg.norm(self.omega_kplus2 - self.omega_k)
        d2 = np.linalg.norm(omega_S2_tC2 - omega_S2_tS2)
        d3 = np.linalg.norm(omega_S3_tC3 - omega_S3_tS3)
        l1 = (d2**2 - d3**2 + d1**2)/(2*d1)
        if(l1 > d2):
            raise Exception
        l2 = np.sqrt(d2**2 - l1**2)
        signVal = 1 if (k_C2 < 0) else -1
        self.omega_kplus1 = self.omega_k + l1*(self.omega_kplus2 - self.omega_k)/d1 + signVal*l2*np.matmul(
            np.array([[0, -1], [1, 0]]), (self.omega_kplus2 - self.omega_k))/d1 #center of turning of C2

        return self.omega_kplus1

    def plan(self, dT):

        self.path_is_not_feasible = True
        self.dT = dT

        while self.path_is_not_feasible:

            ################ generate first connecting spiral ################
            K_S1 = self.generateCurvatureTrajectory(
                self.k_C0, self.k_C1)
            v_S1 = self.generateSpeedTrajectory(
                self.initialState.v, self.headlandSpeed)

            trajectories = self.makeTrajectoriesEqualLength(K_S1, v_S1, False)

            self.cut_v_S1 = trajectories['cutV']
            v_C1 = trajectories['leftover']
            
            xo = [self.initialState.x, self.initialState.y,self.initialState.theta]
            self.S1 = SpiralSegment(trajectories['kTraj'], trajectories['vTraj'], xo, self.dT)
            
            ################ generate last connecting spiral ################
            K_S4 = self.generateCurvatureTrajectory(self.k_C3, self.k_C4)
            v_S4 = self.generateSpeedTrajectory(self.headlandSpeed, self.finalState.v)
            trajectories = self.makeTrajectoriesEqualLength(K_S4, v_S4, True)
            
            self.cut_v_S4 = trajectories['cutV']
            v_C3 = trajectories['leftover']

            xo = [0, 0, 0]
            self.S4 = SpiralSegment(
                trajectories['kTraj'], trajectories['vTraj'], xo, self.dT)

            self.S4.placePath(self.finalState.x, self.finalState.y, self.finalState.theta, False)

            ################ generate second connecting spiral ################
            
            xo = [0, 0, 0]
            if self.cut_v_S1:
                if not(hasattr(self,'vC1_index')): 
                    #initial guess: C1 is long enough that vC1 is used completely
                    self.vC1_index = len(v_C1)
            else: #know VC1 and VC2
                self.vC1_index = len(v_C1)
            
            if self.reverse:
                
                v_S21 = self.generateSpeedTrajectory(v_C1[max(self.vC1_index-1,0)], 0) #max function is in case vC1_index = 0
                K_S21 = self.generateCurvatureTrajectory(self.k_C1, 0)

                if len(v_S21) > len(K_S21):
                    diff = len(v_S21) - len(K_S21)
                    K_S21 = np.append(K_S21, np.zeros(diff))
                elif len(v_S21) < len(K_S21):
                    diff = len(K_S21) - len(v_S21)
                    v_S21 = np.append(v_S21, np.zeros(diff))

                v_S22 = self.generateSpeedTrajectory(0, self.headlandSpeedReverse)
                K_S22 = self.generateCurvatureTrajectory(0, self.k_C2)
                trajectories = self.makeTrajectoriesEqualLength(K_S22,v_S22, False)
                self.cut_v_S2 = trajectories['cutV']
                v_C21 = trajectories['leftover']
                v_S2 = np.append(v_S21, trajectories['vTraj'])
                K_S2 = np.append(K_S21, trajectories['kTraj'])
                self.S2 = SpiralSegment(K_S2, v_S2, xo, self.dT)
            else:
                K_S2 = self.generateCurvatureTrajectory(self.k_C1, self.k_C2)
                v_S2 = v_C1[min(self.vC1_index, len(v_C1) - 1):len(v_C1)] # min function for the case where vC1_index > len(v_C1) - 1
                trajectories = self.makeTrajectoriesEqualLength(K_S2, v_S2, False)
                self.cut_v_S2 = trajectories['cutV']
                v_C21 = trajectories['leftover']
                self.S2 = SpiralSegment(trajectories['kTraj'], trajectories['vTraj'], xo, self.dT)

            ################ generate third connecting spiral ################
            
            if self.cut_v_S4:
                if not(hasattr(self,'vC4_index')):
                    #initial guess: C3 is long enough that vC3 is used completely
                    self.vC3_index = 0
            else: #know value of v_C4
                self.vC3_index = 0
            
            if self.reverse:
                v_S31 = self.generateSpeedTrajectory(self.headlandSpeedReverse, 0)
                K_S31 = self.generateCurvatureTrajectory(self.k_C2, 0)

                trajectories = self.makeTrajectoriesEqualLength(K_S31, v_S31, True)
                self.cut_v_S3 = trajectories['cutV']
                v_C22 = trajectories['leftover']

                v_S32 = self.generateSpeedTrajectory(0, v_C3[min(self.vC3_index, len(v_C3) - 1)])
                K_S32 = self.generateCurvatureTrajectory(0, self.k_C3)
                
                if len(v_S32) > len(K_S32):
                    diff = len(v_S32) - len(K_S32)
                    K_S32 = np.append(np.zeros(diff), K_S32)
                elif len(v_S32) < len(K_S32):
                    diff = len(K_S32) - len(v_S32)
                    v_S32 = np.append(np.zeros(diff),v_S32)

                v_S3 = np.append(v_S31, v_S32)
                K_S3 = np.append(K_S31, K_S32)

                self.S3 = SpiralSegment(K_S3, v_S3, xo, self.dT)
            else:

                K_S3 = self.generateCurvatureTrajectory(self.k_C2, self.k_C3)
                v_S3 = self.generateSpeedTrajectory(self.headlandSpeed, v_C3[min(self.vC3_index, len(v_C3) - 1)])
                trajectories = self.makeTrajectoriesEqualLength(K_S3, v_S3, True)
                self.cut_v_S3 = trajectories['cutV']
                v_C22 = trajectories['leftover']
                self.S3 = SpiralSegment(trajectories['kTraj'], trajectories['vTraj'], xo, self.dT)

            ################ generate center CC segment ################
            if np.abs(self.k_C2) > 0.005:  # center is an arc
                
                try:
                    self.calculateConstantArcs()
                except: #can't use this turning type
                    break 
                
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
                self.C2 = C2ArcSegment(self.k_C2, v_C21, v_C22, self.S2.poses[-1], self.S3.poses[0], self.omega_kplus1, self.dT)

            else: #center is a line

                ################ place S2 and S3 segments ################
                self.phi_C2 = self.calculateCenterLine()
                self.S2.placePath(self.S1.poses[-1][0], self.S1.poses[-1][1], self.S1.poses[-1][2])
                rotAngle = -1*self.phi_C2 - self.S2.poses[-1][2]
                self.S2.rotateAboutPoint(self.omega_k[0], self.omega_k[1], rotAngle)

                self.S3.placePath(self.S4.poses[0][0], self.S4.poses[0][1], self.S4.poses[0][2], False)
                rotAngle = self.phi_C2 - self.S3.poses[0][2]
                self.S3.rotateAboutPoint(self.omega_kplus2[0], self.omega_kplus2[1], rotAngle)

                #check if they intersect
                if self.S2.pathIntersectsWith(self.S3):
                    self.k_C2 = self.k_C2 + self.k_C4*0.1 #increases or decreases the center curvature toward k_end
                    # TODO: implemetnation for RSL LSR cases?
                    continue
                
                ################ make C2 segment ################
                self.C2 = C2LineSegment(v_C21, v_C22, self.S2.poses[-1], self.S3.poses[0], self.dT)

            ################ generate C1 and C3 ################
            self.C1 = CCSegment(self.k_C1, v_C1, self.S1.poses[-1], self.S2.poses[0], self.omega_k, self.dT)

            self.C3 = CCSegment(self.k_C3, v_C3, self.S3.poses[-1], self.S4.poses[0], self.omega_kplus2, self.dT)
            
            self.path_is_not_feasible = False

        # plotting stuff
        plt.plot([i[0] for i in self.S1.poses], [i[1] for i in self.S1.poses])
        plt.plot([i[0] for i in self.S2.poses], [i[1] for i in self.S2.poses])
        plt.plot([i[0] for i in self.S3.poses], [i[1] for i in self.S3.poses])
        plt.plot([i[0] for i in self.S4.poses], [i[1] for i in self.S4.poses])
        plt.plot([i[0] for i in self.C1.poses], [i[1] for i in self.C1.poses])
        plt.plot([i[0] for i in self.C2.poses], [i[1] for i in self.C2.poses])
        plt.plot([i[0] for i in self.C3.poses], [i[1] for i in self.C3.poses])

        plt.plot(self.omega_k[0], self.omega_k[1], 'b^')
        if hasattr(self,'omega_kplus1'):
            plt.plot(self.omega_kplus1[0], self.omega_kplus1[1], 'b^')

        plt.plot(self.omega_kplus2[0], self.omega_kplus2[1], 'b^')

        #plt.plot([(self.k_C1**-1)*cos(theta) + self.omega_k[0] for theta in np.linspace(0, 2*np.pi, 25)],
            #      [(self.k_C1**-1)*sin(theta) + self.omega_k[1] for theta in np.linspace(0, 2*np.pi, 25)], 'r--')
        #if hasattr(self,'omega_kplus1'):
            #plt.plot([(self.k_C2**-1)*cos(theta) + self.omega_kplus1[0] for theta in np.linspace(0, 2*np.pi, 25)],
            #         [(self.k_C2**-1)*sin(theta) + self.omega_kplus1[1] for theta in np.linspace(0, 2*np.pi, 25)], 'r--')
                     
        #plt.plot([(self.k_C3**-1)*cos(theta) + self.omega_kplus2[0] for theta in np.linspace(0, 2*np.pi, 25)],
                #  [(self.k_C3**-1)*sin(theta) + self.omega_kplus2[1] for theta in np.linspace(0, 2*np.pi, 25)], 'r--')
        # plt.arrow(self.S1.poses[-1][0], self.S1.poses[-1][1], 0.1*cos(self.S1.poses[-1][2]), 0.1*sin(self.S1.poses[-1][2]), length_includes_head = True, width = 0.02, head_width = 0.03, color = 'r', alpha = 0.5)
        
        # for i in range(0, len(self.S2.poses), int(len(self.S2.poses)/10)):
        #     plt.arrow(self.S2.poses[i][0], self.S2.poses[i][1], 0.1*cos(self.S2.poses[i][2]), 0.1*sin(self.S2.poses[i][2]), length_includes_head = True, width = 0.01, head_width = 0.03, color = 'r', alpha = 0.5)
        
        # for i in range(0, len(self.S4.poses), int(len(self.S2.poses)/10)):
        #     plt.arrow(self.S4.poses[i][0], self.S4.poses[i][1], 0.1*cos(self.S4.poses[i][2]), 0.1*sin(self.S4.poses[i][2]), length_includes_head = True, width = 0.01, head_width = 0.03, color = 'r', alpha = 0.5)

        plt.xlim([-1, 5])
        plt.ylim([-1, 2.5])
        plt.savefig("trajectory.png")

        return self.path


def main():

    # x pos., ypos., orientation, speed, curvature
    initialState = SmoothPathState(0, 0, 0.5*np.pi, 1, 0)
    finalState = SmoothPathState(4.2, 0, -0.5*np.pi, 1, 0)
    L_w = 1.0
    gamma_max = np.pi/4.0
    
    turningRadius = L_w/tan(gamma_max)  # =L_w/tan(gamma_max)
    dT = 0.005

    kMax = 1/turningRadius
    kMin = -kMax
    kDotMax = 5.0  # max derivative of curvature
    kDotMin = -kDotMax  # min derivative of curvature
    kDDotMax = 1.0
    kDDotMin = -kDDotMax

    kConstraints = [kMax, kMin, kDotMax, kDotMin, kDDotMax, kDDotMin]

    vMax = 1.0
    vMin = -vMax
    vDotMax = 0.9
    vDotMin = -vDotMax
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
    R1L1R = [kMin, kMax, kMin, True]
    # L1R1L = [KMax, KMin, KMax, True]

    planSmoothInst = SmoothPathPlanner()

    pathType = RSR
    planSmoothInst.setConstraints(
        kConstraints, vConstraints, headlandSpeed, headlandSpeedReverse)
    planSmoothInst.setNominalCurvatures(
        pathType[0], pathType[1], pathType[2], pathType[3])
    planSmoothInst.setStartAndGoal(initialState, finalState)

    path = planSmoothInst.plan(dT)

    # print(path)

main()
