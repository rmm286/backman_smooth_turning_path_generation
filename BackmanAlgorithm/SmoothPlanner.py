import numpy as np
from numpy import sin, cos, tan
from ODESolver import RungeKutta4
from StateSpaces import SmoothPathState
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from PathSegment import SpiralSegment, CCSegment, LineSegment

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

        self.vMax = kConstraints[0]
        self.vMin = kConstraints[1]
        self.vDotMax = kConstraints[2]
        self.vDotMin = kConstraints[3]
        self.vDDotMax = kConstraints[4]
        self.vDDotMin = kConstraints[5]

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

    def generateCurvatureTrajectory(self, kInit, kFinal, tInit):
        """
        Generates a curvature trajectory which has starting curvature equal to kInit and final curvature equal to kFinal. Time values start form tInit and increment by self.dT

        Eqn. 4 from paper.

        TODO: implement kDDot constraints
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
        Generates a velocity trajectory which has starting velcoity equal to vInit and final velocity equal to vFinal. Time values start form tInit and increment by self.dT

        Eqn. 7 from paper.

        TODO: implement vDDot constraints
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

    def makeTrajectoriesEqualLength(self, kTraj, vTraj, fromStart=False):
        """
        Takes curvature and velcity trajectories and makes them the same length. Either cutting vTraj from the start or end, or by adding values to the start or end of vTraj.
        """

        if (len(kTraj) < len(vTraj)) and not fromStart:  # cut from end of vTraj
            vTraj = vTraj[0:len(kTraj)]
        elif (len(kTraj) < len(vTraj)) and fromStart:  # cut from start of vTraj
            vTraj = vTraj[len(vTraj) - len(kTraj):len(vTraj)]
        elif (len(kTraj) > len(vTraj)) and not fromStart: # add to end of vTraj
            vTraj = np.append(vTraj, np.array(
                [[vTraj[-1][0], kTraj[i + len(vTraj)][1]]  for i in range((len(kTraj) - len(vTraj)))]), axis=0)
        elif (len(kTraj) > len(vTraj)) and fromStart: # add to start of vTraj
            extension = np.array([[vTraj[0][0], kTraj[i][1]]for i in range(len(kTraj) - len(vTraj))])
            vTrajAdjustedTimes = np.array([[vTraj[i- len(vTraj)][0], kTraj[i + len(vTraj)-1][1]] for i in range(len(vTraj))])
            vTraj = np.append(extension,vTrajAdjustedTimes, axis = 0)

        return {'kTraj': kTraj, 'vTraj': vTraj}

    def makeTrajectoriesEqualLengthAndMatchZeros(self, kTraj, vTraj):
        """
        From section 2.2 curvature and velocity must be zero at point where velocity changes
        From predefined turning types the curvature trajectory will always have a zero point when manuever involved reversing direction. 

        This function matches the zero positions and then calls makeTrajectoriesEqualLength on result.
        """

        # find zero point in kTraj
        kValArray = [np.abs(i[0]) for i in kTraj]
        vValArray = [np.abs(i[0]) for i in vTraj]
        kZeroIndex = np.argmin(kValArray)
        vZeroIndex = np.argmin(vValArray)

        kZeroTime = kTraj[kZeroIndex][1]

        if kZeroIndex < vZeroIndex: #extend kTraj
            diff = vZeroIndex - kZeroIndex
            kTrajExtension = np.array(
                [[0, kTraj[kZeroIndex][1] + (i+1)*self.dT] for i in range(diff)])
            kTraj = np.insert(kTraj, kZeroIndex, kTrajExtension, axis=0)

            return self.makeTrajectoriesEqualLength(kTraj, vTraj)
        else:
            diff = kZeroIndex - vZeroIndex
            vTrajExtension = np.array(
                [[0, vTraj[kZeroIndex][1] + (i+1)*self.dT] for i in range(diff)])
            vTraj = np.insert(vTraj, vZeroIndex, vTrajExtension, axis=0)

            return self.makeTrajectoriesEqualLength(kTraj, vTraj)

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
                self.k_C0, self.k_C1, 0)
            v_S1 = self.generateSpeedTrajectory(
                self.initialState.v, self.headlandSpeed, 0)

            trajectories = self.makeTrajectoriesEqualLength(K_S1, v_S1, False)

            if len(v_S1) > len(trajectories['vTraj']):
                cut_v_S1 = True
                cut_v_S1_index = len(trajectories['vTraj'])-1
            else:
                cut_v_S1 = False
                v_C1 = trajectories['vTraj'][-1][0]
            xo = [self.initialState.x, self.initialState.y,
                  self.initialState.theta]
            self.S1 = SpiralSegment(trajectories['kTraj'], trajectories['vTraj'], xo)
            

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
                v_C3 = trajectories['vTraj'][0][0]

            xo = [0, 0, 0]
            self.S4 = SpiralSegment(
                trajectories['kTraj'], trajectories['vTraj'], xo)

            self.S4.placePath(self.finalState.x, self.finalState.y, self.finalState.theta, False)

            ################ generate second conneting spiral ################
            K_S2 = self.generateCurvatureTrajectory(self.k_C1, self.k_C2, 0)
            xo = [0, 0, 0]
            #TODO: special circumstance where K trajectory is not long enough to allow for changes to velocity, in this case extend K trajectory by adding zeros to end.
            if self.reverse:
                v_S2 = self.generateSpeedTrajectory(self.headlandSpeed, 0, 0)
                v_S2 = np.append(v_S2, self.generateSpeedTrajectory(
                    0, self.headlandSpeedReverse, v_S2[-1][1]), axis=0)
                trajectories = self.makeTrajectoriesEqualLengthAndMatchZeros(
                    K_S2, v_S2)
                self.S2 = SpiralSegment(
                    trajectories['kTraj'], trajectories['vTraj'], xo)
            else:
                if cut_v_S1:
                    v_S2 = v_S1[cut_v_S1_index+1:len(v_S2)]
                    extend_v_S2 = np.array([[v_S1[-1][0], K_S2[i + len(v_S2)][1]]
                                            for i in range(len(K_S2) - len(v_S2))])
                    v_S2 = np.append(v_S2, extend_v_S2, axis=0)
                else:
                    v_S2 = np.array([[v_C1, K_S2[i][1]] for i in range(len(K_S2))])

                self.S2 = SpiralSegment(K_S2, v_S2, xo)

            ################ generate thrid connecting spiral ################
            K_S3 = self.generateCurvatureTrajectory(self.k_C2, self.k_C3, 0)
            #TODO: special circumstance where K trajectory is not long enough to allow for changes to velocity, in this case extend K trajectory by adding zeros to end.
            if self.reverse:
                
                v_S3 = self.generateSpeedTrajectory(self.headlandSpeedReverse, 0, 0)
                v_S3 = np.append(v_S3, self.generateSpeedTrajectory(
                    0, self.headlandSpeed, v_S3[-1][1]), axis=0)
                trajectories = self.makeTrajectoriesEqualLengthAndMatchZeros(
                    K_S3, v_S3)
                self.S3 = SpiralSegment(
                    trajectories['kTraj'], trajectories['vTraj'], xo)
            else:
                if cut_v_S4:
                    v_S3 = v_S4[0:cut_v_S4_index]
                    extend_v_S3 = np.array([[v_S4[-1][0], K_S3[i+len(v_S3)][1]]
                                            for i in range(len(K_S3) - len(v_S3))])
                    v_S3 = np.append(v_S3, extend_v_S3, axis=0)
                else:
                    v_S3 = np.array([[v_C3, K_S3[i][1]] for i in range(len(K_S3))])

                self.S3 = SpiralSegment(K_S3, v_S3, xo)

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

            ################ generate C1 and C3 ################
            v_C1 = np.array([[self.headlandSpeed, self.dT * i] for i in range(100)])
            self.C1 = CCSegment(self.k_C1, v_C1, self.S1.poses[-1], self.S2.poses[0], self.omega_k, self.dT)

            v_C3 = np.array([[self.headlandSpeed, self.dT * i] for i in range(100)])
            self.C3 = CCSegment(self.k_C3, v_C3, self.S3.poses[-1], self.S4.poses[0], self.omega_kplus2, self.dT)


            self.path_is_not_feasible = False

        # plotting stuff
        plt.plot([i[0] for i in self.S1.poses], [i[1] for i in self.S1.poses])
        plt.plot([i[0] for i in self.S2.poses], [i[1] for i in self.S2.poses])
        plt.plot([i[0] for i in self.S3.poses], [i[1] for i in self.S3.poses])
        plt.plot([i[0] for i in self.S4.poses], [i[1] for i in self.S4.poses])
        plt.plot([i[0] for i in self.C1.poses], [i[1] for i in self.C1.poses])
        plt.plot([i[0] for i in self.C3.poses], [i[1] for i in self.C3.poses])

        plt.plot(self.omega_k[0], self.omega_k[1], 'b^')
        if hasattr(self,'omega_kplus1'):
            plt.plot(self.omega_kplus1[0], self.omega_kplus1[1], 'b^')

        plt.plot(self.omega_kplus2[0], self.omega_kplus2[1], 'b^')

        plt.plot([(self.k_C1**-1)*cos(theta) + self.omega_k[0] for theta in np.linspace(0, 2*np.pi, 25)],
                  [(self.k_C1**-1)*sin(theta) + self.omega_k[1] for theta in np.linspace(0, 2*np.pi, 25)], 'r--')
        if hasattr(self,'omega_kplus1'):
            plt.plot([(self.k_C2**-1)*cos(theta) + self.omega_kplus1[0] for theta in np.linspace(0, 2*np.pi, 25)],
                     [(self.k_C2**-1)*sin(theta) + self.omega_kplus1[1] for theta in np.linspace(0, 2*np.pi, 25)], 'r--')
                     
        plt.plot([(self.k_C3**-1)*cos(theta) + self.omega_kplus2[0] for theta in np.linspace(0, 2*np.pi, 25)],
                  [(self.k_C3**-1)*sin(theta) + self.omega_kplus2[1] for theta in np.linspace(0, 2*np.pi, 25)], 'r--')
        # plt.arrow(self.S1.poses[-1][0], self.S1.poses[-1][1], 0.1*cos(self.S1.poses[-1][2]), 0.1*sin(self.S1.poses[-1][2]), length_includes_head = True, width = 0.02, head_width = 0.03, color = 'r', alpha = 0.5)
        
        # for i in range(0, len(self.S2.poses), int(len(self.S2.poses)/10)):
        #     plt.arrow(self.S2.poses[i][0], self.S2.poses[i][1], 0.1*cos(self.S2.poses[i][2]), 0.1*sin(self.S2.poses[i][2]), length_includes_head = True, width = 0.01, head_width = 0.03, color = 'r', alpha = 0.5)
        
        # for i in range(0, len(self.S4.poses), int(len(self.S2.poses)/10)):
        #     plt.arrow(self.S4.poses[i][0], self.S4.poses[i][1], 0.1*cos(self.S4.poses[i][2]), 0.1*sin(self.S4.poses[i][2]), length_includes_head = True, width = 0.01, head_width = 0.03, color = 'r', alpha = 0.5)

        plt.xlim([-1, 2.5])
        plt.ylim([-1, 2.5])
        plt.savefig("trajectory.png")

        return self.path


def main():

    # x pos., ypos., orientation, speed, curvature
    initialState = SmoothPathState(0, 0, 0.5*np.pi, 1, 0)
    finalState = SmoothPathState(1, 0, -0.5*np.pi, 1, 0)
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
    R1L1R = [kMin, kMax, kMin, True]
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
