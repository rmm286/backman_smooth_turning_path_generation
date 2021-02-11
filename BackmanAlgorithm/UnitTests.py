#from SmoothPlanner import SmoothPathPlanner
from PathSegment import SpiralSegment, LineSegment, CCSegment, C2ArcSegment, C2LineSegment
import numpy as np
import matplotlib.pyplot as plt

def TestC2ArcSegment():
    v1 = np.array([1, 1.1,1.2, 1.3, 1.4, 1.5])
    v2 = np.array([1.5, 1.4, 1.3, 1.2, 1.1, 1])
    start = np.array([np.cos((1.3)*np.pi), np.sin((1.3)*np.pi), (2.3)*np.pi-(1/2)*np.pi])
    end = np.array([np.cos((1.0/4.0)*np.pi), np.sin((1.0/4.0)*np.pi), (1.0/4.0)*np.pi-(1/2)*np.pi])
    center = np.array([0,0])
    C2 = C2ArcSegment(0.05, -1, v1, v2, start, end, center, 0.005)
    plt.plot([i[0] for i in C2.poses], [i[1] for i in C2.poses])
    plt.show()

def TestC2LineSegment():
    v1 = np.array([1, 1.1,1.2, 1.3, 1.4, 1.5])
    v2 = np.array([1.5, 1.4, 1.3, 1.2, 1.1, 1])
    start = np.array([np.cos((1.3)*np.pi), np.sin((1.3)*np.pi), (2.3)*np.pi-(1/2)*np.pi])
    end = np.array([np.cos((1.0/4.0)*np.pi), np.sin((1.0/4.0)*np.pi), (1.0/4.0)*np.pi-(1/2)*np.pi])
    C2 = C2LineSegment(v1, v2, start, end, 0.005)
    plt.plot([i[0] for i in C2.poses], [i[1] for i in C2.poses])
    plt.show()

def TestGenerateCurvatureTraj():
    
    K1 = 1
    K2 = -1

    k = K1
    kFinal = K2
    k0 = K1
    kDot0 = 0
    kDotf = 0

    kDotMax = 0.5
    kDotMin = -0.5
    kDDotMax = 1.2
    kDDotMin = -0.8
    kDot = kDot0
    dT = 0.01

    kTrajectory = np.array([k])
    kDotTrajectory = np.array([kDot])

    if k0 < kFinal: #increase k to kf
        riseTime = (kDotMax - k0)/(kDDotMax)
        fallTime = (k0 - kDotMax)/(kDDotMin)
        kRT = k0 + (kDDotMax/2.0)*riseTime**2
        kFT = kRT + (kDDotMax*riseTime)*fallTime + (kDDotMin/2.0)*fallTime**2

        if kFT < kFinal:
            diff = kFinal - kFT
            tTop = diff/kDotMax
            k0toRT = np.array([k0 + (kDDotMax/2.0)*t**2 for t in np.linspace(0,riseTime,int(riseTime/dT))])
            kRTtoDiff = np.array([k0toRT[-1] + kDDotMax*riseTime*t for t in np.linspace(dT, tTop, int(tTop/dT))])
            kDifftoFT = np.array([kRTtoDiff[-1] + kDDotMax*riseTime*t + 0.5*kDDotMin*t**2 for t in np.linspace(dT,fallTime,int(fallTime/dT))])
            kTrajectory = np.append(k0toRT,np.append(kRTtoDiff,kDifftoFT,axis=0),axis=0)
        else:
            t1 = np.sqrt((kFinal - k0)/((kDDotMax/2.0)*(1-kDDotMax/kDDotMin)))
            t2 = -1*(kDDotMax/kDDotMin)*t1
            k0tot1 = np.array([k0 + kDDotMax/2.0*t**2 for t in np.linspace(0,t1,int(t1/dT))])
            kt1tot2 = np.array([k0tot1[-1] + (kDDotMax*t1)*t + (kDDotMin/2.0)*t**2 for t in np.linspace(dT,t2,int(t2/dT))])
            kTrajectory = np.append(k0tot1,kt1tot2,axis=0)
    if k0 > kFinal: #decrease k to kf
        fallTime = (kDotMin - k0)/(kDDotMin)
        riseTime = (k0 - kDotMin)/(kDDotMax)
        kFT = k0 + (kDDotMin/2.0)*fallTime**2
        kRT = kFT + (kDDotMin*fallTime)*riseTime + (kDDotMax/2.0)*riseTime**2

        if kRT > kFinal:
            diff = kFinal - kRT
            tBottom = diff/kDotMin
            k0toFT = np.array([k0 + (kDDotMin/2.0)*t**2 for t in np.linspace(0,fallTime,int(fallTime/dT))])
            kFTtoDiff = np.array([k0toFT[-1] + kDDotMin*fallTime*t for t in np.linspace(dT, tBottom, int(tBottom/dT))])
            kDifftoRT = np.array([kFTtoDiff[-1] + kDDotMin*fallTime*t + 0.5*kDDotMax*t**2 for t in np.linspace(dT,riseTime,int(riseTime/dT))])
            kTrajectory = np.append(k0toFT,np.append(kFTtoDiff,kDifftoRT,axis=0),axis=0)
        else:
            t1 = np.sqrt((kFinal - k0)/((kDDotMin/2.0)*(1-kDDotMin/kDDotMax)))
            t2 = -1*(kDDotMin/kDDotMax)*t1
            k0tot1 = np.array([k0 + kDDotMin/2.0*t**2 for t in np.linspace(0,t1,int(t1/dT))])
            kt1tot2 = np.array([k0tot1[-1] + (kDDotMin*t1)*t + (kDDotMax/2.0)*t**2 for t in np.linspace(dT,t2,int(t2/dT))])
            kTrajectory = np.append(k0tot1,kt1tot2,axis=0)

    print(len(kTrajectory))        

    plt.figure(0)
    plt.clf
    plt.plot([dT*i for i in range(len(kTrajectory))],kTrajectory)
    plt.show()       

def main():
    TestGenerateCurvatureTraj()
    return 0
    
main()
