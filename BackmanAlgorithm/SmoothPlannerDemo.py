from SmoothPlannerClass import SmoothPathPlanner
import matplotlib.pyplot as plt
import numpy as np
from numpy import sin, cos, tan
import time

def main():

    # x pos., ypos., orientation, speed, curvature
    initialState = [0.0, 0.0, 0.5*np.pi, 1, 0]
    finalState = [1.0, 0.0, -0.5*np.pi, 1, 0]
    L_w = 1.0
    gamma_max = np.pi/4.0
    turningRadius = L_w/tan(gamma_max)
    dT = 0.005

    kMax = 0.785
    kMin = -0.785
    kDotMax = 2.0  # max derivative of curvature
    kDotMin = -kDotMax  # min derivative of curvature
    kDDotMax = 3.0
    kDDotMin = -3.0
    vMax = 1.0
    vMin = -vMax
    vDotMax = 2.0
    vDotMin = -vDotMax
    vDDotMax = 3.0
    vDDotMin = -3.0
    headlandSpeed = vMax
    headlandSpeedReverse = vMin

    vConstraints = [vMax, vMin, vDotMax, vDotMin, vDDotMax, vDDotMin]
    kConstraints = [kMax, kMin, kDotMax, kDotMin, kDDotMax, kDDotMin]

    planSmoothInst = SmoothPathPlanner(dT)
    planSmoothInst.setConstraints(kConstraints, vConstraints, headlandSpeed, headlandSpeedReverse)
    planSmoothInst.setStartAndGoal(initialState, finalState)

    t = time.time()
    shortestPath = planSmoothInst.planShortest()
    elapsed = time.time() - t

    # t = time.time()
    # planSmoothInst.setNominalCurvatures(kMin, kMax, kMin, True)
    # shortestPath = planSmoothInst.plan()
    # elapsed = time.time() - t

    planSmoothInst.plotPaths()
    planSmoothInst.plotControls()
    
    plt.figure(3)
    plt.clf()
    plt.title("Final Path")
    text = "Final Time: " + str(shortestPath.finalTime)
    plt.annotate(text, [0, 0])
    plt.plot([i[0] for i in shortestPath.poses], [i[1] for i in shortestPath.poses])
    plt.savefig("finalPath.png")
    print("Shortest path has a final time of: ", shortestPath.finalTime)
    print("Time to calculate final path: ", elapsed)

main()