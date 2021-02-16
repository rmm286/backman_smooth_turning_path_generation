from SmoothPlannerClass import SmoothPathPlanner
import matplotlib.pyplot as plt
import numpy as np
from numpy import sin, cos, tan
import time

def main():

    # x pos., ypos., orientation, speed, curvature
    initialState = [0.0, 0.0, 0.5*np.pi, 0.8, 0]
    finalState = [20, 0.0, -0.5*np.pi, 0.8, 0]
    L_w = 1.0
    gamma_max = np.pi/4.0
    turningRadius = L_w/tan(gamma_max)
    dT = 0.005
    kMax = 1/turningRadius
    kMin = -kMax
    kDotMax = 5.0  # max derivative of curvature
    kDotMin = -kDotMax  # min derivative of curvature
    kDDotMax = 5.0
    kDDotMin = -5.0
    vMax = 1.0
    vMin = -vMax
    vDotMax = 1.0
    vDotMin = -vDotMax
    vDDotMax = 5.0
    vDDotMin = -5.0
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