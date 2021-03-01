#!/usr/bin/env python3
from SmoothPlannerClass import SmoothPathPlanner, planShortest
import matplotlib.pyplot as plt
import numpy as np
import time

def testSingleSourceGoal():
    dT = 0.0005
    initialState = [0.0, 0, 0.5*np.pi, 1, 0]
    nIteration = 6
    finalStates = [[2*(np.cos(i)), 2*(np.sin(i)), i, 1, 0] for i in np.linspace(0,(2-1/nIteration)*np.pi,nIteration-1)]
    vConstraints = [1.0, -1.0, 2.0, -2.0, 5.0, -5.0]
    kConstraints = [0.785, -0.785, 30.0, -30.0, 30.0, -30.0]
    headlandSpeed = vConstraints[0]
    headlandSpeedReverse = vConstraints[1]
    i = 0
    for finalState in finalStates:
        i = i + 1
        if i == 5:
            path = planShortest(kConstraints, vConstraints, headlandSpeed, headlandSpeedReverse, initialState, finalState, dT, 2)
    
    plt.figure()
    plt.clf()
    plt.title("Fwd Motion Allowed")
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")
    plt.plot([i[0] for i in path.poses], [i[1] for i in path.poses])
    plt.savefig("./logs/singleSourceGoal.png")

def goalStateIterationCircle():
    #circle
    dT = 0.001
    nIteration = 6
    initialState = [0.0, 0, 0.5*np.pi, 1, 0]
    finalStates = [[2*(np.cos(i)), 2*(np.sin(i)), i, 1, 0] for i in np.linspace(0,(2-1/nIteration)*np.pi,nIteration-1)]
    vConstraints = [1.0, -1.0, 2.0, -2.0, 5.0, -5.0]
    kConstraints = [0.785, -0.785, 30.0, -30.0, 30.0, -30.0]
    plt.figure()
    plt.clf()
    plt.title("Goal State Iteration")
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")

    for finalState in finalStates:
    
        planner = SmoothPathPlanner(dT)
        planner.setConstraints(kConstraints, vConstraints, vConstraints[0], vConstraints[1])
        planner.setStartAndGoal(initialState, finalState)
        path = planner.planShortest()
        plt.plot([i[0] for i in path.poses], [i[1] for i in path.poses])

    plt.savefig("./logs/goalStateIterationCircle.png")

def goalStateIterationLine():
    #line
    dT = 0.001
    nIteration = 6
    initialState = [0.0, -1.0, 0.5*np.pi, 1, 0]
    finalStates = [[1+0.5*i, 0.0, -0.5*np.pi, 1, 0] for i in range(nIteration)]
    vConstraints = [1.0, -1.0, 2.0, -2.0, 5.0, -5.0]
    kConstraints = [0.785, -0.785, 30.0, -30.0, 30.0, -30.0]
    plt.figure()
    plt.clf()
    plt.title("Goal State Iteration")
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")

    for finalState in finalStates:
    
        planner = SmoothPathPlanner(dT)
        planner.setConstraints(kConstraints, vConstraints, vConstraints[0], vConstraints[1])
        planner.setStartAndGoal(initialState, finalState)
        path = planner.planShortest()
        plt.plot([i[0] for i in path.poses], [i[1] for i in path.poses])

    plt.savefig("./logs/goalStateIterationLine.png")

def goalStateIterationOrient():
    #final orientation

    dT = 0.001
    nIteration = 12
    initialState = [0.0, -1.0, 0.5*np.pi, 1, 0]
    finalStates = [[1, 0.0, i, 1, 0] for i in np.linspace(0,(2-1/nIteration)*np.pi,nIteration-1)]
    vConstraints = [1.0, -1.0, 2.0, -2.0, 5.0, -5.0]
    kConstraints = [0.785, -0.785, 30.0, -30.0, 30.0, -30.0]
    plt.figure()
    plt.clf()
    plt.title("Goal State Iteration")
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")

    for finalState in finalStates:
    
        planner = SmoothPathPlanner(dT)
        planner.setConstraints(kConstraints, vConstraints, vConstraints[0], vConstraints[1])
        planner.setStartAndGoal(initialState, finalState)
        path = planner.planShortest()
        plt.plot([i[0] for i in path.poses], [i[1] for i in path.poses])

    plt.savefig("./logs/goalStateIterationOrient.png")


def goalStateIterationCurvature():
    #final orientation

    dT = 0.001
    initialState = [0.0, -1.0, 0.5*np.pi, 1, 0]
    finalStates = [[1, 0.0, -0.5*np.pi, 1, 0] for i in np.linspace(-0.785,0.785,5)]
    vConstraints = [1.0, -1.0, 2.0, -2.0, 5.0, -5.0]
    kConstraints = [0.785, -0.785, 30.0, -30.0, 30.0, -30.0]
    plt.figure()
    plt.clf()
    plt.title("Goal State Iteration")
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")

    for finalState in finalStates:
    
        planner = SmoothPathPlanner(dT)
        planner.setConstraints(kConstraints, vConstraints, vConstraints[0], vConstraints[1])
        planner.setStartAndGoal(initialState, finalState)
        path = planner.planShortest()
        plt.plot([i[0] for i in path.poses], [i[1] for i in path.poses])

    plt.savefig("./logs/goalStateIterationCurvature.png")

def main():

    #goalStateIterationLine()
    #goalStateIterationCircle()
    #goalStateIterationOrient()
    #goalStateIterationCurvature()
    testSingleSourceGoal()

main()