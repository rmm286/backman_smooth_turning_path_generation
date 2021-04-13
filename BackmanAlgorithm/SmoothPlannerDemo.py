from SmoothPlannerClass import SmoothPathPlanner, planShortest
import matplotlib.pyplot as plt
import numpy as np
from numpy import sin, cos, tan
import time
from scipy.io import savemat

def singleSourceGoal():
    dT = 0.1
    initialState = [0.0, 0.0, 0.5*np.pi, 0, 0.3]
    finalState =  [0.0, 15.0, 0.5*np.pi, 0, 0.3]
    vConstraints = [1.0, -1.0, 1.0, -1.0, 10.0, -10.0]
    kConstraints = [0.6, -0.6, 1.0, -1.0, 10.0, -10.0]
    headlandSpeed = vConstraints[0]
    headlandSpeedReverse = vConstraints[1]

    t = time.time()
    path = planShortest(kConstraints, vConstraints, headlandSpeed, headlandSpeedReverse, initialState, finalState, dT)
    elapsed = time.time() - t

    print("Shortest path has a final time of: ", path.finalTime)
    print("Time to calculate final path: ", elapsed)

    ######################### Plotting Data ####################################
    plotPath(path)
    plotControls(path, dT)
    plotCurveAndSpeed(path, dT)

    ######################### Output to Txt ######################################
    controls_file = open('./logs/control_out.txt', 'w')
    pose_file = open('./logs/pose_out.txt', 'w')
    for row in path.controls:
        np.savetxt(controls_file, row)
    for row in path.poses:
        np.savetxt(pose_file, row)
    
    pose_file.close()
    controls_file.close()

    #controls = np.loadtxt('./logs/control_out.txt').reshape(-1,2)

    print(path.poses)
    print(path.controls) 

def compareHighLowDim():
    dT = 0.005
    initialState = [0.0, 0, 0.5*np.pi, 1, 0]
    finalState =  [1.0, 1.0, -0.5*np.pi, 1, 0]
    vConstraints = [1.0, -1.0, 2.0, -2.0, 5.0, -5.0]
    kConstraints = [0.785, -0.785, 5.0, -5.0, 30.0, -30.0]
    headlandSpeed = vConstraints[0]
    headlandSpeedReverse = vConstraints[1]

    t = time.time()
    pathHighDim = planShortest(kConstraints, vConstraints, headlandSpeed, headlandSpeedReverse, initialState, finalState, dT)
    elapsedHighDim = time.time() - t
    print("High dim. path has a final time of: ", pathHighDim.finalTime)
    print("Time to calculate high dim. path: ", elapsedHighDim)
    
    t = time.time()
    pathLowDim = planShortest(kConstraints, vConstraints, headlandSpeed, headlandSpeedReverse, initialState, finalState, dT, useDDotArg= False)
    elapsedLowDim = time.time() - t
    print("Low dim. path has a final time of: ", pathLowDim.finalTime)
    print("Time to calculate low dim. path: ", elapsedLowDim)

    path = pathHighDim
    plt.figure()
    plt.arrow(path.poses[0][0], path.poses[0][1], 0.1*np.cos(path.poses[0][2]), 0.1*np.sin(path.poses[0][2]), length_includes_head = True, width = 0.01, head_width = 0.03, color = 'r', alpha = 0.5)
    plt.arrow(path.poses[-1][0], path.poses[-1][1], 0.1*np.cos(path.poses[-1][2]), 0.1*np.sin(path.poses[-1][2]), length_includes_head = True, width = 0.01, head_width = 0.03, color = 'b', alpha = 0.5)
    plt.plot([i[0] for i in path.poses], [i[1] for i in path.poses], 'k')
    path = pathLowDim
    plt.plot([i[0] for i in path.poses], [i[1] for i in path.poses], '--k')
    plt.legend(['Constraints on $\ddot{\kappa}$, $\ddot{v}$', 'Constraints only on $\dot{\kappa}$, $\dot{v}$'])
    plt.savefig('./logs/ComparePaths.png')

    fig, (ax1, ax2) = plt.subplots(2)
    fig.suptitle('Curvature Profile')
    ax1.plot([i*dT for i in range(len(path.controls.T[1]))],path.controls.T[1], 'k')
    ax1.set_ylabel('$\kappa$ $(m^{-1})$')
    ax2.plot([i*dT for i in range(len(path.controls.T[1])-1)], np.diff(path.controls.T[1])/dT, 'k')
    ax2.set_ylabel('$\dot{\kappa}$ $(m \cdot s)^{-1}$')
    kappddot = np.diff(np.diff(path.controls.T[1])/dT)/dT
    for i in range(len(kappddot)-1):
        if np.abs(kappddot[i+1] - kappddot[i]) > 1 and np.abs(kappddot[i-1] - kappddot[i]) > 1:
            kappddot[i] = kappddot[i-1]
    ax2.set_xlabel('Time (s)')
    plt.savefig('./logs/' + 'ComparePaths' + 'finalCurvature.png')

    fig, (ax1, ax2) = plt.subplots(2)
    fig.suptitle('Speed Profile')
    ax1.plot([i*dT for i in range(len(path.controls.T[0]))],path.controls.T[0], 'k')
    ax1.set_ylabel('Speed $(m/s)$')
    acc = np.diff(path.controls.T[0])/dT
    for i in range(len(acc)):
        if np.abs(acc[i] - acc[i -1]) > 1 and np.abs(acc[i] - acc[i + 1]) > 1:
            acc[i] = acc[i -1]
    ax2.plot([i*dT for i in range(len(path.controls.T[0])-1)], acc, 'k')
    ax2.set_ylabel('Accel. $(m/s^2)$')
    ax2.set_xlabel('Time (s)')
    plt.savefig('./logs/'+ 'ComparePaths' + 'finalSpeed.png')

    mdic_HighDim = {"path": pathHighDim.poses, "controls": pathHighDim.controls}
    mdic_LowDim = {"path": pathLowDim.poses, "controls": pathLowDim.controls}

    savemat('./logs/HighDim.mat', mdic_HighDim)
    savemat('./logs/LowDim.mat', mdic_LowDim)

def plotCurveAndSpeed(path, dT, filePrefixStr = ''):
    fig, (ax1, ax2) = plt.subplots(2)
    fig.suptitle('Speed and Curvature Profiles')
    ax1.plot([i*dT for i in range(len(path.controls.T[1]))],path.controls.T[1], 'k')
    ax1.set_ylabel('Curvature $(m^{-1})$')
    ax2.plot([i*dT for i in range(len(path.controls.T[0]))],path.controls.T[0], 'k')
    ax2.set_ylabel('Speed $(m/s)$')
    ax2.set_xlabel('Time (s)')
    plt.savefig('./logs/' + filePrefixStr + 'finalCurveSpeed.png')

def plotControls(path, dT, filePrefixStr = ''):
    fig, (ax1, ax2, ax3) = plt.subplots(3)
    fig.suptitle('Curvature Profile')
    ax1.plot([i*dT for i in range(len(path.controls.T[1]))],path.controls.T[1], 'k')
    ax1.set_ylabel('$\kappa$ $(m^{-1})$')
    ax2.plot([i*dT for i in range(len(path.controls.T[1])-1)], np.diff(path.controls.T[1])/dT, 'k')
    ax2.set_ylabel('$\dot{\kappa}$ $(m \cdot s)^{-1}$')
    kappddot = np.diff(np.diff(path.controls.T[1])/dT)/dT
    for i in range(len(kappddot)-1):
        if np.abs(kappddot[i+1] - kappddot[i]) > 1 and np.abs(kappddot[i-1] - kappddot[i]) > 1:
            kappddot[i] = kappddot[i-1]
    ax3.plot([i*dT for i in range(len(path.controls.T[1])-2)], kappddot, 'k')
    ax3.set_xlabel('Time (s)')
    ax3.set_ylabel('$\ddot{\kappa}$ $(m \cdot s^2)^{-1}$')
    plt.savefig('./logs/' + filePrefixStr + 'finalCurvature.png')

    fig, (ax1, ax2, ax3) = plt.subplots(3)
    fig.suptitle('Speed Profile')
    ax1.plot([i*dT for i in range(len(path.controls.T[0]))],path.controls.T[0], 'k')
    ax1.set_ylabel('Speed $(m/s)$')
    ax2.plot([i*dT for i in range(len(path.controls.T[0])-1)], np.diff(path.controls.T[0])/dT, 'k')
    ax2.set_ylabel('Accel. $(m/s^2)$')
    vddot = np.diff(np.diff(path.controls.T[0])/dT)/dT
    for i in range(len(vddot)-1):
        if np.abs(vddot[i+1] - vddot[i]) > 1 and np.abs(vddot[i-1] - vddot[i]) > 1:
            vddot[i] = vddot[i-1]
    ax3.plot([i*dT for i in range(len(path.controls.T[0])-2)], vddot, 'k')
    ax3.set_xlabel('Time (s)')
    ax3.set_ylabel('Jerk $(m/s^3)$')
    plt.savefig('./logs/'+ filePrefixStr + 'finalSpeed.png')    

def plotPath(path, filePrefixStr = ''):
    plt.figure()
    plt.clf()
    plt.title('Final Path')
    plt.arrow(path.poses[0][0], path.poses[0][1], 0.1*np.cos(path.poses[0][2]), 0.1*np.sin(path.poses[0][2]), length_includes_head = True, width = 0.01, head_width = 0.03, color = 'r', alpha = 0.5)
    plt.arrow(path.poses[-1][0], path.poses[-1][1], 0.1*np.cos(path.poses[-1][2]), 0.1*np.sin(path.poses[-1][2]), length_includes_head = True, width = 0.01, head_width = 0.03, color = 'b', alpha = 0.5)
    plt.plot([i[0] for i in path.poses], [i[1] for i in path.poses], 'k')
    plt.savefig('./logs/' + filePrefixStr + 'finalPath.png')

def main():
    singleSourceGoal()
    #compareHighLowDim()



main()
