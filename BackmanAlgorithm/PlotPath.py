import numpy as np
import matplotlib.pyplot as plt

def ParseRSPath(fullFilename):

    RSPathString = open(fullFilename).read().strip().split('\n')
    RSPathArray = np.zeros((len(RSPathString), 3))

    for i in range(len(RSPathString)):
        RSPathArray_i = RSPathString[i].strip().split(' ')
        for j in range(len(RSPathArray_i)):
            RSPathArray[i][j] = float(RSPathArray_i[j])

    return RSPathArray

def plotSE2Path(pathArray):

    RSPointsXArray = np.zeros((len(pathArray)))
    RSPointsYArray = np.zeros((len(pathArray)))
    RSPointsThArray = np.zeros((len(pathArray)))

    for i in range(len(pathArray)):
        RSPointsXArray[i] = pathArray[i][0]
        RSPointsYArray[i] = pathArray[i][1]
        RSPointsThArray[i] = pathArray[i][2]

    plt.plot(RSPointsXArray, RSPointsYArray)
    for i in range(0,len(pathArray),int(max(len(pathArray)/50, 5))):
        plt.arrow(RSPointsXArray[i], RSPointsYArray[i], 0.1*np.cos(RSPointsThArray[i]),0.1*np.sin(RSPointsThArray[i]), length_includes_head = True, width = 0.002, head_width = 0.003, color = 'r', alpha = 0.5)
    plt.ylabel('y (m)')
    plt.xlabel('x (m)')
    plt.show()

def main():

    RSPathFileName = "./../PlanRS/rspath.txt"

    pathArray = ParseRSPath(RSPathFileName)
  
    plotSE2Path(pathArray)

main()