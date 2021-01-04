import numpy as np
import matplotlib.pyplot as plt

# Read RSpath.txt into array

RSPathString = open("./../PlanRS/rspath.txt").read().strip().split('\n')

RSPathArray = np.zeros((len(RSPathString), 3))

for i in range(len(RSPathString)):
    RSPathArray_i = RSPathString[i].strip().split(' ')
    for j in range(len(RSPathArray_i)):
        RSPathArray[i][j] = float(RSPathArray_i[j])

RSPointsXArray = np.zeros((len(RSPathString)))
RSPointsYArray = np.zeros((len(RSPathString)))
RSPointsThArray = np.zeros((len(RSPathString)))

for i in range(len(RSPathArray)):
    RSPointsXArray[i] = RSPathArray[i][0]
    RSPointsYArray[i] = RSPathArray[i][1]
    RSPointsThArray[i] = RSPathArray[i][2]

#plot path
plot_path = True
if plot_path:
    plt.plot(RSPointsXArray, RSPointsYArray)
    for i in range(0,len(RSPathArray),max(len(RSPathArray)/50, 5)):
        plt.arrow(RSPointsXArray[i], RSPointsYArray[i], 0.1*np.cos(RSPointsThArray[i]),0.1*np.sin(RSPointsThArray[i]), length_includes_head = True, width = 0.002, head_width = 0.003, color = 'r', alpha = 0.5)
    plt.ylabel('some numbers')
    plt.show()

def check_curve(points, tolerance, curvature):
    """ 
    input(points: 3xN np.array of SE2 poses
          tolerance: float tolerance of checking)
    
    ouput(
    """

    straight = True
    right = True
    left = True

    #lastX = points[0][0]
    #lastY = points[0][1]
    lastTheta = points[0][2]

    #test if a solid curve
    for i in range(len(points)):
        #dX = points[i][0] - lastX
        #dY = points[i][1] - lastY
        dTheta = points[i][2] - lastTheta
        
        print(dTheta)

        if np.absolute(dTheta) < tolerance*2*3.141519265:
            straight = False
        
        if dTheta - tolerance > 0:
            left = False
        
        if dTheta + tolerance < 0:
            right = False
    
        #lastX = points[i][0]
        #lastY = points[i][1]
        lastTheta = points[i][2]

    if straight + right + left > 1: 
        #raise Exception("More than one curve detected, try changing tolerance.")
        print("More than one curve detected, try changing tolerance.")
        return [straight,right,left]
    
    if straight:
        return "straight"
    
    if right:
        return "right"

    if left:
        return "left"
    
    else:
        return 0

for i in range(0,len(RSPathArray),5):
    check_curve(RSPathArray[i:i+5], 0.00001)
