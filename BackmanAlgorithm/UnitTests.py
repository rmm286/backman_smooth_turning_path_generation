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
    C2 = C2ArcSegment(-1, v1, v2, start, end, center, 0.005)
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

def main():
    TestC2LineSegment()
    return 0
    
main()
