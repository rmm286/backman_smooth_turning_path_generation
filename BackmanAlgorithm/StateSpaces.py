class SmoothPathState:
    """
    Data structure for state of robot in smooth path planner.

    x position, y position, orientation (theta), speed, and curvature
    """

    def __init__(self, x, y, theta, v, k):

        self.x = x
        self.y = y
        self.theta = theta
        self.v = v
        self.k = k

class SE2State:
    """
    Data structure for state of robot in SE2.

    x position, y position, orientation (theta)
    """

    def __init__(self, x, y, theta):

        self.x = x
        self.y = y
        self.theta = theta