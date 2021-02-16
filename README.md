# Continuous and Differentiable Path Planner

* Author: Rayne Milner
* Maintainer: Rayne Milner
* Contact: rmmilner@ucdavis.edu
* Org: University of California at Davis

## Overview

***

This Path planning algorithm is a procedure for generating a path from start to goal state that respects constraints on the 0th, 1st and 2nd derivative of curvature and velocity. The algorithm is taken from Juha Backman's 2015 paper "Smooth turning path generation for agricultural vehicles in headlands". 

## Dependencies

***

* Scipy: <https://www.scipy.org/install.html>

* Python3.5+: <https://www.python.org/downloads/release/python-350/>

## Installation

***
Files are ready for use on download.

Usage only requires one import command `from SmoothPlannerClass import SmoothPathPlanner`.

## Usage

***

The planner can be used as a generic planner to find the shortest path between two states with elements: [x, y, theta, velocity, curvature].

The typical sequence of commands to call the planner method is as follows:

1. Instantiate Planner Class: `inst = SmoothPathPlanner(timeStep)`
2. Set Constraints and desired speeds: `inst.setConstraints(kConstraints, vConstraints, headlandSpeed, headlandSpeedReverse)`
3. Set Start and Goal States: `inst.setStartAndGoal(initialState, finalState)`
4. Call Shortest planner method `shortestPath = planSmoothInst.planShortest()`

The returned path then has member data which contains the sequence of states from start to goal as well as controls.

A few other functionalities are provided by the class:

* If the user desired to plot a path with a specific turning type instead of the shortest path of any type, the user can choose to specify the path parameters with `setNominalCurvatures` and then instead of calling `planShortest()` call `plan()` 

* The path can be plotted using `plotPath()`

* Control inputs can be plotted with `plotControls()` once the final path has been returned (plan or planShortest has been called).

## External Links

***

* <https://www.sciencedirect.com/science/article/abs/pii/S1537511015001397>

## License

***

MIT License

Copyright (c) [2021] [Rayne Michael Milner]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.