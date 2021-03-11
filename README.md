# Kidnapped Vehicle Project
[![Udacity - Self-Driving Car NanoDegree](https://s3.amazonaws.com/udacity-sdc/github/shield-carnd.svg)](http://www.udacity.com/drive)

In this project of the Udacity [Self-Driving Car NanoDegree](https://www.udacity.com/course/self-driving-car-engineer-nanodegree--nd013) program, I implement a 2-dimensional particle filter in order to predict and track the position of a moving object.

A [map of landmarks](data/map_data.txt) and a noisy GPS estimate of the object's initial position is given. Using the Udacity [Term 2 Simulator](https://github.com/udacity/self-driving-car-sim/releases) sensor and control data is fed to the filter logic continuously. The particle with highest weight (best fit) is provided back to the simulator for evaluation. 

## Resources
* [Self-Driving Car NanoDegree](https://www.udacity.com/course/self-driving-car-engineer-nanodegree--nd013) course description at Udacity
* [Kidnapped Vehicle Project Starter Code](https://github.com/udacity/CarND-Kidnapped-Vehicle-Project) on Github

## Prerequisites
You can find detailed instructions for setting up your system on the starter code [project page](https://github.com/udacity/CarND-Kidnapped-Vehicle-Project). It involves installing a Unity-based [simulator tool](https://github.com/udacity/self-driving-car-sim/releases) and a C++ webserver library [uWebSocketIO](https://github.com/uWebSockets/uWebSockets).

## Building and Running the Code
1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make` 
   * On windows, you may need to run: `cmake .. -G "Unix Makefiles" && make`
4. Execute: `./particle_filter`

## Structure
The original starter code was modified in the following files to complete the project assignment:
* [particle_filter.h](src/particle_filter.h) / [particle_filter.cpp](src/particle_filter.cpp): Particle filter implementation (initialization, predication, weights, resampling)
* [helper_functions.h](src/helper_functions.h): Helper functions (Euclidean distance, multivariate Gaussian, helper structs, etc.)

## License
The contents of this repository are covered under the [MIT License](LICENSE).
