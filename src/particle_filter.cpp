/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::numeric_limits;
using std::normal_distribution;
using std::default_random_engine;
using std::discrete_distribution;

constexpr double epsilon = 1E-7;

static default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[])
{
    /**
     * TODO: Set the number of particles. Initialize all particles to
     *   first position (based on estimates of x, y, theta and their uncertainties
     *   from GPS) and all weights to 1.
     * TODO: Add random Gaussian noise to each particle.
     * NOTE: Consult particle_filter.h for more information about this method
     *   (and others in this file).
     */

    constexpr size_t num_particles = 10000;

    normal_distribution<double> dist_x(x, std[0]);
    normal_distribution<double> dist_y(y, std[1]);
    normal_distribution<double> dist_theta(theta, std[2]);

    particles.resize(num_particles);
    for (size_t i = 0; i < num_particles; i++) {
        particles[i].id = i;
        particles[i].x = dist_x(gen);
        particles[i].y = dist_y(gen);
        particles[i].theta = dist_theta(gen);
        particles[i].weight = 1.0;
    }
}

void ParticleFilter::prediction(double delta_t, double std_pos[],
    double velocity, double yaw_rate)
{
    /**
     * TODO: Add measurements to each particle and add random Gaussian noise.
     * NOTE: When adding noise you may find std::normal_distribution
     *   and std::default_random_engine useful.
     *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
     *  http://www.cplusplus.com/reference/random/default_random_engine/
     */

     normal_distribution<double> dist_x(0.0, std_pos[0]);
     normal_distribution<double> dist_y(0.0, std_pos[1]);
     normal_distribution<double> dist_theta(0.0, std_pos[2]);

     for (Particle& p : particles) {

         if (yaw_rate < epsilon) {
            p.x += velocity * cos(p.theta) * delta_t;
	        p.y += velocity * sin(p.theta) * delta_t;
         }
         else {
             p.x += velocity / yaw_rate * (sin( p.theta + yaw_rate*delta_t) - sin(p.theta));
             p.y += velocity / yaw_rate * (cos(p.theta) - cos( p.theta + yaw_rate*delta_t));
             p.theta += yaw_rate * delta_t;
         }

         p.x += dist_x(gen);
         p.y += dist_y(gen);
         p.theta += dist_theta(gen);
     }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted,
    vector<LandmarkObs>& observations)
{
    /**
     * TODO: Find the predicted measurement that is closest to each
     *   observed measurement and assign the observed measurement to this
     *   particular landmark.
     * NOTE: this method will NOT be called by the grading code. But you will
     *   probably find it useful to implement this method and use it as a helper
     *   during the updateWeights phase.
     */

    for (LandmarkObs& obs: observations) {
        double min_distance = numeric_limits<double>::max();

        for (const LandmarkObs& pred: predicted) {
            const double distance = dist(obs.x, obs.y, pred.x, pred.y);
            if (distance < min_distance) {
                min_distance = distance;
                obs.id = pred.id;
            }
        }
    }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
    const vector<LandmarkObs>& observations,
    const Map& map_landmarks)
{
    /**
     * TODO: Update the weights of each particle using a mult-variate Gaussian
     *   distribution. You can read more about this distribution here:
     *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
     * NOTE: The observations are given in the VEHICLE'S coordinate system.
     *   Your particles are located according to the MAP'S coordinate system.
     *   You will need to transform between the two systems. Keep in mind that
     *   this transformation requires both rotation AND translation (but no scaling).
     *   The following is a good resource for the theory:
     *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
     *   and the following is a good resource for the actual equation to implement
     *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
     */

    for (auto& p: particles) {

        // identify landmarks which are in sensor range of this particle
        vector<LandmarkObs> landmarks_in_range;
        for (const auto& lm: map_landmarks.landmark_list) {
            const double distance = dist(p.x, p.y, lm.x_f, lm.y_f);
            if (distance <= sensor_range)
            {
                landmarks_in_range.push_back(LandmarkObs{lm.id_i, lm.x_f, lm.y_f});
            }
        }

        // convert the observations from vehicle to map coordinate system (assuming this particle as reference)
        vector<LandmarkObs> observations_map;

        for (const auto& obs: observations) {
            observations_map.push_back({
                0,
                obs.x * cos(p.theta) - obs.y * sin(p.theta) + p.x,
                obs.x * sin(p.theta) + obs.y * cos(p.theta) + p.y
            });
        }

        // associate each observation with the closest landmark
        dataAssociation(landmarks_in_range, observations_map);

        // calculate the particle weight, based on the distance between observations and associated landmark
        p.weight = 1.0;

        for (const auto& obs: observations_map) {
            const auto& associated_landmark = map_landmarks.landmark_list.at(obs.id-1);
            const double w = normal(obs.x, obs.y, associated_landmark.x_f, associated_landmark.y_f, std_landmark[0], std_landmark[1]);
            p.weight *= w;
        }
    }

}

void ParticleFilter::resample()
{
    /**
     * TODO: Resample particles with replacement with probability proportional
     *   to their weight.
     * NOTE: You may find std::discrete_distribution helpful here.
     *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
     */

    vector<double> weights(particles.size(), 0.0);
    for (size_t i = 0; i < particles.size(); i++) {
        weights[i] = particles[i].weight;
    }

    discrete_distribution<> dist(weights.begin(), weights.end());

    vector<Particle> resampled_particles;
    resampled_particles.resize(number_of_particles());

    for (size_t i = 0; i < resampled_particles.size(); i++) {
        const size_t idx = dist(gen);
        resampled_particles[i] = particles[idx];
    }

    particles = std::move(resampled_particles);
}

void ParticleFilter::SetAssociations(Particle& particle,
    const vector<int>& associations,
    const vector<double>& sense_x,
    const vector<double>& sense_y)
{
    // particle: the particle to which assign each listed association,
    //   and association's (x,y) world coordinates mapping
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates
    particle.associations = associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
    vector<int> v = best.associations;
    std::stringstream ss;
    copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord)
{
    vector<double> v;

    if (coord=="X") {
        v = best.sense_x;
    }
    else {
        v = best.sense_y;
    }

    std::stringstream ss;
    copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
