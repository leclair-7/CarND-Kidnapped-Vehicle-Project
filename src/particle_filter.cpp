/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

#define EPSILON .001
using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 100;
	std_vel = 1.0; 
	std_yaw = 1.0;
	double std_x, std_y, std_theta; // Standard deviations for x, y, and theta
	// TODO: Set standard deviations for x, y, and theta
	std_x = std[0];
	std_y = std[1];
	std_theta = std[2];

	normal_distribution<double> dist_x(x, std_x);
	normal_distribution<double> dist_y(y, std_y);
	normal_distribution<double> dist_theta(theta, std_theta);

	default_random_engine gen;
	Particle p;

	for(int i=0; i < num_particles;i++)
	{
		double sample_x, sample_y, sample_theta;
		sample_x = dist_x(gen);
		sample_y = dist_y(gen);
		sample_theta = dist_theta(gen);

		p.id = i;
		p.x = sample_x;
		p.y = sample_y;
		p.theta = sample_theta;
		p.weight = 1;
		particles.push_back( p );
	}

	is_initialized = true;

	/* This tests init when running on the sim, checks if particles were initialized
	for(int i=0; i < num_particles;i++)
	{
		
		cout<< i << "\t" << particles[i].id << endl;
		cout<< i << "\t" << particles[i].x << endl;
		cout<< i << "\t" << particles[i].y << endl;
		cout<< i << "\t" << particles[i].theta << endl;
		
	}
	 string frog;
	 cout<< "Got to the frog line\n";
	 cin >> frog;
	 */
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	
	/*
	cout<< "Pre-predict" << endl;
	for(int i=0; i < 3;i++)
	{
		
		cout<< i << "\t" << particles[i].id << endl;
		cout<< i << "\t" << particles[i].x << endl;
		cout<< i << "\t" << particles[i].y << endl;
		cout<< i << "\t" << particles[i].theta << endl;
		
	}
	*/

	double v_over_theta;
	double d_theta_dt;
	double theta_0;
	double xf, yf, thetaf;

	/* Noise */
	double std_x, std_y, std_theta; // Standard deviations for x, y, and theta
	// TODO: Set standard deviations for x, y, and theta
	std_x = std_pos[0];
	std_y = std_pos[1];
	std_theta = std_pos[2];	
	default_random_engine gen;

	normal_distribution <double> dist_velocity(velocity, std_vel);
	normal_distribution <double> dist_yawd(yaw_rate, std_yaw);	
	
	for(int i=0; i < num_particles; i++)
	{
		//add noise to velocity and yaw_rate
		velocity = dist_velocity(gen); 
		yaw_rate = dist_yawd(gen);

		double x = particles[i].x;
		double y = particles[i].y;
		theta_0  = particles[i].theta;
		
		normal_distribution<double> dist_x(x, std_x);	
		normal_distribution<double> dist_y(y, std_y);
		normal_distribution<double> dist_theta(theta_0, std_theta);

		x  		= dist_x(gen);
		y 		= dist_y(gen);
		theta_0 = dist_theta(gen);
		
		if ( theta_0 < EPSILON )
		{
			xf     = x + velocity * delta_t * cos(theta_0);
			yf     = y + velocity * delta_t * sin(theta_0);
			thetaf = theta_0;

		} else{

			d_theta_dt   = yaw_rate * delta_t;
			v_over_theta = velocity / yaw_rate;

			xf     = x + v_over_theta * (sin(theta_0 + d_theta_dt) - sin( theta_0 )); 
			yf     = y + v_over_theta * (cos(theta_0 )             - cos( theta_0 + d_theta_dt )); 
			thetaf = theta_0 + v_over_theta;

		}

		particles[i].x     = xf;
		particles[i].y     = yf;
		particles[i].theta = thetaf;

	}

	/*
	cout<< "Post Predict" << endl;
	for(int i=0; i < 3;i++)
	{
		
		cout<< i << "\t" << particles[i].id << endl;
		cout<< i << "\t" << particles[i].x << endl;
		cout<< i << "\t" << particles[i].y << endl;
		cout<< i << "\t" << particles[i].theta << endl;
		
	}
	 string frog;
	 cout<< "Got to the frog line\n";
	 cin >> frog;
	 */
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	

	/* Beginning of the coordinate transformation which does rotation and translation via */
	/* this gets us from particle coordinates to car's map coordinates 
	    __								__	__ _
	xm	|cos(theta)		-sin(theta)		xp| |xc |
	ym	|sin(theta)		cos(theta)		yp| |yc |
	1	|1				1				1 |	|1  |
		--								--- -- --
	*/	
	/* End of the coordinate transformation */

	/* Gauss Norm for weight update*/
	double sig_x, sig_y, x_obs, y_obs, mu_x, mu_y, gauss_norm, exponent;
	
	/*
	sig_x= 0.3
	sig_y= 0.3
	x_obs= 6
	y_obs= 3
	mu_x= 5
	mu_y= 3
	*/

	// calculate normalization term
	gauss_norm= (1./(2. * 3.141592653 * sig_x * sig_y));

	// calculate exponent
	exponent= ( pow((x_obs - mu_x),2 ) )/(2. * pow(sig_x,2) ) + (pow((y_obs - mu_y),2 )/(2. * pow(sig_y,2)));

	// calculate weight using normalization terms and exponent
	double weight= gauss_norm * exp(-1 * exponent);
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
