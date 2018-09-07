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
		//add noise to velocity and yaw_rate, comment the next 2 lines down to make it noiseless
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

	/*
		Still don't know what exactly this association does, however, 
			it does associate what the particle thinks it sees with what it observes,
			if we do that with all particles, that would n^3 runtime.. baaad..
	*/
	double currdist;
	int id_pred;
	
	for ( int i = 0; i < predicted.size(); i++){
		double mindist = 999.0;
		
		for ( int j=0; j < observations.size(); j++){
			currdist = dist(observations[j].x, observations[j].y, predicted[i].x, predicted[i].y);
			if ( currdist < mindist){
				mindist = currdist;
				id_pred = i;
			}
		}

		/*
		//Notice the association made in the printline, may put that in some sort of map soon

		cout<< "The closest to: " << predicted[i].id << " " << predicted[i].x << " " << predicted[i].y  << " is "
								  << observations[id_pred].id << " " << observations[id_pred].x << " " << observations[id_pred].y << endl; 
		*/
	}
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
	///////////////////////////////////////////////////////////////////////////////////////
		for( int p_index = 0; p_index < particles.size(); p_index++){
	
	Particle p = particles[p_index];	
	double particle_weight = 1.0;
	
	//observations
	for ( int i=0; i < observations.size(); i++){
		
		double x_observe_c = observations[i].x;
		double y_observe_c = observations[i].y;

		double xm = p.x + cos(p.theta) * x_observe_c  + -1 * sin ( p.theta) * y_observe_c; 
		double ym = p.y + sin(p.theta) * x_observe_c + cos(p.theta) * y_observe_c;

		double minsize = 9999;
		int id_curr = 9999;
		for( int l=0; l < map_landmarks.landmark_list.size(); l++){
			
			double obs_to_lm = dist(map_landmarks.landmark_list[l].x_f, map_landmarks.landmark_list[l].y_f, xm, ym);
			
			if ( obs_to_lm < sensor_range && obs_to_lm < minsize ){
				minsize = obs_to_lm;
				id_curr = l;
			}			
		}

		//cout << "Observation " << i << " " << x_observe_c << "," << y_observe_c << " most likely landmark " << landmarks[id_curr].id  << endl; 
	
		/* Gauss Norm for weight update*/
		double sig_x = std_landmark[0]; 
		double sig_y = std_landmark[1];
		double x_obs = xm;
		double y_obs = ym;
		double mu_x = map_landmarks.landmark_list[id_curr].x_f;
		double mu_y = map_landmarks.landmark_list[id_curr].y_f;

		// calculate exponent
		double exponent = ( pow((x_obs - mu_x),2 ) )/(2. * pow(sig_x,2) ) + (pow((y_obs - mu_y),2 )/(2. * pow(sig_y,2)));
		
		// calculate normalization term
		double gauss_norm = (1./(2. * 3.141592653 * sig_x * sig_y));

		// calculate weight using normalization terms and exponent
		double weight_curr = gauss_norm * exp(-1 * exponent);
		//cout<< i << " " << weight_curr << endl;
		particle_weight = particle_weight * weight_curr;
	}

	particles[p_index].weight = particle_weight;	

	}
	/////////////////////////////////////////////////////////////    __								__	__ _
	/*  Normalizing particles weights--> otherwise the resampling breaks (probably) */

	double accumulator = 0.0;
	for(int i=0;i< particles.size();i++){
		accumulator += particles[i].weight;
	}
	for(int i=0;i< particles.size();i++){
		particles[i].weight = particles[i].weight / accumulator;
	}

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	//start with std::vector<Particle> particles and num_particles
	vector<Particle>result;
	//outer ring is a cumulative distribution function, cdf
	double c[num_particles];
	c[0] = particles[0].weight;

	for( int k = 1; k < num_particles; k++){
		//think of c's as what percentage around the clock have you gone?
		c[k] = c[k-1] + particles[k].weight;
	}

	// make a random number between 0 and 1 / num_particles
	random_device rd;
	default_random_engine generator(rd()); // rd() provides a random seed
	uniform_real_distribution<double> distribution(0.0, 1./num_particles);	
	double u1 = distribution(generator);

	double u[num_particles];
	u[0] = u1;

	int i = 0;
	for( int j = 0; j < num_particles; j++){
		//think of c's as what percentage around the clock have you gone?
		while( u[j] > c[i] ){
			i += 1;
		}
		// The line below is the link that this function based upon, question is, should we re-do particle's weight?
		//https://classroom.udacity.com/courses/ud810/lessons/3353208568/concepts/33538586070923#
		particles[i].weight = 1.0 / num_particles;
		result.push_back( particles[i] );
		u[j + 1] = u[j] + 1.0 / num_particles;
	}

	particles = result;
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
    return particle;
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
