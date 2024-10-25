#ifndef FUNCTION_H
#define FUNCTION_H
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>
#include <algorithm>
#include <sstream>
#include <armadillo>
#include "../random_number/random.h"
using namespace std;
using namespace arma;

double error(double AV,double AV2,int n);
/**/
ofstream openofile(const string &filename);
/**/
void changeSimType(const std::string& filePath, const std::string& newSimulationType);
/**/
void change_file(const string &filename,string keyword,string new_value);
/**/
pair<double, double> calculateMeanAndStdDev(const double accumulate, const double sum_prog_helper,const double su2_prog_helper, const int dim, const int block);
/**/
vec discrete_walk(Random &rnd,const vec &p, const double a);
/**/
vec continuous_walk(Random &rnd, const vec &p, const double a);
/**/
void print(const vector<double> &p,ofstream &ofile);
/**/
double distance2(const vec &p);
/**/
void print(const vec &p);
/**/
double spot_price(const double S0, const double mu, const double sigma, const double t,const double rnd);
/**/
double spot_price(const double S0, const double mu, const double sigma, const double t,Random &rnd,const int N_steps);
/**/
double call_option_price(const double r, const double T, const double K, const double S);
/**/
double put_option_price(const double r, const double T, const double K, const double S);
/**/
void cartesianToSpherical(const vec& cartesian, vec& spherical);
/**/
vector<double> getMeasure(const string& filename);
/**/
void change_Mag_field(const string& filename, double newValue);
/**/
vector<string> StringVec(double start, double end, double increment);
/**/
int pbc(int position, int size);
/**/

#endif // FUNCTION_H