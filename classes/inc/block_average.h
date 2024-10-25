#include <fstream>
#include <string>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include "../../functions/functions.hpp"
#include "../../random_number/random.h"
#include "../inc/func.h"
using namespace std;
using namespace arma;

#ifndef __blockave__
#define __blockave__
class BlockMeanCalculator {
    public:
    //constructor
    BlockMeanCalculator();
    // Pure virtual method to calculate a quantity
    virtual double calculate(double value,Random &rnd) const = 0;
    // Method to compute block average and its uncertainty
    void Average(const int N_Blocks,const int BlockSize,string filename);
    //default destructor
    virtual ~BlockMeanCalculator();
    // Method to get the average
    vec Get_average(){return _average;};
    private:
    vec _average;
};
#endif // __blockave__

#ifndef __blockave_cop__
#define __blockave_cop__
// ---------------------------------------
// CALL OPTION PRICE COMPUTATION EX 3.1
// ---------------------------------------
class call_price : public BlockMeanCalculator {
    public:
    //constructor
    call_price(double S0, double T, double K, double r, double sigma, double mu, bool direct);
    // method to set Nstep
    void Set_Nstep(int step);
    // method to get Nstep
    double Get_Nstep();
    // Method to calculate the call option price
    double calculate(double value,Random &rnd) const override;
    //default destructor
    ~call_price();
    
    private:
    double m_S0;
    double m_T;
    double m_K;
    double m_r;
    double m_sigma;
    double m_mu;
    bool m_direct;
    int m_Nstep;
};
#endif // __blockave_cop__

// ---------------------------------------
// PUT OPTION PRICE COMPUTATION EX 3.1
// ---------------------------------------
#ifndef  __blockave_pot__
#define __blockave_pot__
class put_price : public BlockMeanCalculator {
    public:
    //constructor
    put_price(double S0, double T, double K, double r, double sigma, double mu, bool direct);
    //method to set Nstep
    void Set_Nstep(int step);
    //method to get Nstep
    double Get_Nstep();
    //method to calculate the put option price
    double calculate(double value,Random &rnd) const override;
    //destructor
     ~put_price();

    private:
    double m_S0;
    double m_T;
    double m_K;
    double m_r;
    double m_sigma;
    double m_mu;
    bool m_direct;
    int m_Nstep;
};
#endif  //__blockave_pot__

// ---------------------------------------
// FIRST METROPOLIS ALGORITHM EX 5.1
// ---------------------------------------
#ifndef __metropolis__
#define __metropolis__
    class Metropolis : public BlockMeanCalculator  {
    public:
    //constructor
    Metropolis(Function &f,vec &r_start,Function &g): 
        m_f(f), m_Nstep(0), m_a(0), m_r(r_start),m_rnew(r_start),m_gaussian(false), m_spherical_coordinates(false), m_save_traj(false),_file(null_stream), m_g(g) {};

    Metropolis(Function &f,vec &r_start,ofstream &file,Function &g): 
        m_f(f), m_Nstep(0), m_a(0), m_r(r_start),m_rnew(r_start),m_gaussian(false), m_spherical_coordinates(true), m_save_traj(true), _file(file), m_g(g) {};
    //method to calculate metropolis
    double calculate(double value,Random &rnd) const override;
    //method to equilibrate the system
    void Equilibrate(int Nstep,Random &rnd);
    //method to set Nstep
    void Set_Nstep(int step);
    //method to get Nstep
    double Get_Nstep();
    //method to get a
    double Get_a();
    //method to set a
    void Set_a(int Nstep,double x_start ,double initial_increment,double target_ratio,double tolerance,Random &rnd,int dimension);
    // method to set the string
    void Set_string(string filename){_string=filename;};
    // method to set the gaussian
    void Set_gaussian(bool gaussian){m_gaussian=gaussian;}
    // method to set the spherical coordinates
    void Set_spherical_coordinates(bool spherical){m_spherical_coordinates=spherical;}
    //destructor
    ~Metropolis();
    
    private:
    Function &m_f; // probability distribution
    double m_Nstep;
    double m_a; // length of the step
    vec &m_r;vec &m_rnew;
    bool m_gaussian;
    bool m_spherical_coordinates;
    bool m_save_traj;
    string _string;
    ofstream &_file;
    static ofstream null_stream; // Null stream to handle default file operations when no file is provided
    Function &m_g; // Observable to measure
};
#endif // __metropolis__
/*
#ifndef __Simulated_Annealing__
#define __Simulated_Annealing__
class Simulated_Annealing : public Metropolis {
    public:
    //constructor
    Simulated_Annealing(Function &f,vec &r_start,ofstream &file_traj,ofstream &file_H,ofstream &file_psi,Function &g): 
        m_f(f),m_g(g),m_r(r_start), _file_traj(file_traj),_file_H(file_H),_file_psi(file_psi), m_save_traj(true), m_save_H(true), m_save_psi(true), m_gaussian(true) {};
    //method to calculate simulated annealing
    double calculate(double value,Random &rnd) const override;
    void calc(double beta,Random &rnd) const;
    //destructor
    ~Simulated_Annealing();
    // functions to set and get values:
    void set_gaussian(bool gaussian){m_gaussian=gaussian;}

    private:
    double _beta=0.;
    bool m_save_traj,m_save_H,m_save_psi;
    vec &m_r;
    ofstream &_file_traj, &_file_H, &_file_psi;
    string _string_traj, _string_H, _string_psi;
    Function &m_f; Function &m_g; 
    bool m_gaussian;
    // Null stream to handle default file operations when no file is provided
    static ofstream null_stream;

};

#endif // __Simulated_Annealing__
*/