#include "../inc/func.h"
using namespace std;

//---------------------------
// FUNCTION CLASS
//---------------------------
//constructor
Function::Function() {}

//---------------------------
// PSI_100 CLASS
//---------------------------
//constructor
psi_100::psi_100() {}
//destructor
psi_100::~psi_100() {}
//method to calculate the function
double psi_100::operator()(arma::vec r) const {return pow(m_a0,-1.5)*pow(M_PI,-0.5)*exp(-r(0)/m_a0);}

//---------------------------
// PSI_210 CLASS
//---------------------------
//constructor
psi_210::psi_210() {}
//destructor
psi_210::~psi_210() {}
//method to calculate the function
double psi_210::operator()(arma::vec r) const {return pow(m_a0,-2.5)*0.125*pow(2./M_PI,0.5)*r(0)*exp(-(r(0))/(2*m_a0))*cos(r(1));}

//---------------------------
// norm_function CLASS  
//---------------------------
//constructor
norm2::norm2() {}
//destructor
norm2::~norm2() {}
//method to calculate the function
double norm2::operator()(arma::vec r) const {return arma::norm(r, 2);}


//---------------------------
// PSI_8 CLASS
//---------------------------
//constructor
psi_8::psi_8() {}
psi_8::psi_8(double mu, double sigma) : m_mu(mu), m_sigma(sigma) {}
//destructor
psi_8::~psi_8() {}
//method to calculate the function
double psi_8::operator()(arma::vec r) const {return exp(-0.5*pow(r(0)-m_mu,2)*pow(m_sigma,-2))+exp(-0.5*pow(r(0)+m_mu,2)*pow(m_sigma,-2));}


//---------------------------
// particle in a potential
//---------------------------
//constructor
Hpsi::Hpsi() {}
Hpsi::Hpsi(double mu, double sigma) : m_mu(mu), m_sigma(sigma) {}
//destructor
Hpsi::~Hpsi() {}
//method to calculate the function
double Hpsi::operator()(arma::vec r) const {
    double d2psi= 0.5*pow(m_hbar,2)/m_m*(((pow(m_sigma,2)-pow(r(0)-m_mu,2))*pow(m_sigma,-4))*exp(-0.5*pow(r(0)-m_mu,2)*pow(m_sigma,-2))+((pow(m_sigma,2)-pow(r(0)+m_mu,2))*pow(m_sigma,-4))*exp(-0.5*pow(r(0)+m_mu,2)*pow(m_sigma,-2)));
    double psi=exp(-0.5*pow(r(0)-m_mu,2)*pow(m_sigma,-2))+exp(-0.5*pow(r(0)+m_mu,2)*pow(m_sigma,-2));
    double V=pow(r(0),4)-2.5*pow(r(0),2);
    return d2psi/psi+V;
    //return d2psi;
    }
