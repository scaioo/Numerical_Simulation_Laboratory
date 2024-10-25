#include <iostream>
#include <cmath>
#include <armadillo>
using namespace std;    

#ifndef __func__
#define __func__
class Function{
    public:
    //constructor
    Function();
    //destructor
    virtual ~Function() = default;
    //method to calculate the function
    //virtual double operator()(double r,double theta,double phi) const = 0;
    virtual double operator()(arma::vec) const = 0;
};
#endif // __func__

#ifndef __psi_100__
#define __psi_100__
class psi_100 : public Function{
    public:
    //constructor
    psi_100();
    //destructor
    ~psi_100();
    //method to calculate the function
    //double operator()(double r,double theta,double phi) const override;
    double operator()(arma::vec) const override;
    private:
    double m_a0=1.;
};
#endif // __psi_100__

#ifndef __norm2__
#define __norm2__
class norm2 : public Function{
    public:
    //constructor
    norm2();
    //destructor
    ~norm2();
    //method to calculate the function
    //double operator()(double r,double theta,double phi) const override;
    double operator()(arma::vec) const override;
    private:
    arma::vec m_r;
};
#endif // __norm2__

#ifndef __psi_210__
#define __psi_210__
class psi_210 : public Function{
    public:
    //constructor
    psi_210();
    //destructor
    ~psi_210();
    //method to calculate the function
    //double operator()(double r,double theta,double phi) const override;
    double operator()(arma::vec) const override;
    private:
    double m_a0=1.;
};
#endif // __psi_210__

#ifndef __psi_8__
#define __psi_8__
class psi_8 : public Function{
    public:
        //constructor
        psi_8(double mu,double sigma);
        psi_8();
        //destructor
        ~psi_8();
        //method to calculate the function
        double operator()(arma::vec) const override;
        // set sigma
        void Set_sigma(double sigma){m_sigma=sigma;}
        // set mu
        void Set_mu(double mu){m_mu=mu;}
    private:
        double m_sigma=1.,m_mu=1.;
};
#endif // __psi_8__

#ifndef __Hpsi__
#define __Hpsi__
class Hpsi : public Function{
    public:
        //constructor
        Hpsi(double mu,double sigma);
        Hpsi();
        //destructor
        ~Hpsi();
        //method to calculate the function
        //double operator()(double r,double theta,double phi) const override;
        double operator()(arma::vec) const override;
        // set sigma
        void Set_sigma(double sigma){m_sigma=sigma;}
        // set mu
        void Set_mu(double mu){m_mu=mu;}
    private:
        double m_a=1.;
        double m_hbar=1.;
        double m_m=1.;
        double m_sigma=1.,m_mu=1.;
};
#endif // __Hpsi__
