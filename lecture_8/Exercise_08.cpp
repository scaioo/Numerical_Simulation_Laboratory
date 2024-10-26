#include<iostream>
#include<armadillo>
#include"../classes/inc/block_average.h"
#include"../classes/inc/func.h"
#include"../random_number/random.h"
using namespace std;
using namespace arma;

int main(int argc, char* argv[]){
    // start random number generator
    Random rnd;
    rnd.start();
    int dimension=2;
    int Nstep=10000;
    int start=5;
    vec x_start={1.},m_r={1.,1.},vec_rnew{0.,0.},psi_value{0.,0.};
    double initial_increment=0.3;
    double ratio_target=0.5;
    double tolerance=0.001;
    int Nblock=1000; int Blocksize=500;
    ofstream _file_traj,_file_H,_file_psi;
    static ofstream null_stream;
    string _string_traj="Results/traj.dat",_string_H="Results/H.dat",_string_psi="Results/psi.dat",filename="helper.dat",filename1="helper_1.dat",filename2="helper_2.dat";
    bool m_save_traj=true,m_save_H=true,m_save_psi=true,m_gaussian=true;
    psi_8 psi(1.,1.);
    Hpsi hpsi(1.,1.);
    double _beta=0.1,mu,sigma;
    do{
        cout << _beta << endl;
        if(m_save_traj) _file_traj.open(_string_traj,ios::app);
        if(m_save_H) _file_H.open(_string_H,ios::app);
        // Determine the dimensionality from m_r (mu and sigma) 
        double T = 1./_beta;
        vec vec_rand(dimension);
        if (m_gaussian) {
            for (int i = 0; i < dimension; ++i) {
                vec_rand(i) = rnd.Gauss(0, T);
            }
        } else {
            for (int i = 0; i < dimension; ++i) {
                vec_rand(i) = rnd.Rannyu(-T, T);
            }
        }
        vec r_new=m_r+vec_rand;
        vec H={0.,0.};
        // calculate the meanvalue of the energy
        mu=m_r(0); sigma=m_r(1);
        psi.Set_mu(mu); psi.Set_sigma(sigma);
        hpsi.Set_mu(mu); hpsi.Set_sigma(sigma); 
        if(_beta>=1000){
            _file_psi.open(_string_psi);
            _file_psi.close();
            Metropolis metro(psi,x_start,_file_psi,hpsi);
            metro.Set_string(_string_psi);
            metro.Set_a(Nstep,start,initial_increment,ratio_target,tolerance,rnd,dimension);
            metro.Average(Nblock,Blocksize,filename);
            H=metro.Get_average();
        }else{
            Metropolis metro(psi,x_start,hpsi);
            metro.Set_a(Nstep,start,initial_increment,ratio_target,tolerance,rnd,dimension);
            metro.Average(Nblock,Blocksize,filename);
            H=metro.Get_average();
        }

        //modify beta (mu and sigma)
        mu=r_new(0); sigma=r_new(1);
        psi.Set_mu(mu); psi.Set_sigma(sigma);
        hpsi.Set_mu(mu); hpsi.Set_sigma(sigma);
        Metropolis metro_new(psi,x_start,hpsi);
        metro_new.Set_a(Nstep,start,initial_increment,ratio_target,tolerance,rnd,dimension);
        metro_new.Average(Nblock,Blocksize,filename2);
        vec H_new=metro_new.Get_average();

        double p=exp(-_beta*H(0));
        double pnew=exp(-_beta*H_new(0));
        double A = min(1.,pnew/p);
        double accept = rnd.Rannyu(0,1);
        if (accept <= A) {
                _beta++;
                //update m_r
                m_r = r_new;
                // save measures on files
                if (m_save_traj) {
                    for (int i = 0; i < dimension; ++i) {
                         _file_traj << m_r(i) << (i == dimension - 1 ? "\n" : "\t");
                    }
                    _file_traj.close();
                }
                if (m_save_H) {
                    for (int i = 0; i < 2; ++i) {
                         _file_H << H(i) << (i == dimension - 1 ? "\n" : "\t");
                    }
                    _file_H.close();
                }
            } else {
                if (m_save_traj) _file_traj.close();
                if (m_save_H) _file_H.close();
            }
    }while(_beta<=1001);
    return 0;
}