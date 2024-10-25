#include <iostream>
#include <iomanip>
#include <fstream>
#include <numeric>
#include <vector>
#include <cmath>
#include "../random_number/random.h"
#include "../functions/functions.hpp"
using namespace std;

int main() {
    //setting the random number generator
    Random rnd;
    rnd.start();
    //declaration of variables
    int N_blocks=100;                                          //number of throws and blocks
    int L=1000000;                                            //number of thorws per block
    long double d=1.1;                                              //distance between two lines
    long double y=0.;                                                //position of the beginning of the needle
    long double l=1.;                                               //length of the needle       //random number generators
    long double cosphi=0.,ly=0.,lx=0.;                                 //cosine of the angle between the needle and the y-axis and the x and y components to generate a random angle
    long double pi=0.,pi_helper=0.;                                   //value of pi and helper variable
    int Nhit=0;                                                //number of intersections

    // blocking average variables
    long double ave=0.;
    long double av2=0.;
    long double sum_prog=0.;
    long double su2_prog=0.;
    long double sum_prog_helper=0.;
    long double su2_prog_helper=0.;
    long double err_prog=0.;

    // open a file to save the results
    ofstream outputfile=openofile("Results/ex_01.3_buffon.txt");
    // writing the first line of the file
    outputfile << "Number_block;\tpigreco;\terror" << endl;
    //throwing the needle
    for(int i=0;i<N_blocks;i++){
        Nhit=0;
        for(int j=0;j<L;j++){
            //generating the position of the end of the needle
            y=rnd.Rannyu(0,d);
            //generating the angle of the needle
            do{
                lx=rnd.Rannyu(-1,1);
                ly=rnd.Rannyu(-1,1);
            }while((pow(lx,2)+pow(ly,2))>1);
            //calculating the cosine of the angle
            cosphi=ly/sqrt(pow(lx,2)+pow(ly,2));
            //checking if the needle intersects a line
            if((y<d && d<=l*cosphi+y) || (y>0 && y+l*cosphi<=0)){Nhit++;}
        }
        //calculating the value of pi
        pi_helper+=2*l*(long double)L/(d*Nhit);
        pi=pi_helper/(i+1);

        //blocking average
        ave=pi;
        av2=pow(ave,2);
        sum_prog_helper+=ave;                 // Sum j{0..i} r_j
        su2_prog_helper+=av2;                 // Sum j{0..i} (r_j)^2
        sum_prog=sum_prog_helper/(i+1);       // Cumulative average
        su2_prog=su2_prog_helper/(i+1);       // Cumulative square average
        err_prog=error(sum_prog,su2_prog,i);  // Statistical uncertainty

        // write the results on the file
        outputfile << i << ";\t" << sum_prog << ";\t" << err_prog << endl;
    }
    return 0;
}