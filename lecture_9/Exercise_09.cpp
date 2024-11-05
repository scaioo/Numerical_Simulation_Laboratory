#include<iostream>
#include<fstream> 
#include<armadillo>
#include <unordered_set> 
#include <utility>
#include"../classes/inc/genetic_algorithm_def.h"
#include "../random_number/random.h" 
using namespace std;
using namespace arma;

int main(){

    // -----------------------------------------------------------//
    //                        Circle                              //
    // -----------------------------------------------------------//

    cout << "-------------------------------------------------"<< endl;
    cout << "                   CIRCLE " << endl;
    cout << "-------------------------------------------------"<< endl; 
    
    // Variables declaration
    arma_rng::set_seed(12345);   
    Random rnd;  
    int seed1[4]={1234,4567,9876,5432};
    rnd.start();
    rnd.SetRandom(seed1,2958,1795);
    Population pop(rnd);  
    int indi_dim=34; 
    int pop_dim=250; 
    int iterations=500; 
    pop.initialize_circle(indi_dim, pop_dim,1., rnd);
    ofstream circle_losses; 
    string circle_losses_path="Results/circle_losses.csv";
    circle_losses.open(circle_losses_path);
    circle_losses << "Generation;\tMean Loss;\tBest Loss" << endl;
    mat A=pop[0].Get_positions();    

    // Save initial matrix in CSV format
    A.save("Results/matrix_circle_0.csv", arma::csv_ascii);

    // Loop over the generations
    for(int i=0;i<iterations;i++){  

        if(i%100==0) cout << "--------------------------------------------------------"<< endl;
        if(i%100==0) cout << "-------------- Generation " << i <<"--------------------"<< endl;
        if(i%100==0) cout << "--------------------------------------------------------"<< endl;

        pop.Selection();
        pop.Mutation(0.65, 0.1, 0.1, 0.1, 0.1);

        if(i%10==0){  
            // save results
            pop.sort();
            A=pop[0].Get_positions();
            A.save("Results/matrix_circle_"+to_string(i)+".csv", arma::csv_ascii);
            cout << "mean loss: "<<pop.get_pop_mean_Loss()<< "\tbest Loss: " << pop[0].get_indi_Loss()<< endl;
            circle_losses << i <<";"<<pop.get_pop_mean_Loss()<<";" << pop[0].get_indi_Loss()  << endl;
        }
    }

    cout << "-------------- Final Generation --------------------"<< endl;
    cout << "Best individual: ";
    pop.sort();
    // Print the best individual 
    pop[0].print_idxs();
    cout <<" loss: " <<pop[0].get_indi_Loss()<< endl; 
    A=pop[0].Get_positions();
    // Save matrix in CSV format
    A.save("Results/matrix_circle_"+to_string(iterations)+".csv", arma::csv_ascii);
    circle_losses.close();

    // -----------------------------------------------------------//
    //                        Square                              //
    // -----------------------------------------------------------//

    cout << "-------------------------------------------------"<< endl;
    cout << "                   SQUARE " << endl;
    cout << "-------------------------------------------------"<< endl; 

    // Variables declaration
    pop.initialize_square(indi_dim, pop_dim,1., rnd);
    ofstream square_losses;  
    string square_losses_path="Results/square_losses.csv";
    square_losses.open(square_losses_path);
    square_losses << "Generation;\tMean Loss;\tBest Loss" << endl;
    A=pop[0].Get_positions();   

    // Save initial matrix in CSV format
    A.save("Risultati/matrix_square_0.csv", arma::csv_ascii);

    // Loop over the generations
    for(int i=0;i<iterations;i++){  
        if(i%100==0) cout << "--------------------------------------------------------"<< endl;
        if(i%100==0) cout << "-------------- Generation " << i <<"--------------------"<< endl;
        if(i%100==0) cout << "--------------------------------------------------------"<< endl;
        
        pop.Selection();
        pop.Mutation(0.65, 0.1, 0.1, 0.1, 0.1);
         
        if(i%10==0){  
            // save results
            pop.sort();
            A=pop[0].Get_positions(); 
            A.save("Results/matrix_square_"+to_string(i)+".csv", arma::csv_ascii);
            cout << "mean loss: "<<pop.get_pop_mean_Loss()<< "\tbest Loss: " << pop[0].get_indi_Loss()<< endl;
            square_losses << i <<";"<<pop.get_pop_mean_Loss()<<";" << pop[0].get_indi_Loss()  << endl;
        }
    }

    cout << "-------------- Final Generation --------------------"<< endl;
    cout << "Best individual: ";
    
    // Print the best individual
    pop.sort(); 
    pop[0].print_idxs();
    std::cout <<" loss: " <<pop[0].get_indi_Loss()<< endl; 
    A=pop[0].Get_positions();
    // Save matrix in CSV format
    A.save("Results/matrix_square_"+to_string(iterations)+".csv", arma::csv_ascii);
    square_losses.close();
    return 0;
}
