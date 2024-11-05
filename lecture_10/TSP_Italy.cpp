#include<iostream>
#include<fstream>
#include<armadillo>
#include <unordered_set>
#include <utility>
#include <stdio.h>
#include <mpi.h>
#include"../classes/inc/genetic_algorithm_def.h"
#include "../random_number/random.h"

using namespace std;
using namespace arma;

int main(int argc, char *argv[]) { 

    arma::mat coordinate=get_coordinates("cap_prov_ita.dat");

        // ----------------------------- //
    // SET PARALLELIZATION VARIABLES //
    // ----------------------------- //
    int size, rank;
    #define MASTER 0
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    if (rank > 3 ) {
    printf("Need 4 MPI tasks, not %d. Quitting...\n",rank);
    MPI_Abort(MPI_COMM_WORLD, argc);
    exit(1);
    }
    MPI_Status stat1, stat2;
    MPI_Request req;
    // ------------------------------- //
    // SET GENETIC ALGORITHM VARIABLES //
    // ------------------------------- // 
    Random rnd;
    rnd.start();    
    Population pop(rnd);
    int indi_dim=coordinate.n_rows;
    int pop_dim=700;
    int iterations=1550; 
    vector<int> buffer;
    vector<int> global_buffer;
    bool migration=false;
    int seed1[4]={1234,4567,9876,5432};
    int seed2[4]={1423,6764,9385,2557};
    int seed3[4]={8238,7394,9865,5231};
    int seed4[4]={5678,2458,8769,7068};
    
    // Check if the feature is enabled
    for (int i = 1; i < argc; ++i) {  // argv[0] is the name of the program
        if(string(argv[i]) == "migration"){
            migration = true;
        }
    }


    // File to save the square losses
    ofstream losses;
    string losses_path;
    losses_path="Results/city_losses_migration_"+to_string(migration)+".csv";
    losses.open(losses_path);
    losses << "Generation;\tMean Loss;\tBest Loss" << endl;

    // Initialization     
    pop.initialize_city(indi_dim, pop_dim,coordinate,rnd);
    pop.checkPopulation();
    mat A=pop[0].Get_positions();
    if(rank==MASTER) rnd.SetRandom(seed1,2893,195);
    if(rank==1) rnd.SetRandom(seed2,2892,3659);
    if(rank==2) rnd.SetRandom(seed3,2894,3575);
    if(rank==3) rnd.SetRandom(seed4,2895,1097);
 
    // Save matrix in CSV format
    A.save("Results/matrix_circle_0.csv", arma::csv_ascii);
    
    for(int i=0;i<iterations;i++){   
        MPI_Barrier(MPI_COMM_WORLD);
        pop.Selection();
        pop.Mutation(0.65, 0.1, 0.2, 0.1, 0.2);
        pop.checkPopulation();
        if(i%50==0 and rank==MASTER){
            pop.sort();
            cout << "-------------- Generation " << i <<"--------------------"<< endl;
            cout << "Mean Loss: " << pop.get_pop_mean_Loss() << " Best Loss: " << pop[0].get_indi_Loss() << endl;
            losses << i <<";"<<pop.get_pop_mean_Loss()<<";" << pop[0].get_indi_Loss()  << endl;
            A=pop[0].Get_positions();
            A.save("Results/matrix_city_"+to_string(i)+"migration_"+to_string(migration)+".csv", arma::csv_ascii);
            } 
        if(i%50==0 and migration){

            buffer.clear();
            global_buffer.clear(); 
            pop.sort();
            for(int j=0;j<indi_dim;j++){buffer.push_back(pop[0].get_idxs(j));}
            global_buffer.resize(size*buffer.size());
            // Gather the buffer
            MPI_Gather(&buffer[0], buffer.size(), MPI_INT, &global_buffer[0], buffer.size(), MPI_INT, MASTER, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
            
            if(rank==MASTER){
                // Create a vector of indexes
                ivec idxs(size);
                for(int j=0;j<size;j++){idxs(j)=j;}
                idxs=shuffle(idxs);
                // Create a matrix to store the best individuals
                vector<vector<int>> matrix(size, vector<int>(indi_dim));

                // Store the best individuals in the matrix with the indexes shuffled
                for(int j=0;j<size;j++){
                    for(int k=0;k<indi_dim;k++){
                        matrix[idxs[j]][k]=global_buffer[j*indi_dim+k];
                    }
                }
                // Reshape the matrix into global buffer
                for(int j=0;j<size;j++){
                    for(int k=0;k<indi_dim;k++){
                        global_buffer[j*indi_dim+k]=matrix[j][k];
                    }
                }
            }
            // Broadcast the global buffer
            MPI_Bcast(&global_buffer[0], global_buffer.size(), MPI_INT, MASTER, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
            // Update the population
            if(rank==MASTER)for(int j=0;j<indi_dim;j++){pop[0].set_idxs(j,global_buffer[j]);}
            if(rank==1)for(int j=0;j<indi_dim;j++){pop[0].set_idxs(j,global_buffer[indi_dim+j]);}
            if(rank==2)for(int j=0;j<indi_dim;j++){pop[0].set_idxs(j,global_buffer[2*indi_dim+j]);}
            if(rank==3)for(int j=0;j<indi_dim;j++){pop[0].set_idxs(j,global_buffer[3*indi_dim+j]);}
        }
    }
    cout << "-------------- Final Generation --------------------"<< endl;
    cout << "Best individual: ";
    pop[0].print_idxs();
    cout << endl;
    A=pop[0].Get_positions();
    // Save the matrix in CSV format
    cout <<"Rank: " <<rank <<" Mean Loss: " << pop.get_pop_mean_Loss() << " Best Loss: " << pop[0].get_indi_Loss() << endl;
    A.save("Risultati/matrix_square_"+to_string(iterations)+".csv", arma::csv_ascii);
    losses.close();
    MPI_Finalize();
    return 0;  
}

  
