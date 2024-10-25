#include<iostream>
#include"simulator_classes/inc/system.h"
#include "../functions/functions.hpp"
#include <filesystem>
#include <cstdlib>
using namespace std;
namespace fs = std::filesystem;
int main(int argc,char* argv[]){
    cout << "start" << endl;
    string input_path;
    string path="Results/output_07.2/simulation_equilibrated/";
    // define a string matrix with the parameters i want to use
    vector<vector<string>> parameters(4, vector<string>(3)); // 2D vector with 4 rows and 3 columns
    parameters[0] = {"solid","liquid", "gas"}; // phase
    parameters[1] = {"1.1", "2.2", "0.8"};     // solid
    parameters[2] = {"0.8", "2.5", "1.1"};     // liquid
    parameters[3] = {"0.05","5.0", "1.2"};     // gas
    cout << "parameters defined" << endl;
    for(int i=0; i<3; i++){
        // set the parameters
        string Rho = parameters[i+1][0];
        string R_cut = parameters[i+1][1];
        string Temp = parameters[i+1][2];
        change_file("input_07.2/"+parameters[0][i]+"/INPUT/input.dat","RHO",Rho);
        change_file("input_07.2/"+parameters[0][i]+"/INPUT/input.dat","R_CUT",R_cut);
        change_file("input_07.2/"+parameters[0][i]+"/INPUT/input.dat","TEMP",Temp);
        cout << "run simulation for phase: " << parameters[0][i] << " with Rho: " << Rho << " R_cut: " << R_cut << " Temp: " << Temp << endl;
        input_path="input_07.2/"+parameters[0][i]+"/";
        string output_path= path+parameters[0][i]+"/";
        // create directories' paths for the simulation performed
        string folderName =output_path;
        string folderName1 =folderName+"OUTPUT";
        string folderName2 =folderName1+"/CONFIG";
        try {
            if (fs::create_directory(folderName) and fs::create_directory(folderName1) and fs::create_directory(folderName2) ) {
              cout << "Cartella creata con successo: " << folderName << endl;
              cout << "Cartella creata con successo: " << folderName1 << endl;
              cout << "Cartella creata con successo: " << folderName2 << endl;
            }else{
              cout << "La cartella esiste già o non può essere creata: " << folderName << endl;
              cout << "La cartella esiste già o non può essere creata: " << folderName1 << endl;
              cout << "La cartella esiste già o non può essere creata: " << folderName2 << endl;
            }
        }catch (const fs::filesystem_error& e) {cerr << "Errore: " << e.what() << endl;}
        // run the simulation 
        cout << "output path: " << output_path << endl;
        System SYS;
        cout << "system created" << endl;
        SYS.initialize(input_path,output_path); 
        cout << "initialize done" << endl;
        SYS.initialize_properties(input_path,output_path);
        cout << "initialize properties done" << endl;
        SYS.findelta(0.5,0.01);
        cout << "findelta done" << endl;
        SYS.block_reset(0);
        cout << "block reset done" << endl;
        for(int i=0; i < SYS.get_nbl(); i++){ //loop over blocks
          for(int j=0; j < SYS.get_nsteps(); j++){ //loop over steps in a block
            SYS.step();
            SYS.measure();
          }
          SYS.averages(output_path,i+1);
          SYS.block_reset(output_path,i+1);
        }
        SYS.finalize(output_path);
        }
    return 0;
}
