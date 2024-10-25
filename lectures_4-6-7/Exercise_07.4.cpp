#include<iostream>
#include"simulator_classes/inc/system.h"
#include "../functions/functions.hpp"
#include <filesystem>
#include <cstdlib>
using namespace std;
namespace fs = std::filesystem;

int main(int argc,char*argv[]){
    string input_path;
    string output_path;
    string path="Results/output_07.4/";
    // define a string matrix with the parameters i want to use
    vector<vector<string>> parameters(4, vector<string>(3)); // 2D vector with 4 rows and 3 columns
    parameters[0] = {"solid","liquid", "gas"}; // phase
    parameters[1] = {"1.1", "2.2", "0.8"};     // solid:  rho rcut T
    parameters[2] = {"0.8", "2.5", "1.1"};     // liquid: rho rcut T
    parameters[3] = {"0.05","5.0", "1.2"};     // gas:    rho rcut T
    cout << "parameters defined" << endl;
    //cycle over the NVT/NVE simualtion
    for(int j=0;j<2;j++){
    //cycle over the phase
    for(int i=0; i<3; i++){
        // set the parameters
        string Rho = parameters[i+1][0];
        string R_cut = parameters[i+1][1];
        string Temp = parameters[i+1][2];
        if(j==0) {
            input_path="input_07.4/NVE/"+parameters[0][i]+"/";
            output_path= path+"NVE/"+parameters[0][i]+"/";
            change_file("input_07.4/NVE/"+parameters[0][i]+"/INPUT/input.dat","RHO",Rho);
            change_file("input_07.4/NVE/"+parameters[0][i]+"/INPUT/input.dat","R_CUT",R_cut);
            change_file("input_07.4/NVE/"+parameters[0][i]+"/INPUT/input.dat","TEMP",Temp);
            change_file("input_07.4/NVE/"+parameters[0][i]+"/INPUT/input.dat","SIMULATION_TYPE",to_string(j));
            change_file("input_07.4/NVE/"+parameters[0][i]+"/INPUT/input.dat","NSTEPS","1000");
            change_file("input_07.4/NVE/"+parameters[0][i]+"/INPUT/input.dat","NBLOCKS","100");
        }else{
            input_path="input_07.4/NVT/"+parameters[0][i]+"/";
            output_path= path+"NVT/"+parameters[0][i]+"/";
            change_file("input_07.4/NVT/"+parameters[0][i]+"/INPUT/input.dat","RHO",Rho);
            change_file("input_07.4/NVT/"+parameters[0][i]+"/INPUT/input.dat","R_CUT",R_cut);
            change_file("input_07.4/NVT/"+parameters[0][i]+"/INPUT/input.dat","TEMP",Temp);
            change_file("input_07.4/NVT/"+parameters[0][i]+"/INPUT/input.dat","SIMULATION_TYPE",to_string(j));
            change_file("input_07.4/NVT/"+parameters[0][i]+"/INPUT/input.dat","NSTEPS","1000");
            change_file("input_07.4/NVT/"+parameters[0][i]+"/INPUT/input.dat","NBLOCKS","100");
                }
        cout << "run simulation for phase: " << parameters[0][i] << " with Rho: " << Rho << " R_cut: " << R_cut << " Temp: " << Temp << " Sim type: "<< j << endl;
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
        cout << "start simulation" << endl;
        cout << "input path: " << input_path << endl;
        cout << "output path: " << output_path << endl;
        System SYS;
        cout << "system created" << endl;
        SYS.initialize(input_path,output_path); 
        cout << "system initialized" << endl;
        SYS.initialize_properties(input_path,output_path);
        cout << "properties initialized" << endl;
        if(j==1)SYS.findelta(0.5,0.01);
        SYS.block_reset(0);
        cout << "start loop simulation" << endl;
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
    }
    return 0;
}