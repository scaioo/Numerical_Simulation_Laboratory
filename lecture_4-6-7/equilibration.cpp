#include <iostream>
#include <filesystem>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include "../functions/functions.hpp"
#include "simulator_classes/inc/system.h"
using namespace std;
namespace fs = std::filesystem;

int main (int argc, char *argv[]){
  // Check if the number of arguments is correct
  if (argc != 5) {
    cerr << "Usage: " << argv[0] << " Rho R_cut T_start output_path" << endl;
    return 1;
  }
  // Parse command-line arguments
  string Rho = argv[1];
  string R_cut = argv[2];
  double temp = atof(argv[3]);
  string output_path = argv[4];
  
  change_file("input_equilibration/INPUT/input.dat","RHO",Rho);
  change_file("input_equilibration/INPUT/input.dat","R_CUT",R_cut);
  
  for (int i=0;i<5;i++){
    ostringstream temp_stream;
    temp_stream << fixed << setprecision(2) << temp;
    change_file("input_equilibration/INPUT/input.dat","TEMP",temp_stream.str());

    // Create directories' paths for the equilibration performed
    string folderName =output_path+"temp_"+temp_stream.str();
    string folderName1 =folderName+"/OUTPUT";
    string folderName2 =folderName1+"/CONFIG";
    try {
      if (fs::create_directory(folderName) and fs::create_directory(folderName1) and fs::create_directory(folderName2)) {
        cout << "Cartella creata con successo: " << folderName << endl;
        cout << "Cartella creata con successo: " << folderName1 << endl;
        cout << "Cartella creata con successo: " << folderName2 << endl;
      }else{
        cout << "La cartella esiste già o non può essere creata: " << folderName << endl;
        cout << "La cartella esiste già o non può essere creata: " << folderName1 << endl;
        cout << "La cartella esiste già o non può essere creata: " << folderName2 << endl;
      }
    }catch (const fs::filesystem_error& e) {cerr << "Errore: " << e.what() << endl;}
    // Create and open file
    string output_path= folderName+"/"; 
    System SYS;
    string input_path="input_equilibration/";
    SYS.initialize(input_path,output_path); 
    SYS.initialize_properties(input_path,output_path);
    SYS.block_reset(0);
    for(int i=0; i < SYS.get_nbl(); i++){ //loop over blocks
      for(int j=0; j < SYS.get_nsteps(); j++){ //loop over steps in a block
        SYS.step();
        SYS.measure();
      }
      SYS.averages(output_path,i+1);
      SYS.block_reset(output_path,i+1);
    }
    SYS.finalize(output_path);
    temp+=0.01;
  }
  return 0;
}

