#include <filesystem>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
#include "../functions/functions.hpp"
#include "simulator_classes/inc/system.h"
using namespace std;
namespace fs = std::filesystem;

int main (int argc, char *argv[]){
  vector<string> phase={"solid","liquid","gas"};
  vector<string> temp={"0.8","1.1","1.2"};
  vector<string> rho={"1.1","0.8","0.05"};
  vector<string> R_cut={"2.2","2.5","5.0"};
  vector<string> ipath={"input_04.2/solid/", "input_04.2/liquid/", "input_04.2/gas/"};
  for (int i=0;i<3;i++){
    ostringstream temp_stream;
    string folderName ="Results/output_04.2/equilibrated_simulation/"+phase[i]+"_phase_temp_"+temp[i]+"_rho_"+rho[i]+"_Rcut_"+R_cut[i];
    string folderName1 =folderName+"/OUTPUT";
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
    // Creazione e apertura del file
    string input_path=ipath[i];
    string output_path= folderName+"/"; 
    int nconf = 1;
    System SYS;
    SYS.initialize(input_path,output_path); 
    SYS.initialize_properties(input_path,output_path);
    SYS.block_reset(0);
    for(int i=0; i < SYS.get_nbl(); i++){ //loop over blocks
      for(int j=0; j < SYS.get_nsteps(); j++){ //loop over steps in a block
        SYS.step();
        SYS.measure();
        if(j%10 == 0){
          nconf++;
        }
      }
      SYS.averages(output_path,i+1);
      SYS.block_reset(output_path,i+1);
    }
    SYS.finalize(output_path);
  }
  return 0;
}

