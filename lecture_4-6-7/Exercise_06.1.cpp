#include <filesystem>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
#include "../functions/functions.hpp"
#include "simulator_classes/inc/system.h"
#include <sstream>
using namespace std;
namespace fs = std::filesystem;

int main(int argc, char *argv[]){

    // Variables declaration
    vector<string> temp_string;
    vector<bool> Mag_field{true,false};
    vector<string> field{"Mag_field_on","Mag_field_off"};
    vector<bool> Metro{true,false};
    vector<string> met_gibb{"Metropolis","Gibbs"};
    vector<string> Mfield{"0.02","0.00"};
    temp_string = StringVec(0.1, 3.0, 0.01);
    // Set the restart parameter to 0
    change_file("input_06.1/INPUT/input.dat","RESTART","0");
    change_file("input_06.1/INPUT/input.dat","NBLOCKS","200");
    change_file("input_06.1/INPUT/input.dat","NSTEPS","1000");
    cout << "equilibration phase" << endl;
    // cycle over the magnetic field (on/off)
    for(int i=0; i<2; i++){
    //cycle over the simtype (metropolis/gibbs)
    for(int k=0; k<2; k++){
    //cycle over the temperature
    cout << "Metro=" << Metro[k] << " Field="<< Mfield[i] << endl;
    for (const string& str : temp_string){
        // set temperature
        change_file("input_06.1/INPUT/input.dat","TEMP",str);                                                                                  
        if(Metro[k]){changeSimType("input_06.1/INPUT/input.dat","2\t1.0\t"+Mfield[i]);}else{changeSimType("input_06.1/INPUT/input.dat","3\t1.0\t"+Mfield[i]);} // set simulation type 
        // set magnetic field
        //if(Mag_field[i]){change_Mag_field("input_06.1/INPUT/input.dat",0.02);}else{change_Mag_field("input_06.1/INPUT/input.dat",0.0);} 
        // Create directories' paths
        string folderName ="Results/output_06.1/equilibration/"+met_gibb[k]+"/"+field[i]+"/temp_"+str;
        string folderName1 =folderName+"/OUTPUT";
        string folderName2 =folderName1+"/CONFIG";
        // Create directories
        try {
            if (fs::create_directory(folderName) and fs::create_directory(folderName1) and fs::create_directory(folderName2) ) {
                //cout << "Cartella creata con successo: " << folderName << endl;
                //cout << "Cartella creata con successo: " << folderName1 << endl;
                //cout << "Cartella creata con successo: " << folderName2 << endl;
            }else{
                //cout << "La cartella esiste già o non può essere creata: " << folderName << endl;
                //cout << "La cartella esiste già o non può essere creata: " << folderName1 << endl;
                //cout << "La cartella esiste già o non può essere creata: " << folderName2 << endl;
            }
        }catch (const fs::filesystem_error& e) {cerr << "Errore: " << e.what() << endl;}
        // Create and open file
        string input_path="input_06.1/";
        string output_path= folderName+"/";
        System SYS;
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
    }
    } 
    }

    //------------------------------------------------------------//
    // After equilibration has been made restart every simulation //
    //------------------------------------------------------------//
    cout << "Start Measure" << endl;
    // Set the restart parameter to 1
    change_file("input_06.1/INPUT/input.dat","RESTART","1");
    change_file("input_06.1/INPUT/input.dat","NBLOCKS","500");
    change_file("input_06.1/INPUT/input.dat","NSTEPS","10000");
    // cycle over the magnetic field (on/off)
    for(int i=0; i<2; i++){
    //cycle over the simtype (metropolis/gibbs)
    for(int k=0; k<2; k++){
    //cycle over the temperature
    cout << "Metro=" << Metro[k] << " Field="<< Mfield[i] << endl;
    for (const string& str : temp_string){
        // Setting the input file's parameters
        if(Metro[k]){changeSimType("input_06.1/INPUT/input.dat","2\t1.0\t"+Mfield[i]);}else{changeSimType("input_06.1/INPUT/input.dat","3\t1.0\t"+Mfield[i]);} // set simulation type  // set simulation type
        //if(Mag_field[i]){change_Mag_field("input_06.1/INPUT/input.dat",0.02);}else{change_Mag_field("input_06.1/INPUT/input.dat",0.0);}        // set magnetic field
        change_file("input_06.1/INPUT/input.dat","TEMP",str);                                                                                  // set temperature
        // Create directories' paths
        string folderName ="Results/output_06.1/measure/"+met_gibb[k]+"/"+field[i]+"/temp_"+str;
        string folderName1 =folderName+"/OUTPUT";
        string folderName2 =folderName1+"/CONFIG";
        // Create directories
        try {
            if (fs::create_directory(folderName) and fs::create_directory(folderName1) and fs::create_directory(folderName2) ) {
                //cout << "Cartella creata con successo: " << folderName << endl;
                //cout << "Cartella creata con successo: " << folderName1 << endl;
                //cout << "Cartella creata con successo: " << folderName2 << endl;
            }else{
                //cout << "La cartella esiste già o non può essere creata: " << folderName << endl;
                //cout << "La cartella esiste già o non può essere creata: " << folderName1 << endl;
                //cout << "La cartella esiste già o non può essere creata: " << folderName2 << endl;
            }
        }catch (const fs::filesystem_error& e) {cerr << "Errore: " << e.what() << endl;}
        // Create and open file
        string input_path="input_06.1/";
        string output_path= folderName+"/";
        string config_path= "Results/output_06.1/equilibration/"+met_gibb[k]+"/"+field[i]+"/temp_"+str+"/";
        System SYS;
        SYS.initialize(input_path,output_path,config_path); 
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
    }
    }
    }

    //-----------------------------------------------------------------------------------------------//
    //  After simulation has been made we have to create some good files with all the data we need   //
    //-----------------------------------------------------------------------------------------------//
    
    // cycle over the magnetic field (on/off)
    for(int i=0; i<2; i++){
    //cycle over the simtype (metropolis/gibbs)
    for(int k=0; k<2; k++){
    string filepath="Results/output_06.1/";
    string filename="Results/output_06.1/results_"+met_gibb[k]+"_"+field[i]+".dat";
    ofstream ofile;
    ofile.open(filename);
    // Write the firs line of the file
    ofile << "Temperature;\t" << "Magnetization;\t" << "Error_Mag;\t" << "Specific_Heat;\t" << "Error_SpHe;\t" << "Susceptibility:\t" << "Error_Susc;\t" << "Total_Energy;\t" << "Error_tEn;\t" << endl;
    //cycle over the temperature
    for (const string& str : temp_string){
        // Get the measures: magnetization, specific heat, susceptibility, total energy
        vector<double> mag = getMeasure(filepath+"measure/"+met_gibb[k]+"/"+field[i]+"/temp_"+str+"/OUTPUT/magnetization.dat");
        vector<double> spec_heat = getMeasure(filepath+"measure/"+met_gibb[k]+"/"+field[i]+"/temp_"+str+"/OUTPUT/specific_heat.dat");
        vector<double> suscept = getMeasure(filepath+"measure/"+met_gibb[k]+"/"+field[i]+"/temp_"+str+"/OUTPUT/susceptibility.dat");
        vector<double> tener = getMeasure(filepath+"measure/"+met_gibb[k]+"/"+field[i]+"/temp_"+str+"/OUTPUT/total_energy.dat");
        // Write the data in the file 
        ofile << setw(16) << str << ";\t" << mag[0] << ";\t" << mag[1] << ";\t" << spec_heat[0] << ";\t" << spec_heat[1] << ";\t" << suscept[0] << ";\t" << suscept[1] << ";\t" << tener[0] << ";\t" << tener[1] << endl;
    }
    }
    }
    return 0;
}
