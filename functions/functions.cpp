#include "functions.hpp"
using namespace std;
using namespace arma;

double error(double AV,double AV2,int n){
    if(n==0) 
        return 0;
    else
        return sqrt((AV2 - AV * AV) / n);
} 

ofstream openofile(const string &filename){
   ofstream ofile(filename);
   if (!ofile.is_open()) {throw invalid_argument("Impossibile aprire il file per la scrittura.");}
   return ofile;  
}

void changeSimType(const std::string& filePath, const std::string& newSimulationType) {
    std::ifstream fileIn(filePath);  // Open the file in read mode
    if (!fileIn.is_open()) {
        std::cerr << "Error opening the file: " << filePath << std::endl;
        return;
    }
    std::string tempFilePath = filePath + ".tmp";  // Create a temporary file for writing
    std::ofstream fileOut(tempFilePath);
    std::string line;
    bool found = false;
    while (getline(fileIn, line)) {
        // Search for the line that contains "SIMULATION TYPE"
        size_t pos = line.find("SIMULATION_TYPE");
        if (pos != std::string::npos) {
            // Replace the entire content to the right of "SIMULATION TYPE"
            line = "SIMULATION_TYPE " + newSimulationType;
            cerr << line << endl;
            found = true;
        }
        fileOut << line << std::endl;  // Write the line (modified or not) into the temporary file
    }
    fileIn.close();
    fileOut.close();
    if (found) {
        // Replace the original file with the modified one
        if (std::remove(filePath.c_str()) != 0) {
            std::cerr << "Error deleting the original file" << std::endl;
            return;
        }
        if (std::rename(tempFilePath.c_str(), filePath.c_str()) != 0) {
            std::cerr << "Error renaming the temporary file" << std::endl;
            return;
        }
        //std::cout << "SIMULATION_TYPE updated successfully!" << std::endl;
    } else {
        std::cerr << "SIMULATION_TYPE not found in the file!" << std::endl;
        std::remove(tempFilePath.c_str());  // Remove the temporary file if nothing was found
    }
}


void change_file(const string &filename,string keyword,string new_value){
    ifstream ifile(filename);
    if (!ifile.is_open()) {throw invalid_argument("Impossibile aprire il file per la lettura.");}
    vector<string> lines;
    string line;
    while (getline(ifile, line)) {
        lines.push_back(line);
    }
    ifile.close();
    for (string& line : lines) {
        istringstream iss(line);
        string key;
        string value;
        if (iss >> key >> value) {
            if (key == keyword) {
                line = keyword + " " + new_value;
            }
        }
    }
    ofstream ofile(filename);
    if (!ofile.is_open()) {throw invalid_argument("Impossibile aprire il file per la scrittura.");}
    for (const string& line : lines) {
        ofile << line << endl;
    }
    ofile.close();
}

pair<double, double> calculateMeanAndStdDev(const double accumulate,const double sum_prog_helper, const double su2_prog_helper, const int dim, const int block){
    // calculation of the mean value and the standard deviation
    double sum_prog=sum_prog_helper/(block+1); // Cumulative average
    double su2_prog=su2_prog_helper/(block+1); // Cumulative square average
    double err_prog=error(sum_prog,su2_prog,block); // Statistical uncertainty
    if(err_prog==0 && block!=0) {throw invalid_argument("Errore nullo.");} // if the error is zero, the program stops
    return make_pair(sum_prog,err_prog);
}

vec discrete_walk(Random &rnd,const vec &p, const double a){
    vec p_new=p;
    double random_number=rnd.Rannyu(0,6);
    if(random_number>=0 && random_number<1){p_new[0]+=a;}
    if(random_number>=1 && random_number<2){p_new[0]-=a;}
    if(random_number>=2 && random_number<3){p_new[1]+=a;}
    if(random_number>=3 && random_number<4){p_new[1]-=a;}
    if(random_number>=4 && random_number<5){p_new[2]+=a;}
    if(random_number>=5 && random_number<6){p_new[2]-=a;}
    return p_new;
}

vec continuous_walk(Random &rnd, const vec &p, const double a){
    vec p_new=p;
    // uniform distribution around a sphere
    vec angles=rnd.spherical();
    double theta=angles[0];
    double phi=angles[1];
    p_new[0]+=a*sin(theta)*cos(phi);
    p_new[1]+=a*sin(theta)*sin(phi);
    p_new[2]+=a*cos(theta);
    return p_new;
} 

void print(const vec &p,ofstream &ofile){
    for(long unsigned int i=0;i<p.size()-1;i++){
        ofile << p[i] << ";";
    }
    ofile << p.back() << endl;
    if(!ofile.good()){throw invalid_argument("Errore nella scrittura del file.");}
}

double distance2(const vec &p){
    return p[0]*p[0]+p[1]*p[1]+p[2]*p[2];
}

void print(const vec &p){
    for(long unsigned int i=0;i<p.size()-1;i++){
        cout << p[i] << ";";
    }
    cout << p.back() << endl;
}

double spot_price(const double S0, const double mu, const double sigma, const double t, const double rnd){
    double wt=rnd;
    return S0*exp((mu-0.5*pow(sigma,2))*t+sigma*wt*sqrt(t));
}

double spot_price(const double S0, const double mu, const double sigma, const double t,Random &rnd,const int N_steps){
    double S=S0;
    double dt=t/N_steps;
    for(int i=0;i<N_steps;i++){
        double appoggio;
        appoggio=S*exp((mu-0.5*pow(sigma,2))*dt+sigma*rnd.Gauss(0,1)*sqrt(dt));
        S=appoggio;
    }
    return S;
}

double call_option_price(const double r, const double T, const double K, const double S){
    return exp(-r*T)*max(S-K,0.);
}

double put_option_price(const double r, const double T, const double K, const double S){
    return exp(-r*T)*max(K-S,0.);
}

void cartesianToSpherical(const vec& cartesian, vec& spherical) {
    int dimension = cartesian.n_elem;
    // Implement the conversion from Cartesian to Spherical coordinates for any dimension
    // This is a placeholder for the actual conversion logic.
    // You need to fill in the logic for converting Cartesian to Spherical coordinates
    // based on the dimension.
    
    // Example for 3D:
    if (dimension == 3) {
        double x = cartesian(0);
        double y = cartesian(1);
        double z = cartesian(2);
        double r = sqrt(x*x + y*y + z*z);
        double theta = acos(z / r);
        double phi = atan2(y, x);
        spherical = {r, theta, phi};
    } else {
        // Implement for other dimensions
        // Note: This will vary significantly based on the specific needs and mathematical
        // formula for higher dimensions.
        // Placeholder logic (incorrect for actual spherical coordinates in higher dims):
        spherical = cartesian;
    }
}

vector<double> getMeasure(const string& filename) {
    ifstream inFile(filename);
    string line;
    vector<string> lines;
    vector<double> result;
    if (!inFile) {
        cerr << "Error opening file: " << filename << endl;
        return result;
    }
    // Read all lines from the file and store them in a vector
    while (getline(inFile, line)) {
        lines.push_back(line);
    }
    inFile.close();
    // Check if the file has any lines
    if (lines.empty()) {
        cerr << "File is empty: " << filename << endl;
        return result;
    }
    // Get the last line
    string lastLine = lines.back();
    istringstream lineStream(lastLine);
    vector<double> numbers;
    double number;
    // Parse the last line to extract numbers
    while (lineStream >> number) {
        numbers.push_back(number);
    }
    // Check if there are at least 4 numbers in the last line
    if (numbers.size() < 4) {
        cerr << "The last line of the file " << filename << " does not contain enough numbers." << endl;
        return result;
    }
    // Retrieve the 3rd and 4th numbers and store them in the result vector
    result.push_back(numbers[2]);
    result.push_back(numbers[3]);
    return result;
}

vector<double> getAcceptance(const string& filename) {
    ifstream inFile(filename);
    string line;
    vector<string> lines;
    vector<double> result;
    if (!inFile) {
        cerr << "Error opening file: " << filename << endl;
        return result;
    }
    // Read all lines from the file and store them in a vector
    while (getline(inFile, line)) {
        lines.push_back(line);
    }
    inFile.close();
    // Check if the file has any lines
    if (lines.empty()) {
        cerr << "File is empty: " << filename << endl;
        return result;
    }
    // Get the last line
    string lastLine = lines.back();
    istringstream lineStream(lastLine);
    vector<double> numbers;
    double number;
    // Parse the last line to extract numbers
    while (lineStream >> number) {
        numbers.push_back(number);
    }
    // Check if there are at least 4 numbers in the last line
    if (numbers.size() < 2) {
        cerr << "The last line of the file " << filename << " does not contain enough numbers." << endl;
        return result;
    }
    // Retrieve the 3rd and 4th numbers and store them in the result vector
    result.push_back(numbers[2]);
    result.push_back(numbers[3]);
    return result;
}

void change_Mag_field(const string& filename, double newValue) {
    ifstream inFile(filename);
    ostringstream tempStream;
    string line;
    bool modified = false;
    if (!inFile) {
        cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    // Read and process each line
    while (getline(inFile, line)) {
        istringstream lineStream(line);
        string word;
        vector<string> tokens;
        // Split the line into tokens
        while (lineStream >> word) {
            tokens.push_back(word);
        }
        // Check if the line starts with "SIMULATION_TYPE" and has at least three numbers
        if (tokens.size() > 3 && tokens[0] == "SIMULATION_TYPE") {
            // Modify the third number
            tokens[3] = to_string(newValue);
            modified = true;
        }
        // Reassemble the line
        for (const auto& token : tokens) {
            tempStream << token << " ";
        }
        tempStream << endl;
    }
    inFile.close();
    if (!modified) {
        cerr << "No line starting with 'SIMULATION_TYPE' found or insufficient tokens." << std::endl;
        return;
    }
    // Write the modified content back to the file
    ofstream outFile(filename);
    if (!outFile) {
        cerr << "Error opening file for writing: " << filename << std::endl;
        return;
    }
    outFile << tempStream.str();
    outFile.close();
}

vector<string> StringVec(double start, double end, double increment) {
    vector<string> result;
    ostringstream stream;
    for (double value = start; value <= end; value += increment) {
        stream.str(""); // Clear the stringstream
        stream.clear(); // Clear any error flags
        stream << fixed << setprecision(2) << value;
        result.push_back(stream.str());
    }
    return result;
}

int pbc(int position, int size){
    // periodic boundary conditions
    return position - size * (int)floor((double)position/(double)size);
}

arma::mat get_coordinates(std::string filename) {
    arma::mat coordinates;
    // Open input file
    std::ifstream infile(filename);
    // verify that the file was opened correctly
    if (!infile.is_open()) {
        std::cerr << "Error: unable to open the file " << filename << std::endl;
        exit(-1);
    }
    // count the number of lines in the file
    int n = 0;
    std::string line;
    while (std::getline(infile, line)) {
        n++;
    }
    // reset the file to the beginning
    infile.clear();
    infile.seekg(0, std::ios::beg);
    // create a matrix to store the coordinates
    coordinates.set_size(n, 2);
    // read the coordinates from the file
    int i = 0;
    while (std::getline(infile, line)) {
        std::stringstream ss(line);
        double x, y;
        // read the coordinates
        if (ss >> x >> y) {
            coordinates(i, 0) = x;
            coordinates(i, 1) = y;
            i++;
        } else {
            std::cerr << "Error: bad row in the file: " << line << std::endl;
            exit(-1);
        }
    }
    // close the file
    infile.close();
    return coordinates;
}

void fillVector( ivec& son, const ivec& parent) {
    std::unordered_set<double> existingValues(son.begin(), son.end());
    size_t fillIndex = 1;
    for (size_t i=0;i < parent.size();i++) {
        int value = parent[i];
        if (existingValues.find(value) == existingValues.end()) {
            while (fillIndex < son.size() && son[fillIndex] != 0) {
                fillIndex++;
            }
            if (fillIndex < son.size()) {
                son[fillIndex] = value;
                existingValues.insert(value);
            }
        }
    }
}
