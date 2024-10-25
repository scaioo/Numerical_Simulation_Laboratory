#include <iostream>
#include <armadillo>
#include <unordered_set>
#include <unordered_map>
#include"../../random_number/random.h"
#include"../../funzioni/functions.hpp"
#include <utility>
#include <cmath>
using namespace std;
using namespace arma;

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

#ifndef GENE_H
#define GENE_H

class gene {
public:
    // Default constructor
    gene() : _index(0), _x(0), _y(0), _rnd(nullptr) {}
    
    // Constructor that initializes based on a radius
    void initialize_circle(int index, double r, Random *rnd){ 
        _index=index;
        _rnd=rnd;
        _x = _rnd->Rannyu(-r, r);
        double app = _rnd->Rannyu();
        int sign = (app < 0.5) ? -1 : 1;
        _y = sign * sqrt(r * r - _x * _x);
    }
    // Constructor that initializes based on a square
    void initialize_square(int index, double l, Random *rnd){
        _index=index;
        _rnd=rnd;
        _x = _rnd->Rannyu(-l, l);
        _y = _rnd->Rannyu(-l, l);
    }
    // Constructor that initializes cities
    void initialize_city(int index, double x, double y){
        _index=index;
        _x=x;
        _y=y;
    }
    // Copy constructor
    gene(const gene& other) : _index(other._index), _x(other._x), _y(other._y), _rnd(other._rnd) {}

    // Operators
    gene& operator=(const gene& other) {
        if (this != &other) {
            _index = other._index;
            _x = other._x;
            _y = other._y;
            _rnd = other._rnd;
        }
    return *this;
    }
    
    double operator[](size_t index) const {
        if (index == 0) {
            return (double)_index;
        } else if (index == 1) {
            return _x;
        } else if (index == 2) {
            return _y;
        } else {
            throw std::out_of_range("Index out of range");
        }
    }

    // Set methods
    void SetIndex(int index) { _index = index; }
    
    void SetPosition(double x, double y) { _x = x; _y = y; }
    
    // Get methods
    int GetIndex() const { return _index; }
    
    vec GetPosition() const { return vec{_x, _y}; }
    
    // Print method for debugging
    void print() const {
        cout << "Index: " << _index << ", Position: (" << _x << ", " << _y << ")" << std::endl;
    }

    // Destructor
    ~gene() {}

private:
    int _index;
    double _x, _y;
    Random *_rnd;  // Store a pointer to the Random object
};

#endif // GENE_H

#ifndef INDIVIDUE_H
#define INDIVIDUE_H

class Individue : gene{
public:
    // Default constructor
    Individue() : indi(), _rnd(nullptr) {}
    
    // initialization of individue on a circle
    void initialize_circle(int dim,double r, Random *rnd) {
        _dim=dim;
        _rnd=rnd;
        indi.set_size(_dim, 1);
        for(size_t i = 0; i < _dim; ++i) {
            indi(i, 0).initialize_circle(i, r, _rnd);
        }
    }
    // initialization in a square
    void initialize_square(int dim,double x, Random *rnd) {
        _dim=dim;
        _rnd=rnd;
        indi.set_size(_dim, 1);
        for(size_t i = 0; i < _dim; ++i) {
            indi(i, 0).initialize_square(i, x, _rnd);
        }
    }
    // initialization of the individue with the cities
    void initialize_city(int dim, const mat& coordinates) {
        _dim=dim;
        indi.set_size(_dim, 1);
        for(size_t i = 0; i < _dim; ++i) {
            indi(i, 0).initialize_city(i, coordinates(i, 0), coordinates(i, 1));
        }
    }

    // Operators
    gene& operator[](size_t index) {
        if (index >= indi.n_elem) {
            throw std::out_of_range("Index out of range");
        }
        return indi(index);
    }
    
    Individue& operator=(const Individue& other) {
        if (this != &other) {
            _dim = other._dim;
            _rnd = other._rnd;
            indi = other.indi;
        }
        return *this;
    }

    // Get methods
    ivec Getidxs() const {
        ivec idxs(_dim);
        for (size_t i = 0; i < _dim; ++i) {
            idxs(i) = indi(i, 0).GetIndex();
        }
    return idxs;
    }
    
    vec Get_Position(int i) const {
        return indi(i, 0).GetPosition();
    } 

    mat Get_positions() const{
        int index=0;
        mat positions(2,_dim);
        while(index!=_dim){
            for(size_t i=0;i<_dim;i++){
            if(index==indi(i,0).GetIndex()){
                positions.col(index) = indi(i,0).GetPosition();
                index++;
                break;
                }
            } // Assign the first element of the vec to the positions vector
        }
        return positions; // Add return statement
    }

    double get_indi_Loss() {
        double loss = 0;
        std::unordered_map<int, arma::vec> index_to_position;
        int min_index = std::numeric_limits<int>::max();
        int max_index = std::numeric_limits<int>::min();

        // Build the map from index to position and find min and max indices
        for (size_t i = 0; i < _dim; i++) {
            int index = indi(i, 0).GetIndex();
            arma::vec position = indi(i, 0).GetPosition();
            index_to_position[index] = position;

            if (index < min_index) {
                min_index = index;
            }
            if (index > max_index) {
                max_index = index;
            }
        }

        // Iterate through the indices and calculate the loss
        for (size_t i = 0; i < _dim; i++) {
            int index = indi(i, 0).GetIndex();
            arma::vec position = indi(i, 0).GetPosition();

            // Look for the next index in the map
            auto it = index_to_position.find(index + 1);
            if (it != index_to_position.end()) {
                // Calculate the norm and update the loss
                loss += arma::norm(position - it->second, 2);
            }
        }

        // Add the distance between the final point and the first point
        if (index_to_position.find(min_index) != index_to_position.end() && 
            index_to_position.find(max_index) != index_to_position.end()) {
            arma::vec start_position = index_to_position[min_index];
            arma::vec end_position = index_to_position[max_index];
            loss += arma::norm(end_position - start_position, 2);
        }

        return loss;
    }

    // Indexes operations
    void ChangeIndexes(const ivec& new_idxs) {
        if (new_idxs.n_elem != _dim) {
            throw std::invalid_argument("Invalid size for new indexes");
        }
        for (size_t i = 0; i < _dim; ++i) {
            indi[i].SetIndex(new_idxs(i));
        }
    }

    void shuffleidxs() {
        ivec idxs = Getidxs();
        double first_element = idxs(0);
        arma::ivec sub_vec = idxs.subvec(1, idxs.n_elem - 1);
        sub_vec = arma::shuffle(sub_vec);
        arma::ivec shuffled_vec(idxs.n_elem);
        shuffled_vec(0) = first_element;
        shuffled_vec.subvec(1, idxs.n_elem - 1) = sub_vec;
        ChangeIndexes(shuffled_vec);
    }
    
// print methods
    void print() const {
        for (size_t i = 0; i < indi.n_rows; ++i) {
            indi(i, 0).print();
        }
    }

    void printIdxs() const {
        for (size_t i = 0; i < indi.n_rows; ++i) {
            cout << indi(i, 0).GetIndex() << " ";
        }
    }

    void printPositions() const {
        for (size_t i = 0; i < indi.n_rows; ++i) {
            cout << "Index: " << indi(i, 0).GetIndex() << ", Position: (" << indi(i, 0)[1] << ", " << indi(i, 0)[2] << ")" << std::endl;
        }
    }

    void printIdxs(ofstream &file) const {
        for (size_t i = 0; i < indi.n_rows; ++i) {
            file << indi(i, 0).GetIndex() << " ";
        }
    }

    void print_positions() {
        for (size_t i = 0; i < _dim; ++i) {
            cout << " (" << indi(i, 0)[1] << ", " << indi(i, 0)[2] << ") ";
        }
    }
    
    void print_positions(ofstream &file) {
        for (size_t i = 0; i < _dim; ++i) {
            file << " (" << indi(i, 0)[1] << ", " << indi(i, 0)[2] << ") ";
        }
    }
    
    //----------------------------------------------------------------------------------------------------//
    // Functions to change the representation of the indexes between two different representations: 
    // a representation where the indexes are the step of the trip
    // a representation where the indexes are the position of the cities
    //----------------------------------------------------------------------------------------------------//
    ivec change_rep(ivec &vec){
        // key has the two different representations
        std::vector<std::pair<int, int>> key;
        for (size_t i = 0; i < vec.n_elem; ++i) {
            key.push_back(std::make_pair(vec(i), i));
        }
        // sort over the first element of the pair
        std::sort(key.begin(), key.end());
        // create the new vector
        ivec result(vec.n_elem);
        for (size_t i = 0; i < key.size(); ++i) {
            result(i) = key[i].second;
        }
        return result;
    }
 
    // ----------------------------------------------------- //
    //                 SERIALIZE INDIVIDUE                   //
    // ----------------------------------------------------- //

    void serialize_individue(const Individue &ind, std::vector<double> &buffer) {
        buffer.clear();
        // Serialize each gene (index, x, y)
        for (size_t i=0;i<indi.n_elem;i++) {
            buffer.push_back(static_cast<double>(indi(i, 0).GetIndex())); // Cast index to double
            buffer.push_back(indi(i, 0).GetPosition()(0));
            buffer.push_back(indi(i, 0).GetPosition()(1));
        }
    }
    
    Individue deserialize_individue(const std::vector<double> &buffer, int num_genes,const Individue &indi) {
        Individue ind=indi;
        int buffer_index=0;
        for (int i = 0; i < num_genes; ++i) {
            int idx = static_cast<int>(buffer[buffer_index]);  // Index (cast from double)
            double x = buffer[1 + buffer_index];                   // x coordinate
            double y = buffer[2 + buffer_index];                   // y coordinate
            ind[i].SetIndex(idx);
            ind[i].SetPosition(x, y);
            buffer_index+=3;
        }
        return ind;
    }

    // Destructor
    ~Individue() {}

// ----------------------------------------------------- //
    //                  MUTATION METHODS                     //
    // ----------------------------------------------------- //

    void pair_permutation(int i, int j){
        
        //if(i==j){cerr << "Error in pair permutation : indexes are the same" << endl; exit(-1);}
        if(i>=_dim-1 || j>=_dim-1){cerr << "Error in pair permutation : indexes out of range" << endl; exit(-1);}
        ivec r = Getidxs();
        // change the representation
        r=change_rep(r);

        // copy the vector without the first element that remains fixed
        ivec r_new = r.subvec(1, r.n_elem - 1);
        ivec r_app=r_new;

        // swap the elements i and j
        r_new(i)=r_app(j);
        r_new(j)=r_app(i);

        // add the first element
        r_new.insert_rows(0,1);
        r_new(0)=r(0);
        // change the representation back
        r_new=change_rep(r_new);
        ChangeIndexes(r_new);
    } 

    void block_shift(int k, int m, int n){

        // ------------------------- PARAMETERS ----------------------------
        // k start of the block, m length of the block, n length of the jump
        // -----------------------------------------------------------------
        ivec r=Getidxs();
        // change the representation
        r=change_rep(r);        
        // check if the block length is greater than the vector size
        if(m+n>r.size()-1){cerr << "Error in Block Shift function : block length and jump big: m="<< m << " n= "<< n << endl; exit(-1);}

        // copy the vector without the first element that remains fixed
        ivec r_new = r.subvec(1, r.n_elem - 1);
        ivec r_app=r_new;

        // move the block of length m of n positions
        for(int i=k;i<k+m;i++){
            r_new(pbc(i+n,r_new.size()))=r_app(pbc(i,r_new.size()));
            } // move the block of length m of n positions 
        for(int i=k+m;i<k+m+n;i++){r_new(pbc(i-m,r_new.size()))=r_app(pbc(i,r_new.size()));} // complete the remaining part of the vector

        // add the first element
        r_new.insert_rows(0,1);
        r_new(0)=r(0);
        // change the representation back
        r_new=change_rep(r_new);        
        ChangeIndexes(r_new);
}

    void block_permutation(int m, int k, int n){
    
    // ------------------------- PARAMETERS ----------------------------
    // m length of the block, k start of the block1, n start of the block2
    // -----------------------------------------------------------------
  
    ivec r=Getidxs();
    // change the representation
    r=change_rep(r);    

    if(((r.size())/2)<(m+1)){cerr << "Error in block permutation : block length or start index too big " << endl; exit(-1);}
    if(fabs(n-k)<(m+1)){cerr << "Error in block permutation : block length or start index too big " << endl; exit(-1);}
    // copy the vector without the first element that remains fixed
    ivec r_new = r.subvec(1, r.n_elem - 1);
    ivec r_app=r_new;
    ivec blk1,blk2;
    for(int i=0;i<m;i++){
        blk1.insert_rows(i,1);
        blk1(i)=r_app(pbc(k+i,r_new.size()));
        blk2.insert_rows(i,1);
        blk2(i)=r_app(pbc(n+i,r_new.size()));
        }
    for(int i=0;i<m;i++){
        r_new(pbc(k+i,r_new.size()))=blk2(i);
        r_new(pbc(n+i,r_new.size()))=blk1(i);
        }
    // add the first element
    r_new.insert_rows(0,1);
    r_new(0)=r(0);
    // change the representation back
    r_new=change_rep(r_new);
    ChangeIndexes(r_new);
}

    void inversion(size_t startIdx, size_t endIdx) {
    ivec r=Getidxs();
    //change the representation
    r=change_rep(r);

     // Ensure the parameters are within valid ranges
     if (startIdx > r.n_elem || endIdx > r.n_elem) {
         cerr << "Error in block inversion: invalid parameters." << endl;
         exit(-1);
         }
     if(startIdx > endIdx) {size_t temp = startIdx; startIdx = endIdx; endIdx = temp;}
     // Copy the original vector
     ivec r_new = r.subvec(1, r.n_elem - 1); 
     // Reverse the subvector between startIdx and endIdx
     while (startIdx < endIdx) {
         std::swap(r_new(startIdx), r_new(endIdx));
         ++startIdx;
         --endIdx;
     }
     // add the first element
     r_new.insert_rows(0,1);
     r_new(0)=r(0);
    
    // change the representation back
    r_new=change_rep(r_new);
    
    ChangeIndexes(r_new);
 }   

    // Check if the individual is valid
    bool checkIndividual() {
    ivec vec = Getidxs();
    if (vec.n_elem == 0) {
        return false; // Empty vector case
    }
    // Check if the first element is zero
    if (vec(0) != 0) {
        return false;
    }
    // Check if all elements are unique
    std::unordered_set<double> uniqueElements;
    for (size_t i = 0; i < vec.n_elem; ++i) {
        if (uniqueElements.find(vec(i)) != uniqueElements.end()) {
            return false; // Found a duplicate
        }
        uniqueElements.insert(vec(i));
    }
    return true;
}

private:
    arma::field<gene> indi;
    int _dim; Random *_rnd;
};

#endif // INDIVIDUE_H

#ifndef POPULATION_H
#define POPULATION_H

class Population : Individue {
public:
    Population() : pop(), _rnd(nullptr) {}

    // Initialization of the population on a circle
    void initialize_circle(int dim, int n, double r, Random *rnd) {
        _dim_indi=dim;
        _dim_pop=n;
        _rnd=rnd;
        pop.set_size(_dim_pop, 1);
        //first individue
        pop(0, 0).initialize_circle(_dim_indi, r, _rnd);
        pop(0,0).shuffleidxs();
        for(size_t i = 1; i < _dim_pop; ++i) {
            pop(i, 0) = pop(0, 0);
            pop(i,0).shuffleidxs();
        }
    }
    // Initialization of the population on a square
    void initialize_square(int dim, int n, double x, Random *rnd) {
        _dim_indi=dim;
        _dim_pop=n;
        _rnd=rnd;
        pop.set_size(_dim_pop, 1);
        pop(0, 0).initialize_square(_dim_indi, x, _rnd);
        for(size_t i = 0; i < _dim_pop; ++i) {
            pop(i, 0) = pop(0, 0);
            pop(i,0).shuffleidxs();
        }
    }
    // Initialization of the population with the cities
    void initialize_city(int dim, int n, const mat& coordinates, Random *rnd) {
        _dim_indi=dim;
        _dim_pop=n;
        _rnd=rnd;
        pop.set_size(_dim_pop, 1);
        for(size_t i = 0; i < _dim_pop; ++i) {
            pop(i, 0).initialize_city(_dim_indi, coordinates);
            pop(i,0).shuffleidxs();
        }
    }

    // Operators
    Individue& operator[](size_t index) {
        if (index >= pop.n_elem) {
            throw std::out_of_range("Index out of range");
        }
        return pop(index);
    }
    
    // Get methods
    double get_pop_mean_Loss(){

        // Sort the population based on the loss
        std::vector<std::pair<double, size_t>> loss_indices;
        for (size_t i = 0; i < pop.n_elem; ++i) {
            loss_indices.push_back(std::make_pair(pop(i, 0).get_indi_Loss(), i));
        }
        std::sort(loss_indices.begin(), loss_indices.end());
        double loss=0;
        int half_pop=(int)_dim_pop/2;
        for(size_t i=0; i<half_pop;i++){
            loss+=loss_indices[i].first;
        }
        return loss/half_pop;
    }
    
    int get_pop_size() const { return _dim_pop;}
    
    int get_indi_size() const { return _dim_indi;}       

    // Print methods
    void print() const {
        for (size_t i = 0; i < pop.n_rows; ++i) {
            cout << "Individue " << i << ":" << endl;
            pop(i, 0).print();
            cout << endl;
        }
    }
    
    void printIdxs() const {
        for (size_t i = 0; i < pop.n_rows; ++i) {
            cout << "Individue " << i << ": ";
            pop(i, 0).printIdxs();
            cout << endl;
        }
    }
     
    void print_population() {
        for (size_t i = 0; i < _dim_pop; ++i) {
            cout << "Individue " << i << ": ";
            pop(i, 0).printIdxs();
            cout << endl;
        }
    }

    // Print the losses of the population
    void print_pop_loss(){
        for (size_t i = 0; i < pop.n_elem; ++i) {
            cout << "Individue " << i << " loss: " << pop(i, 0).get_indi_Loss() << endl;
        }
    }

    void print_full_population(ofstream &file) {
        for (size_t i = 0; i < _dim_pop; ++i) {
            file << " Individue " << i << ": ";
            pop(i, 0).printIdxs(file);
            file << "Loss: " << pop(i, 0).get_indi_Loss() << " ";
            pop(i,0).print_positions(file);
            file << endl;
        }
    }  
    // Sort the population based on the loss
    void sort(){ 
        arma::field<Individue> new_pop;
        new_pop.set_size(_dim_pop,1);
        // Sort the population based on the loss
        std::vector<std::pair<double, size_t>> loss_indices;
        for (size_t i = 0; i < pop.n_elem; ++i) {
            loss_indices.push_back(std::make_pair(pop(i, 0).get_indi_Loss(), i));
        }
        std::sort(loss_indices.begin(), loss_indices.end());
        for(int i=0;i<_dim_pop;i++){
            new_pop(i,0) = pop(loss_indices[i].second,0);
        }
        pop=new_pop;
    }
    // Check if the population is valid
    void checkPopulation() {
        for (size_t i = 0; i < _dim_pop; ++i) {
            if (!pop(i).checkIndividual()) {
                pop(i).printIdxs();
                throw std::runtime_error("Invalid population member: individual "+std::to_string(i)+" is not valid:" );
            }
        }
    }

    // Destructor
    ~Population() {}

    // ----------------------------------------------------- //
    //                  MUTATION METHODS                     //
    // ----------------------------------------------------- //

    //Crossover
    std::pair<Individue,Individue> crossover(const Individue &indi_mother,const Individue &indi_father, double rand){
        
        // ------------------------- PARAMETERS ----------------------------
        // mother and father are two vectors of the same size
        // -----------------------------------------------------------------

        ivec mother=indi_mother.Getidxs();
        ivec father=indi_father.Getidxs();
        //change the representations
        mother=change_rep(mother);
        father=change_rep(father);

        Individue new_mother=indi_mother;
        Individue new_father=indi_father;

        if(mother.size()!=father.size()){cerr << "Error in crossover : vectors of different size: " << endl; exit(-1);}
        ivec child1(mother.n_elem,arma::fill::zeros);
        ivec child2(father.n_elem,arma::fill::zeros);
        //child1.print("child1");
        //child2.print("child2");  
        
        int bar = (int)(mother.n_elem * rand);
        for (int i = 0; i < bar; i++) {
            child1(i) = father(i);
            child2(i) = mother(i); 
        }
        //child1.print("child1");
        //child2.print("child2");        
        fillVector(child1, mother);
        fillVector(child2, father);
        //child1.print("child1");
        //child2.print("child2");
        
        //change the representation back
        child1=change_rep(child1);
        child2=change_rep(child2);

        new_mother.ChangeIndexes(child1);
        new_father.ChangeIndexes(child2);

        return std::make_pair(new_mother,new_father);
    }
    // Selection of the population
    void Selection() {
        int _sel;
        sort();
        arma::field<Individue> new_pop;
        new_pop.set_size(_dim_pop, 1);
        for(int i=0;i<_dim_pop;i++){
            _sel=static_cast<int>((double)_dim_pop * (pow(_rnd->Rannyu(), 10)));
            new_pop(i,0) = pop(_sel,0);
        }
        // Replace the old population with the new one
        pop = new_pop;
        sort();
    }
    // Mutation of the population
    void Mutation(double rnd_cross,double rnd_pair,double rnd_shift,double rnd_blkperm,double rnd_inv){
        int m,n,k;
        double rand_cross,rand_mut;
        int index2=0;
        Individue sindi;
        std::pair<Individue,Individue> new_indi;
        // Perform crossover with a probability
        for(int index1=0;index1<_dim_pop;index1++){
            rand_cross=_rnd->Rannyu();
            index2=index1;
            sindi=pop(index1,0);
            if(rand_cross<rnd_cross){
                //cerr << "crossover" << endl;
                //non far fare crossover tra due individui uguali
                while(arma::approx_equal(pop(index1,0).Getidxs(), pop(index2,0).Getidxs(), "absdiff", 1e-10)){index2=(int)_rnd->Rannyu(0,_dim_pop);}
                double rand=_rnd->Rannyu();
                //cerr << "rand: " << rand << endl;
                //cerr << "old indi:"<< endl;
                //pop(index1,0).printIdxs();
                //pop(index2,0).printIdxs();

                new_indi=crossover(pop(index1,0),pop(index2,0),rand);
                //cerr << "new indi:"<< endl;
                //new_indi.first.printIdxs();
                //new_indi.second.printIdxs();
                if(_rnd->Rannyu()<0.5){sindi=new_indi.first;}else{sindi=new_indi.second;}
            }
            // Perform mutation with a probability
            rand_mut=_rnd->Rannyu();
            if(rand_mut<rnd_pair){
                //cerr << "pair permutation" << endl;
                sindi.pair_permutation((int)_rnd->Rannyu(0,_dim_indi-1),(int)_rnd->Rannyu(0,_dim_indi-1));
                }
            rand_mut=_rnd->Rannyu();
            if(rand_mut<rnd_shift){
                //cerr << "block shift" << endl;
                k=(int)_rnd->Rannyu(0,_dim_indi-1);
                m=(int)_rnd->Rannyu(1,_dim_indi-1);
                n=(int)_rnd->Rannyu(1,_dim_indi-1-m);
                if(_rnd->Rannyu()<0.5){sindi.block_shift(k,m,n);}else{sindi.block_shift(k,n,m);}
                }
            rand_mut=_rnd->Rannyu();
            if(rand_mut<rnd_blkperm){
                //cerr << "block permutation" << endl;
                m=(int)_rnd->Rannyu(0,(int)(_dim_indi)/2);
                k=(int)_rnd->Rannyu(0,_dim_indi-1-m);
                n=(int)_rnd->Rannyu(k+m+1,_dim_indi+k-m-1);
                sindi.block_permutation(m,k,n);
                }
            rand_mut=_rnd->Rannyu();
            if(rand_mut<rnd_inv){       
                //cerr << "inversion" << endl;         
                m=(int)_rnd->Rannyu(0,_dim_indi-1);
                k=(int)_rnd->Rannyu(0,_dim_indi-1);
                sindi.inversion(k,m);
                }
            pop(index1,0)=sindi;
            }
        }
   
    // ----------------------------------------------------- //
    //                 SERIALIZE POPULATION                  //
    // ----------------------------------------------------- //

    void serialize_population(vector<double> &buffer) {
        vector<double> buffer_ind;
        // Serialize each individual
        for (size_t i=0;i<pop.n_elem*0.25;i++) {
            pop(i, 0).serialize_individue(pop(i, 0), buffer_ind);
            // print individue serialized
            buffer.insert(buffer.end(),buffer_ind.begin(),buffer_ind.end());
        }
    }

    Population deserialize_population(const vector<double> &buffer, int num_genes,const Population &popu) {
        // Deserialize each individual
        Population pop=popu; 
        size_t start = 0; 
        for (size_t i = 0; i < pop.get_pop_size(); i++) {
            vector<double> buffer_ind(buffer.begin() + start, buffer.begin() + start + num_genes * 3);
            pop[i] = pop[i].deserialize_individue(buffer_ind, num_genes,pop[i]);
            start += num_genes * 3;
        }
        return pop;
    }

private:
    arma::field<Individue> pop;
    int _dim_indi, _dim_pop,_sel; Random *_rnd;
};

#endif // POPULATION_H
