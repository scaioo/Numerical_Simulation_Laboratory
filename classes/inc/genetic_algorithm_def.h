#include <iostream>
#include <armadillo>
#include <unordered_set>
#include <map>
#include <unordered_map>
#include"../../random_number/random.h"
#include"../../functions/functions.hpp"
#include <utility>
#include <cmath>
using namespace std;
using namespace arma;

#ifndef INDIVIDUE_H
#define INDIVIDUE_H

class Individue{
public:
    // ----------------------------------------------------- //
    //                  CONSTRUCTORS                         //
    // ----------------------------------------------------- //
    // Default constructor
    Individue() : _dim(0), _rnd(nullptr) {}
    // Constructor
    Individue(int dim,Random rnd) : _dim(dim) {
        ivec indi(dim,arma::fill::zeros);
        _indi=indi;
        for(size_t i=0;i<_dim;i++){
            _indi(i)=i;
        }
    }
    // Copy constructor
    Individue& operator=(const Individue& other) {
        _dim = other._dim;
        _rnd = other._rnd;
        _indi=other._indi;
        _indi_map=other._indi_map;
        return *this;
        }
    
    //-------------------------------------------------------//
    //              INITIALIZATION METHODS                   //                    
    //-------------------------------------------------------//
    // initialization of individue on a circle
    void initialize_circle(int dim_indi,double r,Random &rnd) {
        _rnd = &rnd;
        _dim=dim_indi;
        _indi.set_size(_dim);
        for(size_t i=0;i<_dim;i++){_indi(i)=i;}
        double x, y,app;
        int sign;
        for(size_t i = 0; i < _dim; ++i) {
            x=_rnd->Rannyu(-r, r);
            double app = _rnd->Rannyu();
            sign = (app < 0.5) ? -1 : 1;
            y = sign * sqrt(r * r - x * x);
            _indi_map[i]={x,y};
        }
    }
    // initialization in a square
    void initialize_square(int dim_indi,double x,Random &rnd) {
        _rnd = &rnd;
        _dim=dim_indi;
        _indi.set_size(_dim);
        for(size_t i=0;i<_dim;i++){_indi(i)=i;}
        for(size_t i = 0; i < _dim; ++i) {
            _indi_map[i]={_rnd->Rannyu(-x, x),_rnd->Rannyu(-x, x)};
        }
    }
    // initialization of individue given the coordinates
    void initialize_city(int dim_indi,arma::mat coordinates,Random &rnd) {
        _rnd = &rnd;
        _dim=dim_indi;
        _indi.set_size(_dim);
        for(size_t i=0;i<_dim;i++){_indi(i)=i;}
        for(size_t i = 0; i < _dim; ++i) {
            _indi_map[i]={coordinates(i,0),coordinates(i,1)};
        }
    }
    
    // ----------------------------------------------------- //
    //                  GET METHODS                          //
    // ----------------------------------------------------- //
    // Get the positions of the individue
    mat Get_positions() const{
        mat positions(2,_dim);
        for(size_t i=0;i<_dim;i++){
            positions(0,i)=_indi_map.at(_indi(i))[0];
            positions(1,i)=_indi_map.at(_indi(i))[1];
        }
        return positions;
    }
    // Get the loss of the individue
    double get_indi_Loss() const {
        double loss = 0;
        for(size_t i=0;i<_dim-1;i++){
            loss+=arma::norm(_indi_map.at(_indi(i)) - _indi_map.at(_indi(i+1)),2);
        }
        loss+=arma::norm(_indi_map.at(_indi(_dim-1)) - _indi_map.at(_indi(0)),2);
        return loss;
    }
    // Get an index of the individue
    int get_idxs(int i) const {return _indi(i);}
    // Get the indexes of the individue
    ivec get_idxs() const {return _indi;}

    // ----------------------------------------------------- //
    //                 OTHER METHODS                         //
    // ----------------------------------------------------- //
    // shuffle the indexes of the individue
    ivec shuffleidxs() {
        ivec r_new = _indi.subvec(1, _indi.n_elem - 1);
        r_new = arma::shuffle(r_new);
        r_new.insert_rows(0,1);
        r_new(0)=_indi(0);
        _indi=r_new;
        return _indi;}
    // Set the indexes of the individue
    void set_idxs(ivec idxs) {
        if(idxs.n_elem != _dim) {
            cerr << "Error in set_idxs: invalid size of the indexes vector." << endl;
            exit(-1);
        }
        _indi = idxs;
    }

    void set_idxs(int index,int value) {
        if(index>_dim or index < 0){
            cerr << "Error in set_idxs: index bigger than dim indi or smaller than zero" << endl;
            exit(-1);
        }
        _indi(index)=value;

    }
    // Print the positions of the individue
    void print_positions() const {
        for (size_t i = 0; i < _dim; ++i) {
            cout << "Index: " << _indi[i] << ", Position: (" << _indi_map.at(i)[0] << ", " << _indi_map.at(i)[1] << ")" << "\t";
        }
    }
    // Print the indexes of the individue
    void print_idxs() const{
        cout << "Indexes: [";
        for(size_t i=0;i<_dim;i++){cout << _indi(i)<<" ";}
        cout << "]";
    }
    
    // check individue
    bool checkIndividual() {
        std::unordered_set<int> existingIndices;
        for (size_t i = 0; i < _dim; ++i) {
            if (existingIndices.find(_indi(i)) != existingIndices.end()) {
            cerr << "Error in checkIndividual: duplicate index found" << endl; exit(-1);
            }
            if(_indi(0)!=0){cerr << "Error in checkIndividual: first index is not 0" << endl; exit(-1);}
        }
        return true;
    }

    // ----------------------------------------------------- //
    //                  MUTATION METHODS                     //
    // ----------------------------------------------------- //
    // Pair permutation two indexes
    void pair_permutation(int i, int j){
        
        ivec r=_indi;

        if(i>=_dim-1 || j>=_dim-1){cerr << "Error in pair permutation : indexes out of range" << endl; exit(-1);}
        // copy the vector without the first element that remains fixed
        ivec r_new = r.subvec(1, r.n_elem - 1);
        ivec r_app=r_new;

        // swap the elements i and j
        r_new[i]=r_app[j];
        r_new[j]=r_app[i];

        r_new.insert_rows(0,1);
        r_new(0)=r(0);
        _indi=r_new;
    } 
    // Block shift of the indexes
    void block_shift(int k, int m, int n){

        // ------------------------- PARAMETERS ----------------------------
        // k start of the block, m length of the block, n length of the jump
        // -----------------------------------------------------------------
        ivec r=_indi;
        
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
        
        _indi=r_new;
}
    // Block permutation of the indexes
    void block_permutation(int m, int k, int n){
        
        // ------------------------- PARAMETERS ----------------------------
        // m length of the block, k start of the block1, n start of the block2
        // -----------------------------------------------------------------
    
        ivec r=_indi;

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
        _indi=r_new;
    }
    // Inversion of a block of the indexes
    void inversion(size_t startIdx, size_t endIdx) {
        ivec r=_indi;
        for(size_t i=0;i<_dim;i++){r(i)=i;}

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
        _indi=r_new;
    }

    // Destructor
    ~Individue() {}
    
private:
    map<int,arma::vec> _indi_map;
    ivec _indi;
    int _dim; 
    Random *_rnd;
};

#endif // INDIVIDUE_H


#ifndef POPULATION_H
#define POPULATION_H

class Population : Individue {
public:
    // Constructor
    Population(Random rnd) : Individue(), _dim_indi(0),_dim_pop(0),_rnd(&rnd) {}
    
    // ----------------------------------------------------- //
    //              INITIALIZATION METHODS                   //
    // ----------------------------------------------------- //
    // Initialization of the population on a circle
    void initialize_circle(int dim_indi, int dim_pop, double r,Random &rnd) {
        _dim_indi=dim_indi;
        _dim_pop=dim_pop;
        _rnd = &rnd;
        pop.set_size(_dim_pop);
        //first individue
        pop[0].initialize_circle(dim_indi,r,*_rnd);
        pop[0].shuffleidxs();
        for(size_t i = 1; i < _dim_pop; ++i) {
            pop[i] = pop[0];
            pop[i].shuffleidxs();
        }
    }
    // Initialization of the population on a square
    void initialize_square(int dim_indi, int dim_pop, double l,Random &rnd) {
        _dim_indi=dim_indi;
        _dim_pop=dim_pop;
        _rnd=&rnd;
        pop.set_size(_dim_pop);
        pop[0].initialize_square(dim_indi,l,*_rnd);
        for(size_t i = 0; i < _dim_pop; ++i) {
            pop[i] = pop[0];
            pop[i].shuffleidxs();
        }
    }
    // Initialization of the population given the coordinates
    void initialize_city(int dim_indi, int dim_pop, arma::mat &coordinates,Random &rnd) {
        _dim_indi=dim_indi;
        _dim_pop=dim_pop;
        _rnd=&rnd;
        pop.set_size(_dim_pop);
        pop[0].initialize_city(dim_indi,coordinates,*_rnd);
        for(size_t i = 0; i < _dim_pop; ++i) {
            pop[i] = pop[0];
            pop[i].shuffleidxs();
        }
    }
    
    // Operator []
    Individue& operator[](size_t index) {
        if (index >= pop.n_elem) {
            throw std::out_of_range("Index out of range");
        }
        return pop[index];
    }
    
    //-------------------------------------------------------//
    //                   GET METHODS                         //
    //-------------------------------------------------------//

    // Get the mean loss of the population
    double get_pop_mean_Loss(){
        // Sort the population based on the loss
        std::vector<std::pair<double, size_t>> loss_indices;
        for (size_t i = 0; i < pop.n_elem; ++i) {
            loss_indices.push_back(std::make_pair(pop[i].get_indi_Loss(), i));
        }
        std::sort(loss_indices.begin(), loss_indices.end());
        double loss=0;
        int half_pop=(int)_dim_pop/2;
        for(size_t i=0; i<half_pop;i++){
            loss+=loss_indices[i].first;
        }
        return loss/half_pop;
    }
    // Get the size of the population
    int get_pop_size() const { return _dim_pop;}
    // Get the size of the individue
    int get_indi_size() const { return _dim_indi;}       

    //-------------------------------------------------------//
    //                   PRINT METHODS                       //
    //-------------------------------------------------------//
    // Print the positions of the population
    void print_pop_positions() const {
        for (size_t i = 0; i < pop.n_rows; ++i) {
            cout << "Individue " << i << ":" << endl;
            pop[i].print_positions();
            cout << endl;
        }
    }
    // Print the indexes of the population
    void print_pop_idxs() {
        for (size_t i = 0; i < _dim_pop; ++i) {
            cout << "Individue " << i << ": ";
            pop[i].print_idxs();
            cout << endl;
        }
    }
    // Print the losses of the population
    void print_pop_loss(){
        for (size_t i = 0; i < pop.n_elem; ++i) {
            cout << "Individue " << i << " loss: " << pop[i].get_indi_Loss() << endl;
        }
    }
    // Print the full population
    void print_full_population() {
        for (size_t i = 0; i < _dim_pop; ++i) {
            cout << " Individue " << i << ": ";
            pop[i].print_idxs();
            cout << "Loss: " << pop[i].get_indi_Loss() << " ";
            cout << endl;
        }
    }
    
    // Sort the population based on the loss
    void sort(){ 
        arma::field<Individue> new_pop;
        new_pop.set_size(_dim_pop);
        // Sort the population based on the loss
        std::vector<std::pair<double, size_t>> loss_indices;
        for (size_t i = 0; i < pop.n_elem; ++i) {
            loss_indices.push_back(std::make_pair(pop[i].get_indi_Loss(), i));
        }
        std::sort(loss_indices.begin(), loss_indices.end());
        for(int i=0;i<_dim_pop;i++){
            new_pop[i] = pop[loss_indices[i].second];
        }
        pop=new_pop;
    }
    // Check if the population is valid
    void checkPopulation() {
        for (size_t i = 0; i < _dim_pop; ++i) {
            if (!pop[i].checkIndividual()) {
                pop[i].print_idxs();
                throw std::runtime_error("Invalid population member: individual "+std::to_string(i)+" is not valid:" );
            }
        }
    }

    // ----------------------------------------------------- //
    //                  MUTATION METHODS                     //
    // ----------------------------------------------------- //
    //Crossover
    std::pair<ivec ,ivec> crossover(ivec mother_idxs,ivec father_idxs, double rand){
        
        // ------------------------- PARAMETERS ----------------------------
        // mother and father are two vectors of the same size
        // -----------------------------------------------------------------

        ivec new_mother_idxs=mother_idxs;
        ivec new_father_idxs=father_idxs;

        if(mother_idxs.size()!=father_idxs.size()){cerr << "Error in crossover : vectors of different size: " << endl; exit(-1);}
        ivec child1(mother_idxs.n_elem,arma::fill::zeros);
        ivec child2(father_idxs.n_elem,arma::fill::zeros);  
        
        int bar = (int)(mother_idxs.n_elem * rand);
        for (int i = 0; i < bar; i++) {
            child1(i) = father_idxs(i);
            child2(i) = mother_idxs(i); 
        }        
        fillVector(child1, mother_idxs);
        fillVector(child2, father_idxs);


        new_mother_idxs=child1;
        new_father_idxs=child2;

        return std::make_pair(new_mother_idxs,new_father_idxs);
    }
    // Selection of the population
    void Selection() {
        int _sel;
        sort();
        arma::field<Individue> new_pop;
        new_pop.set_size(_dim_pop);
        for(int i=0;i<_dim_pop;i++){
            _sel=static_cast<int>((double)_dim_pop * (pow(_rnd->Rannyu(), 8)));
            new_pop[i] = pop[_sel];
        }
        // Replace the old population with the new one
        pop = new_pop;
        sort();
    }
    // Mutation of the population
    void Mutation(double rnd_cross,double rnd_pair,double rnd_shift,double rnd_blkperm,double rnd_inv){
        
        // variables declaration
        int m,n,k;
        double rand_cross,rand_mut;
        int index2=0;
        Individue sindi;
        std::pair<ivec,ivec> new_indi_idxs;
        
        for(int index1=0;index1<_dim_pop;index1++){
            
            // update indexes and individue
            index2=index1;
            sindi=pop[index1];
            
            // Crossover
            rand_cross=_rnd->Rannyu();
            if(rand_cross<rnd_cross){
                // Do not allow crossover between the same individue
                while(arma::approx_equal(pop[index1].get_idxs(), pop[index2].get_idxs(), "absdiff", 1e-10)){index2=(int)_rnd->Rannyu(0,_dim_pop);}
                double rand=_rnd->Rannyu();
                new_indi_idxs=crossover(pop[index1].get_idxs(),pop[index2].get_idxs(),rand);
                if(_rnd->Rannyu()<0.5){sindi.set_idxs(new_indi_idxs.first);}else{sindi.set_idxs(new_indi_idxs.second);}
            }
            
            // PERFORM MUTATION WITH PROBABILITY rnd_mut

            //Pair permutation
            rand_mut=_rnd->Rannyu();
            if(rand_mut<rnd_pair){
                sindi.pair_permutation((int)_rnd->Rannyu(0,_dim_indi-1),(int)_rnd->Rannyu(0,_dim_indi-1));
            }
            //Block shift
            rand_mut=_rnd->Rannyu();
            if(rand_mut<rnd_shift){
                k=(int)_rnd->Rannyu(0,_dim_indi-1);
                m=(int)_rnd->Rannyu(1,_dim_indi-1);
                n=(int)_rnd->Rannyu(1,_dim_indi-1-m);
                if(_rnd->Rannyu()<0.5){sindi.block_shift(k,m,n);}else{sindi.block_shift(k,n,m);}
            }
            //Block permutation
            rand_mut=_rnd->Rannyu();
            if(rand_mut<rnd_blkperm){
                m=(int)_rnd->Rannyu(0,(int)(_dim_indi)/2);
                k=(int)_rnd->Rannyu(0,_dim_indi-1-m);
                n=(int)_rnd->Rannyu(k+m+1,_dim_indi+k-m-1);
                sindi.block_permutation(m,k,n);
            }
            //Inversion
            rand_mut=_rnd->Rannyu();
            if(rand_mut<rnd_inv){               
                m=(int)_rnd->Rannyu(0,_dim_indi-1);
                k=(int)_rnd->Rannyu(0,_dim_indi-1);
                sindi.inversion(k,m);
            }
            pop[index1]=sindi;
        }
    }

    // Destructor
    ~Population() {}
private:
    arma::field<Individue>  pop;
    int _dim_indi, _dim_pop, _sel; Random *_rnd;
};

#endif // POPULATION_H 
