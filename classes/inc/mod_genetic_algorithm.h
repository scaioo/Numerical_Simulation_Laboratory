#include <iostream>
#include <armadillo>
#include <unordered_set>
#include <map>
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
    gene(int n_genes) : _n_genes(n_genes),_x(0), _y(0), _rnd(nullptr) {generateMap();}
    
    void generateMap(){
        for(size_t i=0;i<_n_genes;i++){
            arma::vec vec={0,0};
            gene_map[i]=vec;
        }
    }
    // Constructor that initializes based on a radius
    void initialize_circle(int index, double r, Random &rnd){ 
        *_rnd = rnd;
        arma:: vec pos_gene;
        for(size_t i=0;i<_n_genes;i++){
            _x = _rnd->Rannyu(-r, r);
            double app = _rnd->Rannyu();
            int sign = (app < 0.5) ? -1 : 1;
            _y = sign * sqrt(r * r - _x * _x);
            gene_map[i]={_x,_y}
        }
        return vec{_x, _y};
    }
    // Constructor that initializes based on a square
    vec initialize_square(int index, double l, Random &rnd) {
        *_rnd = rnd;
        _x = _rnd->Rannyu(-l, l);
        _y = _rnd->Rannyu(-l, l);
        return vec{_x, _y};
    }
    // Destructor
    ~gene() {}

private:
    double _x, _y;
    Random *_rnd;  // Store a pointer to the Random object
    std::map<int, arma::ivec> gene_map;  // Map with int key and position
    int _n_genes;  // Number of genes
};

#endif // GENE_H

#ifndef INDIVIDUE_H
#define INDIVIDUE_H

class Individue : gene{
public:
    // Default constructor
    Individue() : gene(),_dim(0), _rnd(nullptr),_idxs({0,0}) {}
    // initialization of individue on a circle
    void set_map() {
        for(size_t i = 0; i < _dim; ++i) {
            gene new_gene;
            _indi.insert_or_assign(i, new_gene.initialize_circle(i, 1., *_rnd));
            _idxs(i)=i;
        }
    }
    void initialize_circle(int dim,double r,Random &rnd) {
        _dim=dim;
        _rnd = &rnd;
        for(size_t i = 0; i < _dim; ++i) {
            gene new_gene;
            _indi[i] = new_gene.initialize_circle(i, r, *_rnd);
            cerr << _indi[i][0] << " " << _indi[i][1] << endl;
            _idxs(i)=i;
        }
    }
    // initialization in a square
    void initialize_square(int dim,double x,Random &rnd) {
        _dim=dim;
        _rnd = &rnd;
        for(size_t i = 0; i < _dim; ++i) {
            gene new_gene;
            _indi[i]=new_gene.initialize_square(i, x, *_rnd);
            cerr << _indi[i][0] << " " << _indi[i][1] << endl;
            _idxs(i)=i;
    
        }
    }

    Individue& operator=(const Individue& other) {
            if (this != &other) {
                _dim = other._dim;
                _rnd = other._rnd;
                _idxs=other._idxs;
                for(size_t i = 0; i < _dim; ++i) {
                    _indi.insert_or_assign(i, other._indi.at(i));
            }
            }
            return *this;
    
        }

    mat Get_positions() const{
        mat positions(2,_dim);
        for(size_t i=0;i<_dim;i++){
            cerr << _indi.at(_idxs(i))[0] << " " << _indi.at(_idxs(i))[1] << endl;
            positions(0,i)=_indi.at(_idxs(i))[0];
            positions(1,i)=_indi.at(_idxs(i))[1];
        }
        return positions;
    }

    double get_indi_Loss() const {
        double loss = 0;
        for(size_t i=0;i<_dim-1;i++){
            loss+=arma::norm(_indi.at(_idxs(i)) - _indi.at(_idxs(i+1)),2);
        }
        loss+=arma::norm(_indi.at(_idxs(_dim-1)) - _indi.at(_idxs(0)),2);
        return loss;
    }

    ivec get_idxs() const {return _idxs;}

    ivec shuffleidxs() {
        return arma::shuffle(_idxs);
    }

    void set_idxs(ivec idxs) {
        _idxs = idxs;
    }

    void print_positions() const {
        for (size_t i = 0; i < _dim; ++i) {
            cout << "Index: " << _idxs(i) << ", Position: (" << _indi[i][0] << ", " << _indi[i][1] << ")" << std::endl;
        }
    }
    
    void print_idxs() const{
        cout << "Indexes: [";
        for(size_t i=0;i<_dim;i++){cout << _idxs(i)<<" ";}
        cout << "]";
    }
    
    // check individue
    bool checkIndividual() {
        std::unordered_set<int> existingIndices;
        for (size_t i = 0; i < _dim; ++i) {
            if (existingIndices.find(_idxs(i)) != existingIndices.end()) {
            cerr << "Error in checkIndividual: duplicate index found" << endl; exit(-1);
            }
            if(_idxs(0)!=0){cerr << "Error in checkIndividual: first index is not 0" << endl; exit(-1);}
        }
        return true;
    }
    
    // Destructor
    ~Individue() {}

    // ----------------------------------------------------- //
    //                  MUTATION METHODS                     //
    // ----------------------------------------------------- //

    void pair_permutation(int i, int j){
        
        ivec r=_idxs;

        //if(i==j){cerr << "Error in pair permutation : indexes are the same" << endl; exit(-1);}
        if(i>=_dim-1 || j>=_dim-1 || i==0 || i==0 ){cerr << "Error in pair permutation : indexes out of range" << endl; exit(-1);}
        // copy the vector without the first element that remains fixed
        ivec r_new = r.subvec(1, r.n_elem - 1);
        ivec r_app=r_new;

        // swap the elements i and j
        r_new[i]=r_app[j];
        r_new[j]=r_app[i];

        r_new.insert_rows(0,1);
        r_new(0)=r(0);
        _idxs=r;
    } 

    void block_shift(int k, int m, int n){

        // ------------------------- PARAMETERS ----------------------------
        // k start of the block, m length of the block, n length of the jump
        // -----------------------------------------------------------------
        ivec r=_idxs;
        
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
        
        _idxs=r;
}

    void block_permutation(int m, int k, int n){
        
        // ------------------------- PARAMETERS ----------------------------
        // m length of the block, k start of the block1, n start of the block2
        // -----------------------------------------------------------------
    
        ivec r=_idxs;

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
        _idxs=r;
    }

    void inversion(size_t startIdx, size_t endIdx) {
        ivec r=_idxs;
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
        _idxs=r;
    }

private:
    map<int,arma::vec> _indi;
    ivec _idxs;
    int _dim; 
    Random *_rnd;
};

#endif // INDIVIDUE_H

#ifndef POPULATION_H
#define POPULATION_H

class Population : Individue {
public:
    Population() : Individue(), _dim_indi(0),_dim_pop(0),_rnd(nullptr) {}

    // Initialization of the population on a circle
    void initialize_circle(int dim_indi, int dim_pop, double r,Random &rnd) {
        _dim_indi=dim_indi;
        _dim_pop=dim_pop;
        _rnd = &rnd;
        pop.set_size(_dim_pop);
        //first individue
        pop[0].initialize_circle(_dim_indi, r,*_rnd);
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
        pop[0].initialize_square(_dim_indi, l,*_rnd);
        for(size_t i = 0; i < _dim_pop; ++i) {
            pop[i] = pop[0];
            pop[i].shuffleidxs();
        }
    }
    
    // Operators
    Individue& operator[](size_t index) {
        if (index >= pop.n_elem) {
            throw std::out_of_range("Index out of range");
        }
        return pop[index];
    }
    
    // Get methods
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
    
    int get_pop_size() const { return _dim_pop;}
    
    int get_indi_size() const { return _dim_indi;}       

    // Print methods
    void print() const {
        for (size_t i = 0; i < pop.n_rows; ++i) {
            cout << "Individue " << i << ":" << endl;
            pop[i].print_positions();
            cout << endl;
        }
    }
    
    void printIdxs() const {
        for (size_t i = 0; i < pop.n_rows; ++i) {
            cout << "Individue " << i << ": ";
            pop(i, 0).print_idxs();
            cout << endl;
        }
    }
     
    void print_population() {
        for (size_t i = 0; i < _dim_pop; ++i) {
            cout << "Individue " << i << ": ";
            pop[i].print_idxs();
            cout << endl;
        }
    }

    // Print the losses of the population
    void print_pop_loss(){
        for (size_t i = 0; i < pop.n_elem; ++i) {
            cout << "Individue " << i << " loss: " << pop(i, 0).get_indi_Loss() << endl;
        }
    }

    /*void print_full_population(ofstream &file) {
        for (size_t i = 0; i < _dim_pop; ++i) {
            file << " Individue " << i << ": ";
            pop(i, 0).print_idxs(file);
            file << "Loss: " << pop(i, 0).get_indi_Loss() << " ";
            pop(i,0).print_positions(file);
            file << endl;
        }
    } */
    // Sort the population based on the loss
    void sort(){ 
        arma::field<Individue> new_pop;
        new_pop.set_size(_dim_pop);
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
            if (!pop[i].checkIndividual()) {
                pop[i].print_idxs();
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
            _sel=static_cast<int>((double)_dim_pop * (pow(_rnd->Rannyu(), 10)));
            new_pop[i] = pop[_sel];
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
        std::pair<ivec,ivec> new_indi_idxs;
        // Perform crossover with a probability
        for(int index1=0;index1<_dim_pop;index1++){
            rand_cross=_rnd->Rannyu();
            index2=index1;
            sindi=pop[index1];
            if(rand_cross<rnd_cross){
                //cerr << "crossover" << endl;
                //non far fare crossover tra due individui uguali
                while(arma::approx_equal(pop[index1].get_idxs(), pop[index2].get_idxs(), "absdiff", 1e-10)){index2=(int)_rnd->Rannyu(0,_dim_pop);}
                double rand=_rnd->Rannyu();
                //cerr << "rand: " << rand << endl;
                //cerr << "old indi:"<< endl;
                //pop(index1,0).printIdxs();
                //pop(index2,0).printIdxs();

                new_indi_idxs=crossover(pop[index1].get_idxs(),pop[index2].get_idxs(),rand);
                //cerr << "new indi:"<< endl;
                //new_indi.first.printIdxs();
                //new_indi.second.printIdxs();
                if(_rnd->Rannyu()<0.5){sindi.set_idxs(new_indi_idxs.first);}else{sindi.set_idxs(new_indi_idxs.second);}
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


private:
    arma::field<Individue>  pop;
    int _dim_indi, _dim_pop, _sel; Random *_rnd;
};

#endif // POPULATION_H
