# Numerical Simulation Laboratory

This project contains **12 lectures** where are perfomed different types of numerical simulations.
The code is divided in 8 directories named `lecture_<number_lecture>`. In each lecture there are: the `.cpp` files,the `.exe` files the `Makefile`. All the outputs are contained in the relatice `Results` subdirectory.
Lectures 4 6 and 7 are merged into a unique directory because make use of the *NSL_Simulator*, a class to perform MD simulations. 
In the `Jupyter Notebook` direcory are stored all the notebook where are presented the results obtained.
In the `functions` and `classes` directories there are all the classes and functions used to solve the exercises directories there are all the classes and functions used to solve the exercises.

-  `lecture_1` topics:
	- Integral calculation with MonteCarlo sampling
	- Block Average method
	- Sampling different distribution with Random number generator
	- Central limit theorem

- `lecture_2` topics:
	- Importance sampling
	- 3D discrete and continouos Random Walks
	- Diffusive motion

- `lecture_3` topics:
	- spot price modeling
	- geometric brownian motion
	- Black Sholes Theory

- `lecture_4` topics:
	- Virial theorem in MD simulations
	- NVE and NVT MD simulations
	- 1D Ising Model
	- MD equilibration

- `lecture_5` topics:
	- Metropolis algorithm
	- quantum wavefunctions sampling
	- Equilibration in Metropolis

- `lecture_6` topics:
	- 1D Ising model simulation
	- Metropolis and Gibbs sampling for MD simulations
	- magnetization-susceptibility-heat capacity analysis

- `lecture_7` topics:
	- tail corrections
	- Metropolis acceptance ratio in equilibrated Monte Carlo simulation
	- autocorrelation function
	- radial distribution function

- `lecture_8` topics:
	- Hydroge atom with Variational Monte Carlo
	- Ground and Excited state energy 
	- Simulated Annealing

- `lecture_9` topics:
	- Genetic Algorithm
	- Travelling Salesman Problem

- `lecture_10` topics:
	- Parallel computing
	- MPI library

- `lecture_11` topics:
	- FFNN in keras
	- Hyperparameters optimization

- `lecture_12` topics:
	- DNN and CNN in Keras
	- Optimizers Analysis

## Initialization

In each folder there are *.cpp  and *.exe files. In each lecture directory there is a `Results` directory where all the outputs are saved.
to compile a lecture you should follows the following steps:
 ```bash
cd lecture_<lecture_number>/
`make`
./Exercise_<lecture_number>.<exercise_number>
```

### Prerequisites

- `g++` compiler supporting `c++17`
- `Armadillo` library
- `MPI` for parallel computing in `lecture_10`




