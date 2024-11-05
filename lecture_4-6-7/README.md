# lecture_4-6-7

## Structure
This lecture are merged into a unique directory because make use of the `NSL_simulator` to perform MD simulations.
- There are three directory named `input_<number_exercise>`. In this directories there are the input files to start the simulations for the specific exercise. 
- The `equilibration.cpp` file is a script to perform equilibrations. This code is used to equilibrate a sistem and to save the final outputs in a specific output. 
The input files for the equilibration are situated in the `input_equilibration` directory.
- The `Results directory contains 4 subdirectories: `output_04.2` `output_06.1` `output_07.2` and `output_07.4` where are stored the output for the equilibrations and for the simulations.
	- In `output_04.2` can be found two directories: `equilibration` where there are the equilibration outputs at different temperatures for the three phases. In the other directory `equilibrated_simulation` there are the simulation started at the right temperature and after the relaxation of the system.

	- In `output_06.1`can be found two directories: `equilibration` and `measure`. In each directory there will be other subdirectories that divide equilibrations and simulations performed with *Metropolis* or *Gibbs* and with or without Magnetic Field at different temperatures. In the end there are four `files.txt` where all the ouptuts obtained are merged together in a unique file.

	- In `output_07.2`can be found two directories: `equilibration` and `simualtion_equilibrated` where can be found respectively the outpus of the equilibrations and the outputs of simulations for the three phases solid, liquid and gas.

	- In `output_07.2`can be found three directories: `equilibration_NVE` `NVE` and `NVT`. In `equilibration_NVE` there are the equilibration outputs for the three phases solid, liquid and gas in the **NVE** system. In the other two directories there are the results obtained from the simulations perfomed.

### Usage
- To compile all together the comman is `make`. Then the programs can be executed by ./Exercise_<exercise_number>. 
- To compile a single file the comman is `make Exercise_<lecture_number>.<exercise_number>` or `make equilibation`.
- To remove all thee objects files created the command is `make clean` 
