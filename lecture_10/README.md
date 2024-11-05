# Lecture 10

In this lecture, there are two main `.cpp` files:

- **Exercise_10.cpp**: Simulates city configurations either within a square or along a circumference, with options for running the simulation with or without migrations.
- **TSP_Italy.cpp**: Simulates a path optimization across Italian cities, with options for migration.

## Directory Structure

- **Results**: Contains all output data and visualizations:
  - `.csv` files with path coordinates
  - `.png` files of generated visualizations

## Commands

- `make`: Compiles all exercises
- `make Exercise_10`: Compiles `Exercise_10.cpp`
- `make TSP_Italy`: Compiles `TSP_Italy.cpp`
  
### Running Simulations

- `make run_circle`: Executes the circumference simulation without migrations
- `make run_circle_migration`: Executes the circumference simulation with migrations
- `make run_square`: Executes the square simulation without migrations
- `make run_square_migrations`: Executes the square simulation with migrations
- `make run_city`: Executes the Italian cities simulation without migrations
- `make run_city_migration`: Executes the Italian cities simulation with migrations

