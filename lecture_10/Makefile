# Define the target executable name
TARGET = Exercise_10 TSP_Italy
# Define the number of processes to run
NP = 4
# Define the compiler
CC = mpic++
# Define the compiler flags
CFLAGS = -O3 --std=c++17
AFLAGS = -larmadillo
# Define the source, include, object, and function directories
SRC_DIR := ../classes/src
RND_DIR := ../random_number
INC_DIR := ../classes/inc
OBJ_DIR := ../classes/obj
FUNC := ../functions
# list of source files
SRCS := $(wildcard $(SRC_DIR)/*.cpp)
OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRCS))
# Expression to get the list of classes
CLASSES := $(patsubst $(SRC_DIR)/%,$(INC_DIR)/%,$(basename $(SRCS)))

all: $(TARGET)

Exercise_10 : $(OBJS) $(OBJ_DIR)/Exercise_10.o $(RND_DIR)/random.o $(FUNC)/functions.o
	$(CC)  -o $@ $^ $(CFLAGS) $(AFLAGS)
# rule to link the random number generator
$(RND_DIR)/random.o: $(RND_DIR)/random.cpp $(RND_DIR)/random.h
	$(CC) $(CFLAGS) -c -o $@ $<
# rule to link the functions
$(FUNC)/functions.o: $(FUNC)/functions.cpp $(FUNC)/functions.hpp
	$(CC) $(CFLAGS) -c -o $@ $<
# rule to link the file Exercise_10.o
$(OBJ_DIR)/Exercise_10.o : Exercise_10.cpp 
	$(CC) -c Exercise_10.cpp -o $(OBJ_DIR)/Exercise_10.o $(CFLAGS)
# rule to link the object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(INC_DIR)/%.h
	$(CC) $(CFLAGS) -c -o $@ $<

TSP_Italy : $(OBJS) $(OBJ_DIR)/TSP_Italy.o $(RND_DIR)/random.o $(FUNC)/functions.o
	$(CC)  -o $@ $^ $(CFLAGS) $(AFLAGS)
# rule to link the file TSP_Italy.o
$(OBJ_DIR)/TSP_Italy.o : TSP_Italy.cpp 
	$(CC) -c TSP_Italy.cpp -o $(OBJ_DIR)/TSP_Italy.o $(CFLAGS)

run_square: Exercise_10
	mpirun -np $(NP) ./Exercise_10

run_square_migration: Exercise_10
	mpirun -np $(NP) ./Exercise_10 migration

run_circle: Exercise_10
	mpirun -np $(NP) ./Exercise_10 circle

run_circle_migration: Exercise_10
	mpirun -np $(NP) ./Exercise_10 circle migration

run_city: TSP_Italy
	mpirun -np $(NP) ./TSP_Italy

run_city_migration: TSP_Italy
	mpirun -np $(NP) ./TSP_Italy migration
