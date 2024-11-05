CC = g++
CFLAGS = -Wall -O3 --std=c++17
AFLAGS = -larmadillo
RND_DIR := ../random_number
FUNC := ../functions
SRC_DIR := ../classes/src
INC_DIR := ../classes/inc
OBJ_DIR := ../classes/obj

# List of source files
SRCS := $(wildcard $(SRC_DIR)/*.cpp)
OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRCS))
# Rgular expression to get the name of the classes
CLASSES := $(patsubst $(SRC_DIR)/%,$(INC_DIR)/%,$(basename $(SRCS)))

All: Exercise_03.1

Exercise_03.1: $(OBJS) $(OBJ_DIR)/Exercise_03.1.o $(RND_DIR)/random.o $(FUNC)/functions.o
	$(CC) $(CFLAGS) -o $@ $^ $(AFLAGS)

# Rule to compile the functions
$(FUNC)/functions.o: $(FUNC)/functions.cpp $(FUNC)/functions.hpp
	$(CC) $(CFLAGS) -c -o $@ $<
# Rule to compile the random number generator
$(RND_DIR)/random.o: $(RND_DIR)/random.cpp $(RND_DIR)/random.h
	$(CC) $(CFLAGS) -c -o $@ $<
# Rule to compile the Ecercise_03.1
$(OBJ_DIR)/Exercise_03.1.o : Exercise_03.1.cpp
	$(CC) -c Exercise_03.1.cpp -o $(OBJ_DIR)/Exercise_03.1.o $(CFLAGS)
# General rule to compile the objects
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(INC_DIR)/%.h
	$(CC) $(CFLAGS) -c -o $@ $<