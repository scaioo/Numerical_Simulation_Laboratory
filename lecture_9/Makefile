CC = g++
CFLAGS = -O3 --std=c++17
AFLAGS = -larmadillo

SRC_DIR := ../classes/src
RND_DIR := ../random_number
INC_DIR := ../classes/inc
OBJ_DIR := ../classes/obj
FUNC := ../functions

# List of source files
SRCS := $(wildcard $(SRC_DIR)/*.cpp)
OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRCS))
# regular expression to get the name of the classes
CLASSES := $(patsubst $(SRC_DIR)/%,$(INC_DIR)/%,$(basename $(SRCS)))

all: Exercise_09 

Exercise_09 : $(OBJS) $(OBJ_DIR)/Exercise_09.o $(RND_DIR)/random.o $(FUNC)/functions.o
	$(CC)  -o $@ $^ $(CFLAGS) $(AFLAGS)
# Rule to compile the random number generator
$(RND_DIR)/random.o: $(RND_DIR)/random.cpp $(RND_DIR)/random.h
	$(CC) $(CFLAGS) -c -o $@ $<
# Rule to compile the functions
$(FUNC)/functions.o: $(FUNC)/functions.cpp $(FUNC)/functions.hpp
	$(CC) $(CFLAGS) -c -o $@ $<
# Rule to compile the Exercise_09.o
$(OBJ_DIR)/Exercise_09.o : Exercise_09.cpp 
	$(CC) -c Exercise_09.cpp -o $(OBJ_DIR)/Exercise_09.o $(CFLAGS)
# General rule to compile the classes
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(INC_DIR)/%.h
	$(CC) $(CFLAGS) -c -o $@ $<




clean: 
	rm -rf $(OBJ_DIR)/*.o Exercise_09 seed.out
