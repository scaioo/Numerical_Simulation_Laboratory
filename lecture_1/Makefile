CC = g++
CFLAGS = -Wall -O3 --std=c++17
AFLAGS = -larmadillo
RND_DIR := ../random_number
FUNC := ../functions
SRC_DIR := ../classes/src
INC_DIR := ../classes/inc
OBJ_DIR := ../classes/obj

# List for the source files
SRCS := $(wildcard $(SRC_DIR)/*.cpp)
OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRCS))
# Regular exepression to extract classes names
CLASSES := $(patsubst $(SRC_DIR)/%,$(INC_DIR)/%,$(basename $(SRCS)))


All: Exercise_01.1 Exercise_01.2 Exercise_01.3

Exercise_01.1: $(OBJS) $(OBJ_DIR)/Exercise_01.1.o $(RND_DIR)/random.o $(FUNC)/functions.o
	$(CC) $(CFLAGS) -o $@ $^ $(AFLAGS)
# Rule for the compilation of the functions class
$(FUNC)/functions.o: $(FUNC)/functions.cpp $(FUNC)/functions.hpp
	$(CC) $(CFLAGS) -c -o $@ $<
# Rule for the compilation of the random class
$(RND_DIR)/random.o: $(RND_DIR)/random.cpp $(RND_DIR)/random.h
	$(CC) $(CFLAGS) -c -o $@ $<
# Rule for the compilation of the exercise file
$(OBJ_DIR)/Exercise_01.1.o : Exercise_01.1.cpp
	$(CC) -c Exercise_01.1.cpp -o $(OBJ_DIR)/Exercise_01.1.o $(CFLAGS)
# Rule for the compilation of the object file
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(INC_DIR)/%.h
	$(CC) $(CFLAGS) -c -o $@ $<

Exercise_01.2: $(OBJS) $(OBJ_DIR)/Exercise_01.2.o $(RND_DIR)/random.o $(FUNC)/functions.o -larmadillo
	$(CC) $(CFLAGS) -o $@ $^ 
# Rule for the compilation of the exercise file
$(OBJ_DIR)/Exercise_01.2.o : Exercise_01.2.cpp
	$(CC) -c Exercise_01.2.cpp -o $(OBJ_DIR)/Exercise_01.2.o $(CFLAGS)
# Rule for the compilation of the object file
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(INC_DIR)/%.h
	$(CC) $(CFLAGS) -c -o $@ $<

Exercise_01.3: $(OBJS) $(OBJ_DIR)/Exercise_01.3.o $(RND_DIR)/random.o $(FUNC)/functions.o -larmadillo
	$(CC) $(CFLAGS) -o $@ $^ 
# Rule for the compilation of the exercise file
$(OBJ_DIR)/Exercise_01.3.o : Exercise_01.3.cpp
	$(CC) -c Exercise_01.3.cpp -o $(OBJ_DIR)/Exercise_01.3.o $(CFLAGS)
# Rule for the compilation of the object file
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(INC_DIR)/%.h
	$(CC) $(CFLAGS) -c -o $@ $<



clean:
	rm -rf $(OBJ_DIR)/*.o Exercise_01.* seed.out
