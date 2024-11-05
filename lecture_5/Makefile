CC = g++
CFLAGS = -O3 --std=c++17
AFLAGS = -larmadillo

SRC_DIR := ../classes/src
RND_DIR := ../random_number
INC_DIR := ../classes/inc
OBJ_DIR := ../classes/obj
FUNC := ../functions

# list of source and object files
SRCS := $(wildcard $(SRC_DIR)/*.cpp)
OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRCS))
# regular expression to get the name of the classes
CLASSES := $(patsubst $(SRC_DIR)/%,$(INC_DIR)/%,$(basename $(SRCS)))

all: Exercise_05.1 Exercise_05.1_far

Exercise_05.1 : $(OBJS) $(OBJ_DIR)/Exercise_05.1.o $(RND_DIR)/random.o $(FUNC)/functions.o
	$(CC)  -o $@ $^ $(CFLAGS) $(AFLAGS)

# Rule to compile the random.o file
$(RND_DIR)/random.o: $(RND_DIR)/random.cpp $(RND_DIR)/random.h
	$(CC) $(CFLAGS) -c -o $@ $<
# Rule to compile the functions.o file
$(FUNC)/functions.o: $(FUNC)/functions.cpp $(FUNC)/functions.hpp
	$(CC) $(CFLAGS) -c -o $@ $<
# Rule to compile the Exercise_05.1.o file
$(OBJ_DIR)/Exercise_05.1.o : Exercise_05.1.cpp 
	$(CC) -c Exercise_05.1.cpp -o $(OBJ_DIR)/Exercise_05.1.o $(CFLAGS)
# Rule to compile the classes
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(INC_DIR)/%.h
	$(CC) $(CFLAGS) -c -o $@ $<

Exercise_05.1_far : $(OBJS) $(OBJ_DIR)/Exercise_05.1_far.o $(RND_DIR)/random.o $(FUNC)/functions.o
	$(CC)  -o $@ $^ $(CFLAGS) $(AFLAGS)

$(OBJ_DIR)/Exercise_05.1_far.o : Exercise_05.1_far.cpp 
	$(CC) -c Exercise_05.1_far.cpp -o $(OBJ_DIR)/Exercise_05.1_far.o $(CFLAGS)

clean :
	rm $(OBJ_DIR)/*.o $(RND_DIR)/*.o $(FUNC)/*.o Exercise_05.1