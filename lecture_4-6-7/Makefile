CC = g++
CFLAGS = -O3 --std=c++17
AFLAGS = -larmadillo

RND_DIR := ../random_number
SRC_DIR := simulator_classes/src
INC_DIR := simulator_classes/inc
OBJ_DIR := simulator_classes/obj
FUNC := ../functions

# List of source and object files
SRCS := $(wildcard $(SRC_DIR)/*.cpp)
OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRCS))
# Regular expression to get the name of the classes
CLASSES := $(patsubst $(SRC_DIR)/%,$(INC_DIR)/%,$(basename $(SRCS)))

all: equilibration Exercise_04.2 Exercise_06.1 Exercise_07.2 Exercise_07.4

equilibration : $(OBJS) $(OBJ_DIR)/equilibration.o $(FUNC)/functions.o
	$(CC)  -o $@ $^ $(CFLAGS) $(AFLAGS) 
# Rule for the compilation of the object file functions.o
$(FUNC)/functions.o: $(FUNC)/functions.cpp $(FUNC)/functions.hpp
	$(CC) $(CFLAGS) -c -o $@ $<
# Rule for the compilation of the object file equilibration.o
$(OBJ_DIR)/equilibration.o : equilibration.cpp 
	$(CC) -c equilibration.cpp -o $(OBJ_DIR)/equilibration.o $(CFLAGS)
# Generic rule for the compilation of object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(INC_DIR)/%.h
	$(CC) $(CFLAGS) -c -o $@ $<

Exercise_04.2: $(OBJS) $(OBJ_DIR)/Exercise_04.2.o $(FUNC)/functions.o
	$(CC)  -o $@ $^ $(CFLAGS) $(AFLAGS)
# Rule for the compilation of the object file Exercise_04.2.eq.o
$(OBJ_DIR)/Exercise_04.2.o : Exercise_04.2.cpp
	$(CC) -c Exercise_04.2.cpp -o $(OBJ_DIR)/Exercise_04.2.o $(CFLAGS)
# Rule for the compilation of the object file functions.o
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(INC_DIR)/%.h
	$(CC) $(CFLAGS) -c -o $@ $<

Exercise_06.1: $(OBJS) $(OBJ_DIR)/Exercise_06.1.o $(FUNC)/functions.o
	$(CC)  -o $@ $^ $(CFLAGS) $(AFLAGS)
# Rule for the compilation of the object file Exercise_06.1.o
$(OBJ_DIR)/Exercise_06.1.o : Exercise_06.1.cpp
	$(CC) -c Exercise_06.1.cpp -o $(OBJ_DIR)/Exercise_06.1.o $(CFLAGS)
# Rule for the compilation of the object file functions.o
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(INC_DIR)/%.h
	$(CC) $(CFLAGS) -c -o $@ $<

Exercise_07.2: $(OBJS) $(OBJ_DIR)/Exercise_07.2.o $(FUNC)/functions.o
	$(CC)  -o $@ $^ $(CFLAGS) $(AFLAGS)
# Rule for the compilation of the object file Exercise_07.2.o
$(OBJ_DIR)/Exercise_07.2.o : Exercise_07.2.cpp
	$(CC) -c Exercise_07.2.cpp -o $(OBJ_DIR)/Exercise_07.2.o $(CFLAGS)
# Generic rule for the compilation of object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(INC_DIR)/%.h
	$(CC) $(CFLAGS) -c -o $@ $<

Exercise_07.4: $(OBJS) $(OBJ_DIR)/Exercise_07.4.o $(FUNC)/functions.o
	$(CC)  -o $@ $^ $(CFLAGS) $(AFLAGS)
# Rule for the compilation of the object file Exercise_07.4.o
$(OBJ_DIR)/Exercise_07.4.o : Exercise_07.4.cpp
	$(CC) -c Exercise_07.4.cpp -o $(OBJ_DIR)/Exercise_07.4.o $(CFLAGS)
# Generic rule for the compilation of object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(INC_DIR)/%.h
	$(CC) $(CFLAGS) -c -o $@ $<


clean :
	rm -f $(OBJ_DIR)/*.o Exercise_04.2 Exercise_06.1 Exercise_07.2
remove:
	rm ../OUTPUT/*.*
	rm ../OUTPUT/CONFIG/*.*
