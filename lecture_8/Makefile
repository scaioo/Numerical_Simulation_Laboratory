CC = g++
CFLAGS = -O3 --std=c++17
AFLAGS = -larmadillo

SRC_DIR := ../classes/src
RND_DIR := ../random_number
INC_DIR := ../classes/inc
OBJ_DIR := ../classes/obj
FUNC := ../functions

# Lista dei file sorgente per ciascuna classe
SRCS := $(wildcard $(SRC_DIR)/*.cpp)
OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRCS))
# Espressione regolare per estrarre i nomi delle classi
CLASSES := $(patsubst $(SRC_DIR)/%,$(INC_DIR)/%,$(basename $(SRCS)))

all: Exercise_08

Exercise_08 : $(OBJS) $(OBJ_DIR)/Exercise_08.o $(RND_DIR)/random.o $(FUNC)/functions.o
	$(CC)  -o $@ $^ $(CFLAGS) $(AFLAGS)
#regola per la compilazione del file oggetto random.o
$(RND_DIR)/random.o: $(RND_DIR)/random.cpp $(RND_DIR)/random.h
	$(CC) $(CFLAGS) -c -o $@ $<
# Regola per la compilazione delle functions
$(FUNC)/functions.o: $(FUNC)/functions.cpp $(FUNC)/functions.hpp
	$(CC) $(CFLAGS) -c -o $@ $<
# Regola per la compilazione del file oggetto Exercise_04.2.o
$(OBJ_DIR)/Exercise_08.o : Exercise_08.cpp 
	$(CC) -c Exercise_08.cpp -o $(OBJ_DIR)/Exercise_08.o $(CFLAGS)
# Regola generica per compilare i file oggetto
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(INC_DIR)/%.h
	$(CC) $(CFLAGS) -c -o $@ $<

clean: 
	rm -rf $(OBJ_DIR)/*.o Exercise_08 seed.out