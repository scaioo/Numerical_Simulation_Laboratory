CC = g++
CFLAGS = -Wall -O3 --std=c++17
AFLAGS = -larmadillo
RND_DIR := ../random_number
FUNC := ../functions
SRC_DIR := ../classes/src
INC_DIR := ../classes/inc
OBJ_DIR := ../classes/obj

# Lista dei file sorgente per ciascuna classe
SRCS := $(wildcard $(SRC_DIR)/*.cpp)
OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRCS))
# Espressione regolare per estrarre i nomi delle classi
CLASSES := $(patsubst $(SRC_DIR)/%,$(INC_DIR)/%,$(basename $(SRCS)))


All: Exercise_02.1 Exercise_02.2

Exercise_02.1: $(OBJS) $(OBJ_DIR)/Exercise_02.1.o $(RND_DIR)/random.o $(FUNC)/functions.o
	$(CC) $(CFLAGS) -o $@ $^ $(AFLAGS)

# Regola per la compilazione delle functions
$(FUNC)/functions.o: $(FUNC)/functions.cpp $(FUNC)/functions.hpp
	$(CC) $(CFLAGS) -c -o $@ $<
# Regola per la compilazione della classe random
$(RND_DIR)/random.o: $(RND_DIR)/random.cpp $(RND_DIR)/random.h
	$(CC) $(CFLAGS) -c -o $@ $<
# Regola per la compilazione del file oggetto
$(OBJ_DIR)/Exercise_02.1.o : Exercise_02.1.cpp
	$(CC) -c Exercise_02.1.cpp -o $(OBJ_DIR)/Exercise_02.1.o $(CFLAGS)
# Regola generica per compilare i file oggetto
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(INC_DIR)/%.h
	$(CC) $(CFLAGS) -c -o $@ $<

Exercise1: Exercise_02.1
Exercise1_exe: Exercise_02.1
	./Exercise_02.1

Exercise_02.2: $(OBJS) $(OBJ_DIR)/Exercise_02.2.o $(RND_DIR)/random.o $(FUNC)/functions.o -larmadillo
	$(CC) $(CFLAGS) -o $@ $^ 

# Regola per la compilazione del file oggetto
$(OBJ_DIR)/Exercise_02.2.o : Exercise_02.2.cpp
	$(CC) -c Exercise_02.2.cpp -o $(OBJ_DIR)/Exercise_02.2.o $(CFLAGS)
# Regola generica per compilare i file oggetto
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(INC_DIR)/%.h
	$(CC) $(CFLAGS) -c -o $@ $<

Exercise2: Exercise_02.2
Exercise2_exe: Exercise_02.2
	./Exercise_02.2
clean:
	rm -rf $(OBJ_DIR)/*.o Exercise_02.1 seed.out
