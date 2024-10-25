CC = g++
CFLAGS = -O3 --std=c++17
AFLAGS = -larmadillo

RND_DIR := ../random_number
SRC_DIR := simulator_classes/src
INC_DIR := simulator_classes/inc
OBJ_DIR := simulator_classes/obj
FUNC := ../functions

# Lista dei file sorgente per ciascuna classe
SRCS := $(wildcard $(SRC_DIR)/*.cpp)
OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRCS))
# Espressione regolare per estrarre i nomi delle classi
CLASSES := $(patsubst $(SRC_DIR)/%,$(INC_DIR)/%,$(basename $(SRCS)))

all: equilibration Exercise_04.2 Exercise_06.1 Exercise_07.2 Exercise_07.4

equilibration : $(OBJS) $(OBJ_DIR)/equilibration.o $(FUNC)/functions.o
	$(CC)  -o $@ $^ $(CFLAGS) $(AFLAGS) 
# Regola per la compilazione delle functions
$(FUNC)/functions.o: $(FUNC)/functions.cpp $(FUNC)/functions.hpp
	$(CC) $(CFLAGS) -c -o $@ $<
# Regola per la compilazione del file oggetto Exercise_04.2.o
$(OBJ_DIR)/equilibration.o : equilibration.cpp 
	$(CC) -c equilibration.cpp -o $(OBJ_DIR)/equilibration.o $(CFLAGS)
# Regola generica per compilare i file oggetto
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(INC_DIR)/%.h
	$(CC) $(CFLAGS) -c -o $@ $<

Exercise_04.2: $(OBJS) $(OBJ_DIR)/Exercise_04.2.o $(FUNC)/functions.o
	$(CC)  -o $@ $^ $(CFLAGS) $(AFLAGS)
# Regola per la compilazione del file oggetto Exercise_04.2.eq.o
$(OBJ_DIR)/Exercise_04.2.o : Exercise_04.2.cpp
	$(CC) -c Exercise_04.2.cpp -o $(OBJ_DIR)/Exercise_04.2.o $(CFLAGS)
# Regola generica per compilare i file oggetto
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(INC_DIR)/%.h
	$(CC) $(CFLAGS) -c -o $@ $<

Exercise_06.1: $(OBJS) $(OBJ_DIR)/Exercise_06.1.o $(FUNC)/functions.o
	$(CC)  -o $@ $^ $(CFLAGS) $(AFLAGS)
# Regola per la compilazione del file oggetto Exercise_04.2.eq.o
$(OBJ_DIR)/Exercise_06.1.o : Exercise_06.1.cpp
	$(CC) -c Exercise_06.1.cpp -o $(OBJ_DIR)/Exercise_06.1.o $(CFLAGS)
# Regola generica per compilare i file oggetto
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(INC_DIR)/%.h
	$(CC) $(CFLAGS) -c -o $@ $<

Exercise_07.2: $(OBJS) $(OBJ_DIR)/Exercise_07.2.o $(FUNC)/functions.o
	$(CC)  -o $@ $^ $(CFLAGS) $(AFLAGS)
# Regola per la compilazione del file oggetto Exercise_04.2.eq.o
$(OBJ_DIR)/Exercise_07.2.o : Exercise_07.2.cpp
	$(CC) -c Exercise_07.2.cpp -o $(OBJ_DIR)/Exercise_07.2.o $(CFLAGS)
# Regola generica per compilare i file oggetto
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(INC_DIR)/%.h
	$(CC) $(CFLAGS) -c -o $@ $<

Exercise_07.4: $(OBJS) $(OBJ_DIR)/Exercise_07.4.o $(FUNC)/functions.o
	$(CC)  -o $@ $^ $(CFLAGS) $(AFLAGS)
# Regola per la compilazione del file oggetto Exercise_04.2.eq.o
$(OBJ_DIR)/Exercise_07.4.o : Exercise_07.4.cpp
	$(CC) -c Exercise_07.4.cpp -o $(OBJ_DIR)/Exercise_07.4.o $(CFLAGS)
# Regola generica per compilare i file oggetto
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(INC_DIR)/%.h
	$(CC) $(CFLAGS) -c -o $@ $<


clean :
	rm -f $(OBJ_DIR)/*.o Exercise_04.2 Exercise_04.2.eq Exercise_06.1 Exercise_07.2
remove:
	rm ../OUTPUT/*.*
	rm ../OUTPUT/CONFIG/*.*
