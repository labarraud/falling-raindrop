# Compilateur utilise
CC = g++

# options en mode optimise
OPTIM_FLAG = -O3 -std=c++11 -Wall -Woverloaded-virtual

# options en mode debug
DEBUG_FLAG = -DLINALG_DEBUG -g -std=c++11 -Wall -Woverloaded-virtual

# executable produit
PROG = run

# fichier source a compiler
SRC = TimeScheme.o spacescheme.o DiffusionConvectionProblem.o particle.o velocity.o Matrix.o main.o


all:$(PROG)

# par defaut on compile en optimise
$(PROG):$(SRC)
	$(CC) $(SRC) $(OPTIM_FLAG) -o $(PROG)

%.o:%.cxx
	$(CC) $(DEBUG_FLAG) -c $<

main.o:main.cc
	$(CC) $(DEBUG_FLAG) -c main.cc


optim : $(SRC)
	$(CC) $(SRC) $(OPTIM_FLAG) -o test.x
	mv test.x $(PROG)
	
debug : $(SRC)
	$(CC) $(SRC) $(DEBUG_FLAG) -o test.x
	mv test.x $(PROG)
	
clean:
	rm -f *.o
	rm run
