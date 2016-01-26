# Compilateur utilise
CC = g++

# options en mode optimise
OPTIM_FLAG = -O3 -std=c++11 -Wall -Woverloaded-virtual

# options en mode debug
DEBUG_FLAG = -DLINALG_DEBUG -g -std=c++11 -Wall -Woverloaded-virtual

# executable produit
PROG = run

# fichier source a compiler
SRC = main.cxx

# par defaut on compile en optimise
optim : $(SRC)
	$(CC) $(SRC) $(OPTIM_FLAG) -o test.x
	mv test.x $(PROG)

debug : $(SRC)
	$(CC) $(SRC) $(DEBUG_FLAG) -o test.x
	mv test.x $(PROG)
