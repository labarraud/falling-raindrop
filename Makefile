# Compilateur utilise
CC = g++
 
# Les différents répertoires contenant respectivement les fichiers : Sources *.cxx, Headers *.hxx, Objets *.o, l'exécutable
SRCDIR=Src
HEADDIR=Include
LIBDIR=Object
 
# options en mode optimise
OPTIM_FLAG = -O3 -std=c++11 -Wall -Woverloaded-virtual

# options en mode debug
DEBUG_FLAG = -g -std=c++11 -Wall -Woverloaded-virtual

LINKER   = g++ -o

CFLAGS   = -std=c++11 -Wall -I.

LFLAGS   = -Wall -I. -lm
# L'exécutable
BIN=run
SOURCES  := $(wildcard $(SRCDIR)/*.cxx)
INCLUDES := $(wildcard $(HEADDIR)/*.hxx)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.cxx=$(LIBDIR)/%.o)

rm       = rm -f

optim :	$(OBJECTS)
	$(LINKER) $(BIN) $(OPTIM_FLAG) $(OBJECTS)
	@echo "Linking optimum complete!"

debug : $(OBJECTS)
	$(LINKER) $(BIN) $(DEBUG_FLAG) $(OBJECTS)
	@echo "Linking debug complete!"

$(BIN): $(OBJECTS)
	$(LINKER) $(BIN) $(LFLAGS) $(OBJECTS)

$(OBJECTS): $(LIBDIR)/%.o : $(SRCDIR)/%.cxx
	$(CC) $(CFLAGS) -c $< -o $@

# Nettoyage des objets => Tout sera recompiler !
clean:
	rm $(LIBDIR)/*.o

# Nettoyage complet => clean + effacement du l'exécutable
Clean: clean
	rm $(BIN)



