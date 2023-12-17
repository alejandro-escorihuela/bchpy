# 17-12-2023
# alex

OBJ = recurABc.so # Objectes
COM = gcc # Compilador
LIB = -lm # Llibreries (-l*, -L*, -I*)
AVS = -W -Wall -Wextra -ansi -pedantic # Avisos
OPT = -march=native -O2 -pipe -std=c17 -shared -fPIC # Optimitzacio
OPC = $(OPT) $(AVS) $(LIB) # Opcions

all: $(OBJ)

recurABc.so: recurABc.c recurABc.h
	@echo Compilant $<
	$(COM) $(OPC) $< -o $@

clean:
	@echo Netejant...
	rm -rf $(OBJ) *~
