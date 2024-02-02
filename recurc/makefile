# 17-12-2023
# alex

OBJ = metodeABA.so recurAB.so recurABc.so # Objectes
ICC = /opt/intel/oneapi/compiler/2024.0/bin/icx  # Compilador C intel
# Per al compilador d'intel primer:
# $ export LD_LIBRARY_PATH="/opt/intel/oneapi/compiler/latest/lib:/usr/local/lib"
GCC = gcc # Compilador C gnu
COM = $(GCC)
LIB = -lm # Llibreries (-l*, -L*, -I*)
AVS = -W -Wall -Wextra -ansi -pedantic # Avisos
OPT = -march=native -O2 -pipe -std=c17 -shared -fPIC # Optimitzacio
OPC = $(OPT) $(AVS) $(LIB) # Opcions

all: $(OBJ)

metodeABA.so: metodeABA.c metodeABA.h recurAB.h recurABc.h
	@echo Compilant $<
	$(COM) $(OPC) $< -o $@

recurAB.so: recurAB.c recurAB.h
	@echo Compilant $<
	$(COM) $(OPC) $< -o $@

recurABc.so: recurABc.c recurABc.h
	@echo Compilant $<
	$(COM) $(OPC) $< -o $@

clean:
	@echo Netejant...
	rm -rf $(OBJ) *~
