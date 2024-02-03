
# 17-12-2023
# alex

OBJ = metodeABA.so metodeXX.so recurAB.so recurABc.so recurXX.so recurXXc.so esq.so # Objectes
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

metodeABA.so: metodeABA.c metodeABA.h recurAB.h recurABc.h esq.h
	@echo Compilant $<
	$(COM) $(OPC) $< -o $@

recurAB.so: recurAB.c recurAB.h esq.h
	@echo Compilant $<
	$(COM) $(OPC) $< -o $@

recurABc.so: recurABc.c recurABc.h esq.h
	@echo Compilant $<
	$(COM) $(OPC) $< -o $@

metodeXX.so: metodeXX.c metodeXX.h recurXX.h recurXXc.h esq.h
	@echo Compilant $<
	$(COM) $(OPC) $< -o $@

recurXX.so: recurXX.c recurXX.h esq.h
	@echo Compilant $<
	$(COM) $(OPC) $< -o $@

recurXXc.so: recurXXc.c recurXXc.h
	@echo Compilant $<
	$(COM) $(OPC) $< -o $@

esq.so: esq.c esq.h
	@echo Compilant $<
	$(COM) $(OPC) $< -o $@

clean:
	@echo Netejant...
	rm -rf $(OBJ) *~
