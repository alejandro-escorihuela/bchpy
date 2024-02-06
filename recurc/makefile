# 17-12-2023
# alex

OBJ = metodes.so recurAB.so recurABc.so recurXX.so recurXXc.so recurS2.so recurS2c.so recurS4.so recurS4c.so esq.so # Objectes
ICC = /opt/intel/oneapi/compiler/2024.0/bin/icx  # Compilador C intel
# Per al compilador d'intel primer:
# $ export LD_LIBRARY_PATH="/opt/intel/oneapi/compiler/latest/lib:/usr/local/lib"
GCC = gcc # Compilador C gnu
COM = $(GCC)
LIB = -lm # Llibreries (-l*, -L*, -I*)
AVS = -W -Wall -Wextra -ansi -pedantic # Avisos
OPT = -march=native -O2 --fast-math -pipe -std=c17 -shared -fPIC # Optimitzacio
OPC = $(OPT) $(AVS) $(LIB) # Opcions

all: $(OBJ)

metodes.so: metodes.c metodes.h recurAB.h recurABc.h recurXX.h recurXXc.h recurS2.h recurS2c.h recurS4.h recurS4c.h esq.h
	@echo Compilant $<
	$(COM) $(OPC) $< -o $@

recurAB.so: recurAB.c recurAB.h esq.h
	@echo Compilant $<
	$(COM) $(OPC) $< -o $@

recurABc.so: recurABc.c recurABc.h esq.h
	@echo Compilant $<
	$(COM) $(OPC) $< -o $@

recurXX.so: recurXX.c recurXX.h esq.h
	@echo Compilant $<
	$(COM) $(OPC) $< -o $@

recurXXc.so: recurXXc.c recurXXc.h
	@echo Compilant $<
	$(COM) $(OPC) $< -o $@

recurS2.so: recurS2.c recurS2.h esq.h
	@echo Compilant $<
	$(COM) $(OPC) $< -o $@

recurS2c.so: recurS2c.c recurS2c.h
	@echo Compilant $<
	$(COM) $(OPC) $< -o $@

recurS4.so: recurS4.c recurS4.h esq.h
	@echo Compilant $<
	$(COM) $(OPC) $< -o $@

recurS4c.so: recurS4c.c recurS4c.h
	@echo Compilant $<
	$(COM) $(OPC) $< -o $@

esq.so: esq.c esq.h
	@echo Compilant $<
	$(COM) $(OPC) $< -o $@

clean:
	@echo Netejant...
	rm -rf $(OBJ) *~
