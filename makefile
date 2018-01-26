# *****************************************************
# Variables to control Makefile operation

CC = g++
CFLAGS  = -g -Wall

# ****************************************************
# Targets needed to bring the executable up to date

# Nessie 
nessie: main.o FastaUtilities.o Nessie.o LinkedlistKmer.o HashTable.o BitArray/bit_array.o bitscan/tables.o bitscan/bitboards.o bitscan/bitboardn.o bitscan/bitboard.o bitscan/bbsentinel.o
	$(CC) $(CFLAGS) -o nessie main.o Nessie.o FastaUtilities.o LinkedlistKmer.o HashTable.o bit_array.o tables.o bitboards.o bitboardn.o bitboard.o bbsentinel.o
	@echo ' '
	@echo 'Successfully built nessie!'
	@echo ' '
	
main.o: src/main.cpp src/Nessie.h src/FastaUtilities.h
	$(CC) $(CFLAGS) -c src/main.cpp

Nessie.o: src/Nessie.cpp src/Nessie.h src/Functions.h src/LinkedlistKmer.h src/HashTable.h src/BitArray/bit_array.h src/bitscan/tables.h src/bitscan/bitboards.h src/bitscan/bitboardn.h src/bitscan/bitboard.h src/bitscan/bbsentinel.h
	$(CC) $(CFLAGS) -c src/Nessie.cpp
	
FastaUtilities.o: src/FastaUtilities.cpp src/FastaUtilities.h
	$(CC) $(CFLAGS) -c src/FastaUtilities.cpp

LinkedlistKmer.o: src/LinkedlistKmer.cpp src/LinkedlistKmer.h src/Functions.h
	$(CC) $(CFLAGS) -c src/LinkedlistKmer.cpp
	
HashTable.o: src/HashTable.cpp src/HashTable.h src/LinkedlistKmer.h src/Functions.h
	$(CC) $(CFLAGS) -c src/HashTable.cpp
	
#bitscan/bitscan.o
bitscan/tables.o: src/bitscan/tables.cpp src/bitscan/tables.h src/bitscan/bbtypes.h src/bitscan/config.h
	$(CC) $(CFLAGS) -c src/bitscan/tables.cpp
	
bitscan/bitboards.o: src/bitscan/bitboards.cpp src/bitscan/bitboards.h src/bitscan/bbobject.h src/bitscan/bitboard.h
	$(CC) $(CFLAGS) -c src/bitscan/bitboards.cpp

bitscan/bitboardn.o: src/bitscan/bitboardn.cpp src/bitscan/bitboardn.h src/bitscan/bbobject.h src/bitscan/bitboard.h
	$(CC) $(CFLAGS) -c src/bitscan/bitboardn.cpp
	
bitscan/bitboard.o: src/bitscan/bitboard.cpp src/bitscan/bitboard.h src/bitscan/tables.h
	$(CC) $(CFLAGS) -c src/bitscan/bitboard.cpp
	
bitscan/bbsentinel.o: src/bitscan/bbsentinel.cpp src/bitscan/bbsentinel.h src/bitscan/bbintrinsic.h src/bitscan/bbalg.h
	$(CC) $(CFLAGS) -c src/bitscan/bbsentinel.cpp

#BitArray/bit_array.o
BitArray/bit_array.o: src/BitArray/bit_array.c src/BitArray/bit_array.h src/BitArray/bit_macros.h
	$(CC) $(CFLAGS) -c src/BitArray/bit_array.c

# Comments


# Clean
clean: 
	$(RM) *.o 



