CC=gcc
CXX=g++
CFLAGS=-O3 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -Wall -DVERBOSE
BIN=../bin

all: nwalign

nwalign: NW.c
	$(CC) $(CFLAGS) common.c NW.c -lpthread -o $(BIN)/nwalign
clean:
	rm -rf $(BIN)/nwalign
