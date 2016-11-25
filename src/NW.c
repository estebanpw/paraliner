/*********

File		NWscore2rows.c
Author		EPW <estebanpw@uma.es>
Description	Calculates a Needleman-Wunsch scores matrix using two rows and stores the number of identities and gaps, but does not retrieve the alignment.
		It is mainly inteded for use with reads vs genomes.

INPUT	<char * X>	Sequence X translated to numbers stored as one byte each letter
		<uint64_t Xstart>	Position in X sequence to start the alignment
		<uint64_t Xend>		Position in X sequence to end the alignment
		<char * Y>	Sequence Y translated to numbers stored as one byte each letter
		<uint64_t Xstart>	Position in Y sequence to start the alignment
		<uint64_t Xend>		Position in Y sequence to end the alignment
		<int iGap>		Penalty to open gap
		<int eGap>		Penalty to extend gap
		<int **PAM>		Matrix storing the match and miss rewards and penalties for each combination of letters
		<struct cell * ...>	Cell structs allocated from outside this function
		
RETURNS
		<struct cell bc>	Bottom cell with the best score

**********/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <unistd.h>
#include <inttypes.h>
#include <math.h>
#include <string.h>
#include <pthread.h>
#include "structs.h"
#include "common.h"

#define MAX_ROUTE 255
#define MAXALF 30
#define WINDOW 5


/*
gcc NWscore2rows.c  -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -O3 -o nwalign
./nwalign x y 8 0 pamDNAident.txt 0
*/

int main(int ac, char **av){
	
	FILE *f, *g, *out;
	
	if(ac != 7) terror("USE: ./nwalign seqx seqy alignment_file igap egap show");


	f = fopen64(av[1], "rt");
	g = fopen64(av[2], "rt");
	out = fopen64(av[3], "wt");
	if(f == NULL || g == NULL || out == NULL) terror("Could not open files");	

	int iGap = - atoi(av[4]);
	int eGap = - atoi(av[5]);
	int show = atoi(av[6]);
	

	struct Sequence sx, sy;
	sx.data = (char *) malloc(STARTSIZE*sizeof(char));
	sy.data = (char *) malloc(STARTSIZE*sizeof(char));

	if(sx.data == NULL || sy.data == NULL) terror("Could not allocate sequences");

	
	fprintf(stdout, "[INFO] Loading sequences\n");
	uint64_t xlen = loadSeq(f, &sx);
	fprintf(stdout, "[INFO] Loaded X sequence %s\n", sx.id);
	uint64_t ylen = loadSeq(g, &sy);
	fprintf(stdout, "[INFO] Loaded Y sequence %s\n", sy.id);

	fclose(f);
	fclose(g);

	fprintf(stdout, "[INFO] Length of sequence X: %"PRIu64"\n", xlen);
	fprintf(stdout, "[INFO] Length of sequence Y: %"PRIu64"\n", ylen);
	
	uint64_t maximum_len = max(xlen, ylen);

	char * reconstruct_X = (char *) malloc(maximum_len * sizeof(char));
	char * reconstruct_Y = (char *) malloc(maximum_len * sizeof(char));
	if(reconstruct_Y == NULL || reconstruct_X == NULL) terror("Could not allocate output alignment sequences");

	struct positioned_cell * mc = (struct positioned_cell *) malloc(ylen * sizeof(struct positioned_cell));
	struct cell ** table = (struct cell **) malloc(xlen * sizeof(struct cell *));
	uint64_t i,j;


	for(i=0;i<xlen;i++){
		table[i] = (struct cell *) malloc(ylen * sizeof(struct cell));
	}
    if(table == NULL) terror("Could not allocate NW table");



   	NW(sx.data, 0, xlen, sy.data, 0, ylen, iGap, eGap, table, mc, show);
   	backtrackingNW(sx.data, 0, xlen, sy.data, 0, ylen, table, reconstruct_X, reconstruct_Y);

   	i=0;
   	j=0;
   	while(i<maximum_len && j<maximum_len){
   		fprintf(out, "%"PRIu64"\t\t", i);
   		i += fwrite(&reconstruct_X[i], sizeof(char), ALIGN_LEN, out);
   		fprintf(out, " %"PRIu64" \n", i);
   		fprintf(out, "%"PRIu64"\t\t", j);
   		j += fwrite(&reconstruct_Y[j], sizeof(char), ALIGN_LEN, out);
   		fprintf(out, " %"PRIu64" \n", j);

   	}
	
   	free(mc);
	for(i=0;i<xlen;i++){
		free(table[i]);
	}
	free(table);
	free(sx.data);
	free(sy.data);
	free(reconstruct_X);
	free(reconstruct_Y);
	fclose(out);
	return 0;	
}


