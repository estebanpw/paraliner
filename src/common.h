
/*
	Returns +MATCH if symbols are equal, -MATH otherwise
*/
inline int64_t pam_dna_score(char c1, char c2);

/*
	Prints an error message and terminates the program
*/
void terror(const char * m);
/*
	Read a sequence from a fasta file
*/
uint64_t loadSeq(FILE * query, struct Sequence * s);

/*
	Buffered reading from file
*/
char buffered_fgetc(char *buffer, uint64_t *pos, uint64_t *read, FILE *f);

/*
	Calculates NW table with two rows and stores a cellpath of scores, identities, gaps and starting and ending positions
*/
void NW(char * X, uint64_t Xstart, uint64_t Xend, char * Y, uint64_t Ystart, uint64_t Yend, int iGap, int eGap, struct cell ** table, struct positioned_cell * mc, int show);
