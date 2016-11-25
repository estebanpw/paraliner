#include <inttypes.h>

#define READBUF 50000000 // 50 MB
#define MAXLID 50
#define STARTSIZE 50000000 
#define MATCH 4


struct Sequence{ 
  char id[MAXLID];
  char * data;
};

struct cell{
	int64_t score;
	uint64_t xfrom;
	uint64_t yfrom;
};
struct positioned_cell{
	int64_t score;
	uint64_t xpos;
	uint64_t ypos;
};