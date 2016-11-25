#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <inttypes.h>
#include <limits.h>
#include <ctype.h>
#include "structs.h"
#include "common.h"
#define max(a,b)    (((a)>=(b)) ? (a):(b))
#define min(x,y)    (((x) < (y)) ? (x) : (y))



inline int64_t pam_dna_score(char c1, char c2){
    if(c1==c2) return MATCH;
    return (-MATCH);
}

void terror(const char * m) {
    printf("\nERR**%s***\n",m); 
    exit(-1);
}

uint64_t loadSeq(FILE * query, struct Sequence * s){
    char c;
    uint64_t length = 0, k = 0;

    uint64_t idx = 0, r = 0, reallocs = 1;
    char * temp_seq_buffer = NULL;
    if ((temp_seq_buffer = calloc(READBUF, sizeof(char))) == NULL) {
        terror("Could not allocate memory for read buffer");
    }
    //To force reading from the buffer
    idx = READBUF + 1;


    while((c = buffered_fgetc(temp_seq_buffer, &idx, &r, query))!= '>' && (!feof(query) || (feof(query) && idx < r))); //Get to ">"
    while(k < MAXLID && c !='\n' && c != ' '){
        s->id[k++] = c;
        c = buffered_fgetc(temp_seq_buffer, &idx, &r, query);
    }
    s->id[k] = '\0';
    while(c!='\n') c = buffered_fgetc(temp_seq_buffer, &idx, &r, query);
    c = buffered_fgetc(temp_seq_buffer, &idx, &r, query);
    while(c!='>' && (!feof(query) || (feof(query) && idx < r))) {
        c = toupper(c);
        if(c == 'A' || c == 'C' || c == 'T' || c == 'G')  s->data[length++]=c;
        if(length > reallocs*STARTSIZE){
            reallocs++;
            s->data = (char *) realloc(s->data, STARTSIZE*reallocs*sizeof(char));
        }
        c = buffered_fgetc(temp_seq_buffer, &idx, &r, query);
    }
    s->data[length]='\0';

    free(temp_seq_buffer);

    return length;
}
           
char buffered_fgetc(char *buffer, uint64_t *pos, uint64_t *read, FILE *f) {
    if (*pos >= READBUF) {
        *pos = 0;
        memset(buffer, 0, READBUF);
        *read = fread(buffer, 1, READBUF, f);
    }
    *pos = *pos + 1;
    return buffer[*pos-1];
}


void NW(char * X, uint64_t Xstart, uint64_t Xend, char * Y, uint64_t Ystart, uint64_t Yend, int iGap, int eGap, struct cell ** table, struct positioned_cell * mc, int show){
    
    uint64_t i,j;
    int64_t scoreDiagonal,scoreLeft,scoreRight,score;

    
    int percentage=0;
    
    struct positioned_cell mf;
    mf.score = INT64_MIN;
    

    //First row. iCounter serves as counter from zero
    //printf("..0%%");
    for(i=0;i<Yend;i++){
        table[0][i].score = pam_dna_score(X[0], Y[i]);
        //table[Xstart][i].xfrom = Xstart;
        //table[Xstart][i].yfrom = i;
        //Set every column max
        mc[i].score = table[0][i].score;
        mc[i].xpos = 0;
        mc[i].ypos = i;

    }
    
    //Set row max
    mf.score = table[0][0].score;
    mf.xpos = 0;
    mf.ypos = 0;

    //Go through full matrix
    for(i=1;i<Xend;i++){
        //Fill first rowcell
        if(show==1 && ((double)i/(double)Xend)*100 > percentage+1){
            percentage = ((double)i/(double)Xend)*100;
            printf("..%d%%\n", percentage);
        }
        


        table[i][0].score = pam_dna_score(X[i],Y[0]);
        mf.score = table[i][0].score;
        mf.xpos = i;
        mf.ypos = 0;

        for(j=1;j<Yend;j++){

            //Check if max in row has changed
            if(j > 1 && mf.score <= table[i][j-2].score){
                mf.score = table[i][j-2].score;
                mf.xpos = i;
                mf.ypos = j-2;
            }
            
            score = pam_dna_score(X[i], Y[j]);
            scoreDiagonal = table[i][j-1].score + score;

            if(j>1){
                scoreLeft = mf.score + iGap + (j - (mf.ypos+1))*eGap + score;
                }else{
                    scoreLeft = INT_MIN;
                }
                
            if(i>1){
                scoreRight = mc[j-1].score + iGap + (i - (mc[j-1].xpos+1))*eGap + score;
                }else{
                    scoreRight = INT_MIN;
                }
            
            //Choose maximum
            
            if(scoreDiagonal >= max(scoreLeft, scoreRight)){
                //Diagonal
                table[i][j].score = scoreDiagonal;
                table[i][j].xfrom = i-1;
                table[i][j].xfrom = j-1;
                                
            }else if(scoreRight >= scoreLeft){
                //Gap in genome
                table[i][j].score = scoreRight;
                table[i][j].xfrom = mc[j-1].xpos;
                table[i][j].yfrom = mc[j-1].ypos;
                
            }else{
                //Gap in read
                table[i][j].score = scoreLeft;
                table[i][j].xfrom = mf.xpos;
                table[i][j].yfrom = mf.ypos;
            }

            //check if column max has changed
            if(i > 1 && j > 1 && table[i-2][j-1].score > mc[j-1].score){
                mc[j-1].score = table[i-2][j-1].score;
                mc[j-1].xpos = i-2;
                mc[j-1].ypos = j-1;
            }
        }
    }

    fprintf(stdout, "[INFO] Table completed\n");

}
