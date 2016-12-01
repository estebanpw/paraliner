#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <inttypes.h>
#include <limits.h>
#include <ctype.h>
#include "structs.h"
#include "common.h"




inline int64_t pam_dna_score(char c1, char c2){
    if(c1==c2) return MATCH;
    return (-MATCH+2);
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


struct positioned_cell NW(char * X, uint64_t Xstart, uint64_t Xend, char * Y, uint64_t Ystart, uint64_t Yend, int64_t iGap, int64_t eGap, struct cell ** table, struct positioned_cell * mc, int show){
    
    uint64_t i,j;
    int64_t scoreDiagonal,scoreLeft,scoreRight,score;

    struct positioned_cell bc;
    bc.score = INT64_MIN;
    
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
                mf.score = table[i-1][j-2].score;
                mf.xpos = i-1;
                mf.ypos = j-2;
            }
            
            score = pam_dna_score(X[i], Y[j]);
            scoreDiagonal = table[i-1][j-1].score + score;

            if(j>1){
                scoreLeft = mf.score + iGap + (j - (mf.ypos+1))*eGap + score;
                }else{
                    scoreLeft = INT64_MIN;
                }
                
            if(i>1){
                scoreRight = mc[j-1].score + iGap + (i - (mc[j-1].xpos+1))*eGap + score;
                }else{
                    scoreRight = INT64_MIN;
                }
            
            //Choose maximum
            //printf("Score DIAG: %"PRId64"; LEFT: %"PRId64"; RIGHT: %"PRId64"\n", scoreDiagonal, scoreLeft, scoreRight);
            if(scoreDiagonal >= scoreLeft && scoreDiagonal >= scoreRight){
                //Diagonal
                table[i][j].score = scoreDiagonal;
                table[i][j].xfrom = i-1;
                table[i][j].yfrom = j-1;
                                
            }else if(scoreRight > scoreLeft){
                table[i][j].score = scoreRight;
                table[i][j].xfrom = mc[j-1].xpos;
                table[i][j].yfrom = mc[j-1].ypos;
                
            }else{
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
	    if(i == Xend-1 || j == Yend-1){
		//Check for best cell
		if(table[i][j].score >= bc.score){ bc.score = table[i][j].score; bc.xpos = i; bc.ypos = j; }
	    }
        }
    }
    
    fprintf(stdout, "[INFO] Table completed\n");
    return bc;
}



void backtrackingNW(char * X, uint64_t Xstart, uint64_t Xend, char * Y, uint64_t Ystart, uint64_t Yend, struct cell ** table, char * rec_X, char * rec_Y, struct positioned_cell * bc, uint64_t * ret_head_x, uint64_t * ret_head_y){
    uint64_t curr_x, curr_y, prev_x, prev_y, head_x, head_y;
    int64_t k;
    head_x = 2*max(Xend, Yend);
    head_y = 2*max(Xend, Yend);
    curr_x = bc->xpos;
    curr_y = bc->ypos;
    prev_x = curr_x;
    prev_y = curr_y;
   
    
    for(k=Xend-1; k>curr_x; k--) rec_X[head_x--] = '-';
    for(k=Yend-1; k>curr_y; k--) rec_Y[head_y--] = '-';


    //PONER HEAD_X y HEAD_Y con tamaño doble de cinta!!!!!!!
    while(curr_x > 0 && curr_y > 0){
        curr_x = table[prev_x][prev_y].xfrom;
        curr_y = table[prev_x][prev_y].yfrom;
    	//printf("(%"PRIu64", %"PRIu64") ::: ", curr_x, curr_y);
    	//printf("(%"PRIu64", %"PRIu64") ::: ", prev_x, prev_y);
	//getchar();

        if((curr_x == (prev_x - 1)) && (curr_y == (prev_y -1))){
            //Diagonal case
	    //printf("DIAG\n");
            rec_X[head_x--] = X[prev_x];
            rec_Y[head_y--] = Y[prev_y];
        }else if((prev_x - curr_x) > (prev_y - curr_y)){
            //Gap in X
	    //printf("Gap X\n");
            for(k=prev_x;k>=curr_x;k--){
                rec_Y[head_y--] = '-';
		rec_X[head_x--] = X[k];
            }
        }else{
            //Gap in Y
	    //printf("GAP Y\n");
            for(k=prev_y;k>=curr_y;k--){
                rec_X[head_x--] = '-';
		rec_Y[head_y--] = Y[k];
            }
        }
	prev_x = curr_x;
	prev_y = curr_y;
    }
    //printf("curr: %"PRIu64", %"PRIu64"\n", curr_x, curr_y);
    uint64_t huecos_x = 0, huecos_y = 0;
    for(k=(int64_t)curr_x-1; k>=0; k--){ rec_X[head_x--] = '-'; huecos_x++;}
    for(k=(int64_t)curr_y-1; k>=0; k--){ rec_Y[head_y--] = '-'; huecos_y++;}
    
    if(huecos_x >= huecos_y){
	while(huecos_x > 0) {rec_Y[head_y--] = ' '; huecos_x--;}
    }else{
	while(huecos_y > 0) {rec_X[head_x--] = ' '; huecos_y--;}
    }

    *ret_head_x = head_x;
    *ret_head_y = head_y;
}
