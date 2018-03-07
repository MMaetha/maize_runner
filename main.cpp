#include <stdio.h>
#include <stdlib.h>
#include "fasta.h"
#include <map>
#include <cstring>
#include <iostream>

using namespace std;

FILE* outputFile;

multimap<string,int> kmap;
int k = 5;

void build_map(char* seq)
{
    int len = strlen(seq);
    for(int i=0; i<len-k+1; i++) {
        char tmp = seq[i+k];
        seq[i+k] = '\0';

        string st = string(seq+i);
        //cout << i << "--" << st << "--" << endl;
        kmap.insert(make_pair(st,i));
        seq[i+k] = tmp;
       /* if(i % 10000 == 0) {
            printf("%d\n",i);
        } */
        //free(kgram);

        if(i > 100000)
        	break;
    }
    // printf("Complete\n");
}

void search_map(char* seq)
{
    int len = strlen(seq);
    for(int i=0; i<len-k+1; i++) {
        char tmp = seq[i+k];
        seq[i+k] = '\0';
        //char* kgram = strdup(seq+i);

        string kgram = string(seq+i);

        //fprintf(outputFile,"%s\n",kgram);
        for(auto it = kmap.equal_range(kgram).first; it!=kmap.equal_range(kgram).second; ++it){
            fprintf(outputFile,"%d\n", (*it).second-i);
        }
        //fprintf(outputFile,"count na ja : %lu \n",kmap.count(kgram));
        seq[i+k] = tmp;
        //free(kgram);
    }
    //fprintf(outputFile,"Complete\n");
    //fprintf(outputFile,"---------------------------------\n");
}

void upup(char* st)
{
    char* p=st;
    while(*p != '\0') {
        if((*p >= 'a') && (*p <= 'z')) {
            *p += 'A' - 'a';
        }
        p++;
    }
}

int main(int argc, char **argv) {
    if(argc != 3){
        printf("Invalid parameters\n");
        exit(1);
    }
    else printf("bitch!\n");

    outputFile = fopen(argv[2],"a");
    if(outputFile == NULL){
        printf("Cannot open output file.\n");
        exit(1);
    }
    else printf("bitch!\n");

    FASTAFILE *ffp;
    char *seq;
    char *name;
    int L;
    int c;
    /* argv[1] is the name of a FASTA file */
    ffp = OpenFASTA("sample_seq\\IRGSP-1.0_genome.fasta");
    c = 1;
    int targetCms = atoi(argv[1]);
    printf("targetCms = %d\n", targetCms);
    while (ReadFASTA(ffp, &seq, &name, &L)) {
        //printf("%d\n",c);
        //printf("name: %s\n", name);
        //printf("size: %d\n", L);
        //printf("%s\n\n", seq);
        printf("Read chromosome file : c = %d\n", c);
        if(c==targetCms){
            upup(seq);
            build_map(seq);
            printf("--- Read chromosome file : c = %d complete ---\n", c);
        }

        free(seq);
        free(name);
        //getchar();
        c++;
        //break;
    }
    CloseFASTA(ffp);

    //for(auto it = kmap.begin(); it != kmap.end(); it++) {
    //    cout << it->first << it->second << endl;
    //}

    //return 0;

    ffp = OpenFASTA("sample_seq\\TAIR9_TE.fas");
    c = 1;
    while (ReadFASTA(ffp, &seq, &name, &L)) {
        //printf("%d\n",c);
        //printf("name: %s\n", name);
        //printf("size: %d\n", L);
        //printf("%s\n\n", seq);

        // write file
       /* fprintf(outputFile,"%d\n",c);
        fprintf(outputFile,"name: %s\n", name);
        fprintf(outputFile,"size: %d\n", L);
        */
       //build_map(seq);
        upup(seq);
        search_map(seq);
        free(seq);
        free(name);

        if(c>=1) break;
            printf("Search Map : c = %d\n", c);
        //getchar();
        c++;
    }
    CloseFASTA(ffp);

    fclose(outputFile);

    exit(0);
}
