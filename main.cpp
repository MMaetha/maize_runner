#include <stdio.h>
#include <stdlib.h>
#include "fasta.h"
#include <map>
#include <cstring>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>

using namespace std;

FILE* outputFile;

multimap<string,int> kmap;
int k = 5;
int peak=-100,place=-100;


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
       	// if(i % 10000 == 0) {
        //    printf("%d\n",i);
        //} 
        //free(kgram);

        if(i > 100000)
        	break;
    }
    // printf("Complete\n");
}

int min(int x, int y, int z){
   	return min(min(x, y), z);
}

int editDist(string str1 , string str2 , int m ,int n)
{
    if (m == 0) return n;

    if (n == 0) return m;

    if (str1[m-1] == str2[n-1])
        return editDist(str1, str2, m-1, n-1);
    
    return 1 + min ( editDist(str1,  str2, m, n-1),    // Insert
                     editDist(str1,  str2, m-1, n),   // Remove
                     editDist(str1,  str2, m-1, n-1) // Replace
                   );
}

void search_map(char* seq)
{
    int len = strlen(seq);
    vector<int> v;
    for(int i=0; i<len-k+1; i++) {
        char tmp = seq[i+k];
        seq[i+k] = '\0';
        //char* kgram = strdup(seq+i);

        string kgram = string(seq+i);

        //fprintf(outputFile,"%s\n",kgram);
        for(auto it = kmap.equal_range(kgram).first; it!=kmap.equal_range(kgram).second; ++it){
            // printf("%d\n", (*it).second-i);
            v.push_back(it->second - i);
        }
        
        //fprintf(outputFile,"count na ja : %lu \n",kmap.count(kgram));
        seq[i+k] = tmp;
        //free(kgram);
    }
    printf("Mapping completed.\n");
    sort(v.begin(), v.end());
    vector<int>::iterator low,up;
    vector<int> highRank,highFreq;
    int freq;
    int percentile = 85;
    for(int i=0; i<v.size(); i++) {
    		if(v[i]<0){
    			continue;
    		}
    		else{
    			low = lower_bound(v.begin(),v.end(),v[i]);
    			up = upper_bound(v.begin(),v.end(),v[i]);
    			freq = up-low;
    			if(v[i+1]!=v[i]){
     				fprintf(outputFile, "%d , %d\n",v[i],freq);
     				}
     			}
     		if(freq>peak){
     			peak = freq;
     			place = i;
   	   		}       
    	}
    printf("Maximum frequency is %d at %d \n",peak,place);

    for(int i=0; i<v.size(); i++) {
    		if(v[i]<0){
    			continue;
    		}
    		else{
    			low = lower_bound(v.begin(),v.end(),v[i]);
    			up = upper_bound(v.begin(),v.end(),v[i]);
    			freq = up-low;
     			}
     		if(freq>percentile*peak/100){
     			highRank.push_back(v[i]);
     			highFreq.push_back(freq);
   	   		}       
    	}

    printf("High frequency place (> %d %)\n",percentile);
    for(int i=0; i<highRank.size(); i++) {
    	if(highRank[i+1]!=highRank[i])
    	printf("At %d has %d freq.\n", highRank[i], highFreq[i]);
     }

    printf("Sorting completed.\n");
	//printf("Similarity score = %.2f",)
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

    outputFile = fopen(argv[2],"a");
    if(outputFile == NULL){
          printf("Cannot open output file.\n");
        exit(1);
    }

    FASTAFILE *ffp;
    char *seq;
    char *name;
    int L;
    int c;
    string chr,TE;
    /* argv[1] is the name of a FASTA file */
    printf("opening FASTA file.\n");
    ffp = OpenFASTA("sample_seq/IRGSP-1.0_genome.fasta");
    printf("openFASTA completed\n");
    c = 1;
    int targetCms = atoi(argv[1]);
    printf("targetCms = %d\n", targetCms);
    while (ReadFASTA(ffp, &seq, &name, &L)) {
        //printf("%d\n",c);
        //printf("name: %s\n", name);
        //printf("size: %d\n", L);
        //printf("%s\n\n", seq);
        //printf("Read chromosome file : c = %d\n", c);
        if(c==targetCms){
        	chr = seq;
            upup(seq);
            build_map(seq);
        }

        free(seq);
        free(name);
        //getchar();
        c++;
        //break;
    }
    printf("--- Read chromosome file : c = %d completed ---\n",targetCms);
    printf("Reading file completed.\n");
    CloseFASTA(ffp);

    //for(auto it = kmap.begin(); it != kmap.end(); it++) {
    //    cout << it->first << it->second << endl;
    //}

    //return 0;

    ffp = OpenFASTA("sample_seq/TAIR9_TE.fas");
    c = 2;
    while (ReadFASTA(ffp, &seq, &name, &L)) {
        //printf("%d\n",c);
        //printf("name: %s\n", name);
        //printf("size: %d\n", L);
        //printf("%s\n\n", seq);

        // write file
       	// fprintf(outputFile,"%d\n",c);
        // fprintf(outputFile,"name: %s\n", name);
        // fprintf(outputFile,"size: %d\n", L);
        //build_map(seq);
        upup(seq);
        search_map(seq);
        TE = seq;
        free(seq);
        free(name);

        if(c>=2) break;
            printf("Search Map : c = %d\n", c);
        //getchar();
        c++;
    }
    
    CloseFASTA(ffp);
    printf("%d",peak);
    string chr_short = chr.substr(place-200,place+200);
    double matchScore = editDist(chr_short,TE,chr_short.length(),TE.length());
  	printf("Max score is %.2f",matchScore);
    fclose(outputFile);

    printf("Completed.\n");

    exit(0);
}
