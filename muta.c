//
//  muta.c
//
//
//  Copyright © 2014 Andres Breton
//
//  Mutant finder: looks and reports for gene knockouts on mutant candidate sequences
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


int hits(int *a, int *b, int *c, int *d);


int main(int argc, char *argv[]) {
    
    if (argc < 3) {
        fprintf(stderr, "\n[ERROR]: No input provided. Please provide the reference file followed by the query file\n\n.");
        printf("\nThis program is designed to help outmated BLAST execution of query sequences.\n");
        printf("\nSpecifically, it's intended for use with deletion mutants. This program will help you find whether your deletion fragment is present in your mutant candidate sequence.\n");
        printf("\nExecution Example:\n$ muta organism_sequence.fasta mutant_sequence.fasta\n\n");
        return 0;
    }
    
    
    //COMMAND VARIABLES
    
    char BLAST[250];
    char *pwd;
    char *outfmt;
    char blastout[50];
    
    pwd = getenv("PWD");
    outfmt = "\"6 qseqid sseqid pident sstart send length qlen evalue\"";
    
    sprintf(blastout, "%s.blastout", argv[2]);
    sprintf(BLAST, "blastn -db %s -query %s -outfmt %s -out %s", argv[1], argv[2], outfmt, blastout);
    
    
    //EXECUTE BLAST & PRINT RESULTS
    
    system(BLAST);
    
    
    if (system(BLAST) == 0) {
        printf("\n");
        printf("\n**********************************************************\n");
        printf("* BLAST was successful!\n");
        printf("* Output file < %s > created", blastout);
        printf("\n**********************************************************\n");
        printf("\n");
        
        printf("These are the BLAST hits for your sequence...\n\n");
        
        
        //EXTRACT BLAST HIT VALUES
        
        FILE *blastf;
        int a, b, c, d;
        char CUT[150];
        
        sprintf(CUT, "for i in 4 5; do cut -f$i %s ;done", blastout); //4 and 5 are tab locations in file containing hits
    
        system(CUT);
        
        blastf = popen(CUT, "r");
        if (blastf == NULL) {
            printf("Could not open/find BLAST output file %s", blastout);
            return 0;
        }
            fscanf(blastf, "%d %d %d %d", &a, &c, &b, &d);
            pclose(blastf);
       
        
        //EXPLAIN RESULTS
        //printf("\n\nNOTE: If 2 fragments were produced by the BLAST search:\n* The first two numbers (1 and 2) are the start and end, respectively, of your RFR (right-flanking region).\n* The remaining two numbers (3 and 4) are that of the start and end, respectively, of your LFR (left-flanking region)\n\n");
        
        hits(&a, &b, &c, &d);
        
    } else {
        printf("\n\n**********************************************************\n");
        printf("* Sorry, BLAST was UNSUCCESSFUL :( \n");
        printf("* Please re-run program");
        printf("\n**********************************************************\n\n");
      
    }
    
    return 0;
}



int hits(int *a, int *b, int *c, int *d){
    
    int START, END, ps, GSIZE, DELSIZE;
    
    
    if (*a == 0) {
        printf("**********************************************************\n");
        printf("* SORRY: Your mutant did NOT produce any BLAST hits");
        printf("\n**********************************************************\n\n");
        
        return 0;
    }else {
    
    printf("\nPlease enter the starting nucleotide position of your deleted fragment in its genome (ex. 12345): ");
    scanf("%d", & START);
    printf("\nPlease enter the ending nucleotide position of your deleted fragment in its genome (ex. 12589): ");
    scanf("%d", & END);
    //printf("\n\nGene sequence in a mutant can be included in the sequence at times, yet still account as COMPLETE gene deletion.");
    //printf("\n\tYou can select the percentage of the gene's sequence to be included...");
    //printf("\n\t0 = complete gene deletion (No gene sequence allowed)");
    //printf("\nPlease enter percentage (%%) [ex. 15]: ");
    //scanf("%d", & ps);
        
    }
    
    
    GSIZE = END - START +1;
    printf("\n\nGene size: %d\n\n", GSIZE);
    
//   printf("Values= %d %d %d %d", *a,*b,*c,*d); //CHECK
    
    
    //ACCOUNT SINGLE FRAGMENT HITS...STATE IT
    
    while (*b == 0){
        printf("\n**NOTE**: Your mutant contains only 1 of the fragments (LFR or RFR)!");
        if (*a > 0) {
            if (*a <= START) {
                if (*c <= START){
                    printf("\n\tYou have a LFR at %d-%d\n\n", *a, *c);
                    //break;
                }
                if (*c >= START) {
                    printf("\nMutant contains gene sequence!\n\n");
                    //break;
                }
            }
            
            if (*a >= START && *a <= END) {
                printf("\nMutant contains gene sequence!\n\n");
                //break;
            }
            
            if (*a && *c >= END) {
                printf("\n\tYour RFR is at %d-%d\n\n", *a, *c);
                //break;
            }
            
        }
        return 0;
    }
    
    
    //TAKE CARE OF POSITIONING, COMPARE HITS TO GENE LOCATION, AND REPORT RESULTS OF MUTANT/KNOCKOUT GENE
    
    if (*c < *b) {
        if (*a < *b) {
            if (*c < *d) {
                if (*d <= START || *a >= END) {
                    printf("CONGRATULATIONS! You have created a mutant with complete gene deletion!\n\n");
                } else {
                    printf("\n**********************************************************\n");
                    printf("Sorry: Your mutant does not have the gene deletion desired.\n\n");
                    }
            //(*C > *D)
            } else {
                if (*c <= START || *a >= END) {
                    printf("CONGRATULATIONS! You have created a mutant with complete gene deletion!\n\n");
                } else {
                    printf("\n**********************************************************\n");
                    printf("Sorry: Your mutant does not have the gene deletion desired.\n\n");
                }
            }
        
        //(*A > *B)
        } else {
            if (*c < *d) {
                if (*d <= START || *b >= END) {
                    printf("CONGRATULATIONS! You have created a mutant with complete gene deletion!\n\n");
                } else {
                    printf("\n**********************************************************\n");
                    printf("Sorry: Your mutant does not have the gene deletion desired.\n\n");
                }
            //(*C > *D)
            } else {
                if (*c <= START || *b >= END) {
                    printf("CONGRATULATIONS! You have created a mutant with complete gene deletion!\n\n");
                } else {
                    printf("\n**********************************************************\n");
                    printf("Sorry: Your mutant does not have the gene deletion desired.\n\n");
                }
            }
        }
        
     //(*C > *B)
    }else {
        if (*a < *b) {
            if (*c < *d) {
                if (*b <= START || *c >= END) {
                    printf("CONGRATULATIONS! You have created a mutant with complete gene deletion!\n\n");
                } else {
                    printf("\n**********************************************************\n");
                    printf("Sorry: Your mutant does not have the gene deletion desired.\n\n");
                }
            //(*C > *D)
            } else {
                if (*b <= START || *d >= END) {
                    printf("CONGRATULATIONS! You have created a mutant with complete gene deletion!\n\n");
                } else {
                    printf("\n**********************************************************\n");
                    printf("Sorry: Your mutant does not have the gene deletion desired.\n\n");
                }
            }
            
        //(*A > *B)
        } else {
            if (*c < *d) {
                if (*a <= START || *c >= END) {
                    printf("CONGRATULATIONS! You have created a mutant with complete gene deletion!\n\n");
                } else {
                    printf("\n**********************************************************\n");
                    printf("Sorry: Your mutant does not have the gene deletion desired.\n\n");
                }
            //(*C > *D)
            } else {
                if (*a <= START || *d >= END) {
                    printf("CONGRATULATIONS! You have created a mutant with complete gene deletion!\n\n");
                } else {
                    printf("\n**********************************************************\n");
                    printf("Sorry: Your mutant does not have the gene deletion desired.\n\n");
                }
            }
        }
  
    }
    
    return 0;
    
}


