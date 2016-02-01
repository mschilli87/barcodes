/**********************
* general information *
**********************/

/*
 * file:        remap_barcodes.c
 * created:     2016-01-11
 * last update: 2016-02-01
 * author(s):   Marcel Schilling <marcel.schilling@mdc-berlin.de>
 * purpose:     re-map cell barcodes for drop-seq data
 */


/*************************************
* change log (reverse chronological) *
*************************************/

/*
 * 2016-02-01:  added collapsing of barcodes to use that differ in one position only amogst each
 *              other
 * 2016-01-29:  replaced TOP flag with negative index to block top barcodes in hash table
 *              switched values of AMBIGUOUS & TOP flags
 *              switched top barcode indexing from 0-based to 1-based to use the now-free 0 index
 *              to block top barcodes in hash (AMBIGUOUS vs. TOP)
 *              named magic numbers (AMBIGUOUS & FILE_END) to improve clarity
 *              fixed typos in comments ('[un]ambigous[ly]' --> '[un]ambiguous[ly]')
 *              replaced ssize_t by size_t for positive constant (barcode length) & loop indices
 *              (see http://stackoverflow.com/q/15739490 for detailed discussion)
 * 2016-01-13:  added generation of barcodes with 1 deletion
 *              renamed variables ('barcode' --> 'target', 'line' --> 'barcode',
 *              'barcode_len' --> 'correct_barcode_length', 'line_length' --> 'barcode_length') for
 *              clarity / fixed comments & line wrapping / removed trailing whitespace
 *              doubled hash table size to avoid collisions
 *              added optional 3rd command line argument to specify number of barcode to use / added
 *              parameter checking
 *              switched to test for skipping of mutation already before modifying the line buffer
 *              resolved ambiguity for barcodes to use with distance 1 / fixed typos in comments
 *              fixed store() function to store keys by value, not reference (see
 *              http://stackoverflow.com/q/34752934/2451238)
 * 2016-01-12:  replaced nested 'if' with 'continue' statement / standardised library comments /
 *              fixed typo in comment
 *              added barcodes to use to hash table to avoid 1-missmatch collisions / fixed comments
 *              fixed off-by-one error in check for too many barcodes to use in input file
 *              adjusted for-loop to allow -O2 without causing undefined behavior / fixed typo in
 *              comments ('satify' --> 'satisfy')
 *              fixed definition of number of barcodes to use constant / fixed comment alignment
 *              fixed const string definition for command line arguments
 *              adjusted RPAA struct initialization to satisfy -Wpedantic
 *              replaced C++ style comments by ANSI C style comments
 *              replaced test output with reading & re-mapping of the barcodes to map to the ones to
 *              use, if possible
 * 2016-01-11:  replaced example code by reading barcodes to use and generating mappings for all
 *              1-mutation-barcodes (incl. test output with mappings for of all N1A mutations)
 *              labeled Rosettacode POSIX associative arrays (RPAA) lines individually / added
 *              section headers
 *              initial version (POSIX associative array example taken from
 *              http://rosettacode.org/wiki/Associative_arrays/Creation/C)
 */


/*******************
* copyright notice *
*******************/

/*
 * All lines ending with an 'RPAA' (Rosettacode POSIX associative arrays) comment were taken from
 * the POSIX associative array example released under GNU Free Documentation License 1.2
 * (http://www.gnu.org/licenses/fdl-1.2.html) at
 * from http://rosettacode.org/wiki/Associative_arrays/Creation/C.
 * The code was copied from there on 2016-01-11 (site version from 2012-09-29, 23:48).
 */


/************
* libraries *
************/

#include <inttypes.h> /* intptr_t                                      */
#include <search.h>   /* hcreate(), hsearch()                          */
#include <stdio.h>    /* perror(), fopen(), getline(), printf()        */
#include <string.h>   /* strdup(), strcpy()                            */
#include <stdlib.h>   /* exit(), EXIT_FAILURE, EXIT_SUCCESS, strtoul() */
#include <errno.h>    /* errno                                         */


/*************
* data types *
*************/

/* entry in hit list (see below)*/
typedef struct hit_list_entry{
  /*
   * ID of one barcode to use that has hamming distance 1 to the barcode to use corresponding to the
   * hit list the entry belongs to
   */
  size_t hit_id;

  /* position in which the two barcodes to use differ (hamming distance 1) */
  size_t pos_mismatch;

  /* pointer to the next entry in hit list the entry belongs to */
  struct hit_list_entry* next;
} hit_list_entry_t;

/* hit list with barcodes to use that have hamming distance 1 to a certain barcode to use */
typedef struct hit_list{

  /* pointer to the first entry in the list */
  hit_list_entry_t* head;

  /* length of the list */
  size_t length;
} hit_list_t;


/************
* functions *
************/

void exit_with_error(const char* error_message){
  perror(error_message);
  exit(EXIT_FAILURE);
}

/*                                                                                         * RPAA *
 * Must hcreate() the hash table before calling fetch() or store().                        * RPAA *
 *                                                                                         * RPAA *
 * Because p->data is a pointer, fetch() and store() cast between                          * RPAA *
 * void * and intptr_t.                                                                    * RPAA *
 */                                                                                       /* RPAA */

/* Fetch value from the hash table. */                                                    /* RPAA */
int                                                                                       /* RPAA */
fetch(const char *key, intptr_t *value)                                                   /* RPAA */
{                                                                                         /* RPAA */
  /* previously: */ /* ENTRY e = {key: (char *)key}, *p;  */                              /* RPAA */
  ENTRY e,*p;       /* adjusted from RPAA version (see above) to satisfy -Wpedantic */
  e.key=(char*)key; /* adjusted from RPAA version (see above) to satisfy -Wpedantic */
  p = hsearch(e, FIND);                                                                   /* RPAA */
  if (p) {                                                                                /* RPAA */
    *value = (intptr_t)p->data;                                                           /* RPAA */
    return 1;                                                                             /* RPAA */
  } else                                                                                  /* RPAA */
    return 0;                                                                             /* RPAA */
}                                                                                         /* RPAA */

/* Store key-value pair into the hash table. */                                           /* RPAA */
void                                                                                      /* RPAA */
store(const char *key, intptr_t value)                                                    /* RPAA */
{                                                                                         /* RPAA */
  /*                                                                                       * RPAA *
   * hsearch() may insert a new entry or find an existing entry                            * RPAA *
   * with the same key. hsearch() ignores e.data if it finds an                            * RPAA *
   * existing entry. We must call hsearch(), then set p->data.                             * RPAA *
   */                                                                                     /* RPAA */
  /* previously: */  /* ENTRY e = {key: (char *)key}, *p;  */                             /* RPAA */
  ENTRY e,*p;        /* adjusted from RPAA version (see above) to satisfy -Wpedantic   */
  e.key=strdup(key); /* see http://stackoverflow.com/a/34753734/2451238 for discussion */
  p = hsearch(e, ENTER);                                                                  /* RPAA */
  if (p == NULL)                                                                          /* RPAA */
    exit_with_error("hsearch");
  p->data = (void *)value;                                                                /* RPAA */
}                                                                                         /* RPAA */

/* save 1-based indices in 0-based array (index 0 will represent ambiguity) */
char *save(char *const *array, const size_t index, const char* value){
  return strcpy(array[index-1],value);
}

/* get 1-based indices from 0-based array (index 0 will represent ambiguity) */
const char *get(char *const *const array, const size_t index){return array[index-1];}

/* create an empty hit list */
hit_list_t* empty_hit_list(){

  /* create a new hit list object (malloc'ed memory must be free'ed later) */
  hit_list_t* the_list=(hit_list_t*)malloc(sizeof(hit_list_t));

  /* empty list is headless */
  the_list->head=NULL;

  /* empty list has no entries */
  the_list->length=0;

  /* return (pointer to) newly generated empty list */
  return(the_list);
}

/* insert a barcode ID into a hit list */
void insert(hit_list_t* the_list,const size_t the_hit_id,const size_t the_pos_mismatch){

  /* create a new hit_list_entry object (malloc'ed memory must be free'ed later) */
  hit_list_entry_t* the_entry=(hit_list_entry_t*)malloc(sizeof(hit_list_entry_t));

  /* assign given barcode ID to new hit list entry */
  the_entry->hit_id=the_hit_id;

  /* assign given mismatch position to new hit list entry */
  the_entry->pos_mismatch=the_pos_mismatch;

  /* new entry will be inserted before previous list head */
  the_entry->next=the_list->head;

  /* new entry becomes new list head */
  the_list->head=the_entry;

  /* list length increased by 1 */
  ++(the_list->length);
}

/* delete a hit list */
void delete_hit_list(hit_list_t* the_list){

  /* only delete list if it exists */
  if(the_list!=NULL){

    /* pointer to entry to delete next from the list */
    hit_list_entry_t* the_entry=NULL;

    /* as long as the list still has entries */
    for(the_entry=the_list->head;the_entry!=NULL;the_entry=the_list->head){

      /* separate the first entry from the list */
      the_list->head=the_entry->next;

      /* delete the separate entry by free'ing its memory */
      free(the_entry);

      /* un-set pointer to free'ed memory */
      the_entry=NULL;
    }

    /* delete the empty list by free'ing its memory */
    free(the_list);

    /* un-set pointer to free'ed memory */
    the_list=NULL;

  } /* existing list deleted */
}

/* clear a hit list */
void clear_hit_list(hit_list_t* the_list){

  /* delete old hit list */
  delete_hit_list(the_list);

  /* assign new empty hit list */
  the_list=empty_hit_list();
}


/************
* constants *
************/

/* number of nucleotides to be considered for the correction of a cell barcode */
/*
 * Note: The following declaration cannot be used to use n_nucleotides as size for the
 * nucleotides array constant (see below):
 * static const ssize_t n_nucleotides=4; // warning: excess elements in array initializer
 * See http://stackoverflow.com/a/1674459 for a detailed discussion.
 */
enum{n_nucleotides=4};

/* nucleotides to be considered for the correction of a cell barcode */
static const char nucleotides [n_nucleotides]={'A','C','G','T'};

/* default number of barcodes to use (i.e. to map the remaining ones to, if possible) */
/*
 * Note: The following declaration cannot be used to use n_barcodes_use_default to initialize the
 * n_barcodes_use parameter:
 * static const size_t n_barcodes_use_default=1000; // error: initializer element is not constant
 * See http://stackoverflow.com/a/1674459 for a detailed discussion.
 */
enum{n_barcodes_use_default=1000};

/* length of barcodes */
static const size_t correct_barcode_length=12;

/* size factor to scale the hash table with to avoid collisions */
static const unsigned short hash_table_size_factor=2;

/* line length returned by getline on file end */
static const ssize_t FILE_END=-1;

/* barcode index to use for ambiguously mappable barcodes */
static const ssize_t AMBIGUOUS=0;

/* character used to represent 'any nucleotide' */
static const char any_nucleotide='N';


/************
* variables *
************/

/* take command line arguments & return exit status */
int main(int argc,char** argv){

  /* file to read barcodes to use (i.e. to map the remaining ones to, if possible) from */
  const char* file_barcodes_use=NULL;

  /* file to read barcodes to map to the ones to use, if possible, from */
  const char* file_barcodes_remap=NULL;

  /* default number of barcodes to use  */
  size_t n_barcodes_use=n_barcodes_use_default;

  /* file to read from */
  FILE* file;

  /* number of barcodes read from file */
  size_t barcodes_read=0;

  /* barcode read from file */
  char* barcode=NULL;

  /* size of buffer used to read line from file */
  size_t buff_size=0;

  /* length of barcode read from file */
  ssize_t barcode_length=0;

  /* iterators used for for-loops */
  size_t i,j=0;

  /* character array to store barcodes to use */
  char** barcodes_use=NULL;

  /* character used to iterate over strings */
  char orig='\0';

  /* index of barcode to use retrieved from hash table */
  intptr_t target=NULL;

  /* array to count number of top barcodes having hamming distance 1 to each top hit */
  hit_list_t** top_hits=NULL;

  /* pointer to iterate over hit lists */
  hit_list_entry_t* top_hit=NULL;


/*************
* parameters *
*************/

  /* ensure at least the input file names were specified */
  if(argc<3) exit_with_error("missing parameter(s)");

  /* file to read barcodes to use (i.e. to map the remaining ones to, if possible) from */
  file_barcodes_use=argv[1];

  /* file to read barcodes to map to the ones to use, if possible, from */
  file_barcodes_remap=argv[2];

  /* optional: number of barcodes to use (i.e. to map the remaining ones to, if possible) */
  if(argc==4){
    char *p=NULL; errno=0;
    n_barcodes_use=strtoul(argv[3],&p,10);
    if(errno||*p!='\0') exit_with_error("invalid number");
  }

  /* ensure no extra parameters were specified */
  if(argc>4) exit_with_error("too many parameters");


/*****************
* initialization *
*****************/

  /* initialize hash table */
  if(!hcreate(hash_table_size_factor*                      /* hash collision prevention   */
              n_barcodes_use*                              /* for each barcode to use:    */
                (1+                                        /*   original barcode          */
                 correct_barcode_length*(n_nucleotides-1)+ /*   barcodes with 1 missmatch */
                 (correct_barcode_length-1)*n_nucleotides  /*   barcodes with 1 deletion  */
                )
             )
    ) exit_with_error("cannot initialize hash table");

  /* initialize barcode & top hits arrays (malloc'ed memory must be free'ed later) */
  barcodes_use=(char**)malloc(n_barcodes_use*sizeof(char*));
  top_hits=(hit_list_t**)malloc(n_barcodes_use*sizeof(hit_list_t*));
  for(i=0;i<n_barcodes_use;++i){
    barcodes_use[i]=(char*)malloc((correct_barcode_length+1)*sizeof(char));
    top_hits[i]=empty_hit_list();
  }


/***********************
* read barcodes to use *
***********************/

  /* open file with barcodes to use (i.e. to map the remaining ones to, if possible) for reading */
  if((file=fopen(file_barcodes_use,"r"))==NULL) exit_with_error(file_barcodes_use);

  /* read barcodes to use from file & keep track of their count */
  for(barcodes_read=0;(barcode_length=getline(&barcode,&buff_size,file))!=FILE_END;++barcodes_read){

    /* ensure input file does not contain more barcodes than would fit into the hash table */
    if(barcodes_read==n_barcodes_use) exit_with_error("too many barcodes to use in input list");

    /* remove newline from barcode */
    barcode[--barcode_length]='\0';

    /* ensure input file does not contain any barcode of wrong length */
    if(barcode_length!=correct_barcode_length) exit_with_error("wrong barcode length");

    /* store original barcode in array */
    save(barcodes_use,barcodes_read+1,barcode);

    /* block original barcode in hash table */
    store(barcode,-1*(int)(barcodes_read+1));


/*************************************
* generate barcodes with 1 missmatch *
*************************************/

    /* mutate each position in the barcode remembering original nucleotide */
    for(orig=barcode[i=0];i<barcode_length;orig=barcode[++i]){

      /* try all nucleotides as replacement */
      for(j=0;j<n_nucleotides;++j) {

        /* skip replacement with the original itself */
        if(nucleotides[j]==orig) continue;

        /* mutate current position to current nucleotide */
        barcode[i]=nucleotides[j];

        /* add new mapping if current mutated barcode is unassigned */
        if(!fetch(barcode,&target)) store(barcode,barcodes_read+1);

        /* count top hit if current mutated barcode matched unmutated barcode read previously*/
        else if((int)target<0){

          /* add target to top hits of current barcode */
          insert(top_hits[barcodes_read],-1*(int)target,i);

          /* add current barcode to top hits of target */
          insert(top_hits[-1*(int)target-1],barcodes_read+1,i);
        }

        /* otherwise remove ambiguous old mapping */
        else store(barcode,AMBIGUOUS);

      } /* all possible nucleotides exhausted for current position */

      /* restore current position before proceeding with the next position to mutate */
      barcode[i]=orig;

    } /* all possible positions exhausted for current barcode */


/************************************
* generate barcodes with 1 deletion *
************************************/

    /* delete 1st position in the barcode remembering original nucleotide */
    for(orig=barcode[i=0];i<barcode_length-1;++i) barcode[i]=barcode[i+1];

    /* delete each position in the barcode remembering original nucleotide */
    for(i=0;i<barcode_length-1;++i){

      /* try all nucleotides as replacement */
      for(j=0;j<n_nucleotides;++j) {

        /* add current nucleotide to the end (1st UMI nucleotide shifted in by deletion) */
        barcode[barcode_length-1]=nucleotides[j];

        /* add new mapping if current mutated barcode is unassigned */
        if(!fetch(barcode,&target)) store(barcode,barcodes_read+1);

        /* otherwise remove ambiguous old mapping */
        else if((int)target>0) store(barcode,AMBIGUOUS);

      } /* all possible nucleotides exhausted for current position */

      /* store next position at barcode end */
      barcode[barcode_length-1]=barcode[i+1];

      /* shift deletion 1 nucleotide downstream */
      barcode[i+1]=orig;

      /* remember original nucleotide at current position */
      orig=barcode[barcode_length-1];

    } /* all possible positions exhausted for current barcode */

    /* restore last position (not necessary since barcode will be overwritten in next iteration) */
    /* barcode[barcode_length-1]=orig; */

  } /* all possible lines exhausted from file with barcode reads to use */

  /* close file */
  fclose(file);

  /* ensure input file does not contain less barcodes than are supposed to be used */
  if(barcodes_read<n_barcodes_use) exit_with_error("too little barcodes to use in input list");


/********************************************
* re-map barcodes to use amongst each other *
********************************************/

  /* check all barcodes to use for potential collapsing with other barcodes to use */
  for(i=0;i<n_barcodes_use;++i)

    /* skip barcodes to use that have hamming distance >1 to all other barcodes to use */
    if(top_hits[i]->length!=0){

   /* compare hit length lists of all hits to number of hits */
      for(top_hit=top_hits[i]->head
         ;top_hit!=NULL&&top_hits[top_hit->hit_id-1]->length==top_hits[i]->length
         ;top_hit=top_hit->next){}

   /*
    * if all barcodes to use that have hamming distance 1 to the current barcode to use have the
    * same number of barcodes to use with hamming distance 1 to them as the current one, then
    * they all differ in the same postion and can be collapsed by setting this position to N in
    * all of those barcodes
    */
      if(top_hit==NULL){

     /* print mapping of current barcode to use to the collapsed one while remapping it */
        printf("%s\t",get(barcodes_use,i+1));
        barcodes_use[i][top_hits[i]->head->pos_mismatch]=any_nucleotide;
        printf("%s\n",get(barcodes_use,i+1));

     /* remap all barcodes to use with hamming distance 1 to the current one */
        for(top_hit=top_hits[i]->head;top_hit!=NULL;top_hit=top_hit->next){

       /* print mapping of current barcode to use to the collapsed one while remapping it */
          printf("%s\t",get(barcodes_use,top_hit->hit_id));
          barcodes_use[top_hit->hit_id-1][top_hit->pos_mismatch]=any_nucleotide;
          printf("%s\n",get(barcodes_use,top_hit->hit_id));

       /* clear corresponding hit lists to avoid double-processing */
          clear_hit_list(top_hits[top_hit->hit_id-1]);

        } /* all barcodes to use with hamming distance 1 to the current one remapped */

      } /* current collapsing done */

    } /* current barcode checked for potential collapsing */

  /* unset now un-used pointer iterator */
  top_hit=NULL;


/***************************************************************
* re-map barcodes to map to the ones to use to the ones to use *
***************************************************************/

  /* open file with barcodes to map the ones to use, if possible for reading */
  if((file=fopen(file_barcodes_remap,"r"))==NULL) exit_with_error(file_barcodes_remap);

  /* read barcodes to map to the ones to use, if possible, from file */
  while((barcode_length=getline(&barcode,&buff_size,file))!=FILE_END){

    /* remove newline from barcode */
    barcode[--barcode_length]='\0';

    /* ensure input file does not contain any barcode of wrong length */
    if(barcode_length!=correct_barcode_length) exit_with_error("wrong barcode length");

    /* print mapping if current barcode was unambiguously assigned to a barcode to use */
    if(fetch(barcode,&target)&&target!=AMBIGUOUS)
      printf("%s\t%s\n",barcode,get(barcodes_use,target));
  }

  /* close file */
  fclose(file);


/***************
* finalization *
***************/

  /* empty file reading buffer */
  free(barcode);

  /* clean-up by free'ing previously malloc'ed memory */
  for(i=0;i<n_barcodes_use;++i){
    free(barcodes_use[i]);
    delete_hit_list(top_hits[i]);
  }
  free(barcodes_use);
  free(top_hits);

  /*                                                                                       * RPAA *
   * DO NOT CALL hdestroy().                                                               * RPAA *
   *                                                                                       * RPAA *
   * With BSD libc, hdestroy() would call free() with each key in                          * RPAA *
   * table. Our keys are static strings, so free() would crash.                            * RPAA *
   */                                                                                     /* RPAA */

  /* end reporting success */
  return(EXIT_SUCCESS);
}
