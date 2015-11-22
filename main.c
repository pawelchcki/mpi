#include <mpi.h>
#define _XOPEN_SOURCE
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <assert.h>

/* 4 literowe
aa5UYq6trT5u.
bahAZ9Hk7SCf6
ddoo3WocSpthU
jkmD2RlhoMyuA
zzm4NUIIb7VIk
kkv864igyJC9o

5 literowe
aaSPfLTmjh3fU

6 literowe
aaLTdQr7DyHuU
*/

char *stro="aa5UYq6trT5u.";

#define PSIZE 4
#define ROOT 0

void root(){
//    char[PSIZE] passwd;
    char passwd[PSIZE+1];
    MPI_Status status;
    MPI_Recv(passwd, PSIZE+1, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG,MPI_COMM_WORLD, &status);
    printf("Haslo odebrane: %s od %d", passwd, status.MPI_SOURCE);
}

char min(char a, char b){
    return (a>b)?b:a;
}

bool check(char *strno, char *cmp, char *salt){
    char * x=crypt(cmp, salt);

    if ((strcmp(x,stro))==0){
        return true;
    }
    return false;
}

#define FIRST 'A'
#define LAST 'z'

bool crack(char *strno, char *cmp, char *salt, int position){
    if (strno[position] == '\0'){
        return false;
    }
    for(char c = FIRST; c<=LAST; c++ ){
        cmp[position] = c;
        if (check(strno, cmp, salt)){
            return true;
        } else {
            if (crack(strno, cmp, salt, position + 1)){
                return true;
            }
        }
    }
    return false;
}


int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int size,rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    char cmp[5]={0};

    char salt[3]={0};
    salt[0]=stro[0];
    salt[1]=stro[1];

    char all = LAST-FIRST;
    int num_workers = size - 1;
    int worker_rank = rank - 1;

    assert(num_workers > 0);
    printf("saffsa\n");

    char batch_size = ceill((double)all/num_workers);

    if (rank == ROOT){
        printf("afsfasfsa");
        root();
    } else {

        assert(worker_rank >= 0);
        char start = FIRST + worker_rank*batch_size;
        printf("afsfasfsa");

        char end = min(start + batch_size, LAST);
        printf("aaaa - %hhu, %hhu", FIRST, LAST);
fflush(stdout);
        for(char i=start; i<=end; i++){
            cmp[0] = i;
            if (crack(stro, cmp, salt, 1)){
                printf("znaleziono %s\n", cmp);
            }
        }

    }
    MPI_Finalize();
}
