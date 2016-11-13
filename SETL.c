#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <mpi.h>

/***********************************************************
  Helper functions 
***********************************************************/
void distributePattern(char **pattern, int patternSize, int numTasks, int iterations);

void distributeWorld(char** curWorld, int size, int numTasks, int curIteration);

void distributeSearch(char** curWorld, int size, int pSize, int numTasks, int iterations, int curIteration);

//For exiting on error condition
void die(int lineNo);

//For trackinng execution
long long wallClockTime();

int receiveWorld(char **curSmallWorld, int iter, int rank, int size, int numTasks);

int receiveSearch(char **curSmallWorld, int size, int numTasks, int rowsReceived, int iterations, int curIteration, int rank);

int worldRowsToReceive(int size, int numTasks, int rank);

int searchRowsToReceive(int size, int pSize, int numTasks, int rank);

void returnWorld(char** nextSmallWorld, int size, int numTasks, int rowsReceived, int curIteration, int rank);

void receiveReturnWorld(char **nextW, int size, int numTasks, int curIteration);

void receivePattern(char **pattern, int patternSize, int iterations, int numTasks, int rank);
/***********************************************************
  Square matrix related functions, used by both world and pattern
***********************************************************/

char** allocateRektMatrix( int size, int rows, char defaultValue );

char** allocateSquareMatrix( int size, char defaultValue );

void freeSquareMatrix( char** );

void printSquareMatrix( char**, int size );

void printRektMatrix(char** matrix, int size, int row);


/***********************************************************
   World  related functions
***********************************************************/

#define ALIVE 'X' 
#define DEAD 'O'

char** readWorldFromFile( char* fname, int* size );

int countNeighbours(char** world, int row, int col);

void evolveWorld(char** curWorld, char** nextWorld, int size, int rows);


/***********************************************************
   Simple circular linked list for match records
***********************************************************/

typedef struct MSTRUCT {
    int iteration, row, col, rotation;
    struct MSTRUCT *next;
} MATCH;


typedef struct {
    int nItem;
    MATCH* tail;
} MATCHLIST;

MATCHLIST* newList();

void deleteList( MATCHLIST*);

void insertEnd(MATCHLIST*, int, int, int, int);

void printList(MATCHLIST*);

MATCHLIST* joinLists(MATCHLIST* list0, MATCHLIST* list1);

/***********************************************************
   Search related functions
***********************************************************/

//Using the compass direction to indicate the rotation of pattern
#define N 0 //no rotation
#define E 1 //90 degree clockwise
#define S 2 //180 degree clockwise
#define W 3 //90 degree anti-clockwise

char** readPatternFromFile( char* fname, int* size );

void rotate90(char** current, char** rotated, int size);

void searchPatterns(char** world, int wSize, int numRows, int iteration, 
        char** patterns[4], int pSize, int* countBuffer, int* foundCount, int rank, int rowsPerTask);

void searchSinglePattern(char** world, int wSize, int numRows,  int interation,
        char** pattern, int pSize, int rotation, int* countBuffer, int* foundCount, int rank, int rowsPerTask);

/***********************************************************
   Main function
***********************************************************/


int main( int argc, char** argv)
{   
    int numtasks, rank, iter, iterations, task;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    int* paramsRecv = (int*)malloc(sizeof(int) * 3);
    char* patternsRecv;
    char **patterns[4];
    MPI_Status Stat;
    int size, patternSize, rowsReceived;
    long long before, after;    
    MATCHLIST*list;
    char **curW, **nextW, **temp, dummy[20];
    char **curSmallWorld, **nextSmallWorld;
    char **smallSearch;
    int dir, count;


    if(rank == 0){
        int params[3];
        
        if (argc < 4 ){
            fprintf(stderr, 
                "Usage: %s <world file> <Iterations> <pattern file>\n", argv[0]);
            exit(1);
        } 


        curW = readWorldFromFile(argv[1], &size);
        nextW = allocateSquareMatrix(size+2, DEAD);


        printf("World Size = %d\n", size);

        iterations = atoi(argv[2]);
        printf("Iterations = %d\n", iterations);

        patterns[N] = readPatternFromFile(argv[3], &patternSize);
        for (dir = E; dir <= W; dir++){
            patterns[dir] = allocateSquareMatrix(patternSize, DEAD);
            rotate90(patterns[dir-1], patterns[dir], patternSize);
        }
        printf("Pattern size = %d\n", patternSize);

#ifdef DEBUG
        printSquareMatrix(patterns[N], patternSize);
        printSquareMatrix(patterns[E], patternSize);
        printSquareMatrix(patterns[S], patternSize);
        printSquareMatrix(patterns[W], patternSize);
#endif

        params[0] = size;
        params[1] = patternSize;
        params[2] = iterations;
        //param[3] = searchCount;
        //param[4] = evolveCount;
        for(task = 1; task<numtasks; task++){
            MPI_Send(params, 3, MPI_INT, task, task, MPI_COMM_WORLD);
        }
    }

    if(rank!=0){
        MPI_Recv(&paramsRecv[0], 3, MPI_INT, 0, rank, MPI_COMM_WORLD, &Stat);       
    

    size = paramsRecv[0];
    patternSize = paramsRecv[1];
    iterations = paramsRecv[2];

    int *answers = (int*)malloc(sizeof(int) * searchRowsToReceive(size, patternSize, numtasks, rank));

//    printf("Rank: %d - Received size=%d, pSize=%d, iterations=%d\n", rank, size, patternSize, iterations);
    }
    if(rank == 0){
        //printf("Rank: %d - Entering distribute pattern\n", rank);
        distributePattern(patterns[N], patternSize, numtasks, iterations);
        //printf("Rank: %d - Outside distribute pattern, pattern distributed\n", rank);
     
        //Start timer
        before = wallClockTime();

        //Actual work start

    }
    list = newList();

    if(rank !=0){
    patterns[N] = allocateSquareMatrix(patternSize, DEAD);
    
    receivePattern(patterns[N], patternSize, iterations, numtasks, rank);

    for (dir = E; dir <= W; dir++){
        patterns[dir] = allocateSquareMatrix(patternSize, DEAD);
        rotate90(patterns[dir-1], patterns[dir], patternSize);
    }

  //  printf("Rank: %d - Pattern received outside function:\n", rank);

#ifdef DEBUG
 //   printSquareMatrix(patterns[N], patternSize);
#endif
}
    int* countBuffer;

    for (iter = 0; iter < iterations; iter++){
        int foundCount = 0;

        if(rank == 0){
    #ifdef DEBUG
            printf("World Iteration.%d\n", iter);
            printSquareMatrix(curW, size+2);
    #endif       
            distributeSearch(curW, size, patternSize, numtasks, iterations, iter);
        }
        
        if(rank !=0){
        int searchToReceive = searchRowsToReceive(size, patternSize, numtasks, rank);
        smallSearch = allocateRektMatrix(size+2, searchToReceive, DEAD);
        rowsReceived = receiveSearch(smallSearch, size, numtasks, searchToReceive, iterations, iter, rank);
//        printf("Rank: %d - rowsReceived: %d\n", rank, rowsReceived);
        countBuffer = (int*) malloc(sizeof(int*) * rowsReceived * size);

        searchPatterns(smallSearch, size, rowsReceived, iter, patterns, patternSize, countBuffer, &foundCount, rank, (size - patternSize + 1)/(numtasks-1));

        MPI_Send(&countBuffer[0], foundCount, MPI_INT, 0, (iterations+numtasks*2) + (iterations*numtasks) + numtasks + (iter*numtasks)+rank, MPI_COMM_WORLD);
        
    //    printf("RANK: %d ITER: %d SEARCH FINISH LIAO FOUNDCOUNT = %d\n", rank, iter, foundCount );
}
        if(rank == 0){
            int t, number_amount, k;
            MATCHLIST * iterList = newList();
            MATCHLIST* temp0 = newList();
            MATCHLIST* temp1 = newList();
            MATCHLIST* temp2 = newList();
            MATCHLIST* temp3 = newList();
            for(t = 1; t<numtasks; t++){
                number_amount = 0;
                MPI_Probe(t, (iterations+numtasks*2) + (iterations*numtasks) + numtasks + (iter*numtasks) + t, MPI_COMM_WORLD, &Stat);
                MPI_Get_count(&Stat, MPI_INT, &number_amount);
                int* number_buf = (int*)malloc(sizeof(int) * number_amount);
                MPI_Recv(number_buf, number_amount, MPI_INT, t, (iterations+numtasks*2) + (iterations*numtasks) + numtasks + (iter*numtasks) + t,
             MPI_COMM_WORLD, &Stat);
                for(k = 0; k<number_amount; k+=4){
                    if(number_buf[k+3]==0){
                        insertEnd(temp0, number_buf[k], number_buf[k+1], number_buf[k+2], number_buf[k+3]);
                    }else if(number_buf[k+3]==1){
                        insertEnd(temp1, number_buf[k], number_buf[k+1], number_buf[k+2], number_buf[k+3]);
                    }else if(number_buf[k+3]==2){
                        insertEnd(temp2, number_buf[k], number_buf[k+1], number_buf[k+2], number_buf[k+3]);
                    }else if(number_buf[k+3]==3){
                        insertEnd(temp3, number_buf[k], number_buf[k+1], number_buf[k+2], number_buf[k+3]);
                    }
                }
            }
            iterList = joinLists(temp0, temp1);
            iterList = joinLists(iterList, temp2);
            iterList = joinLists(iterList, temp3);
            list = joinLists(list, iterList);
        //    printf("Rank: 0 - Before enter distribute world\n");
            distributeWorld(curW, size, numtasks, iter);
        }
   if(rank!=0){     
        int rowsToReceive = worldRowsToReceive(size, numtasks, rank);
        curSmallWorld = allocateRektMatrix(size+2, rowsToReceive, DEAD);
      //  printf("Rank: %d - Allocated %d rows to memory and entering receiveWorld\n", rank, rowsToReceive);
        rowsReceived = receiveWorld(curSmallWorld, size, numtasks, iter, rank);
      //  printf("Rank: %d - Outside receiveWorld function, task received %d rows from root in iteration %d:\n", rank, rowsReceived,iter);
   //     printRektMatrix(curSmallWorld, size, rowsReceived);

        nextSmallWorld = allocateRektMatrix(size+2, rowsReceived, DEAD);
       // printf("Rank: %d - Finished printing received matrix\n", rank);

        //Generate next generation
  //      printf("Rank: %d - WHAT?!?! Small world is evolving!\n", rank);
        evolveWorld( curSmallWorld, nextSmallWorld, size, rowsToReceive-2 );
       // printf("Rank: %d - World has evolved!\n", rank);

        returnWorld(nextSmallWorld, size, numtasks, rowsReceived, iter, rank);
    }   
        if(rank == 0){
            receiveReturnWorld(nextW, size, numtasks, iter);

            temp = curW;
            curW = nextW;
            nextW = temp;
        }
    }

    if(rank == 0){
        printList( list );

        //Stop timer
        after = wallClockTime();

        printf("Sequential SETL took %1.2f seconds\n", 
            ((float)(after - before))/1000000000);


        //Clean up
        deleteList( list );

        freeSquareMatrix( curW );
        freeSquareMatrix( nextW );

        freeSquareMatrix( patterns[0] );
        freeSquareMatrix( patterns[1] );
        freeSquareMatrix( patterns[2] );
        freeSquareMatrix( patterns[3] );
    
    }
    MPI_Finalize();
    return 0;
}

/***********************************************************
  Helper functions 
***********************************************************/
void receiveReturnWorld(char **nextW, int size, int numTasks, int curIteration){
    int i, rowsPerTask, count;
    MPI_Status Stat;

    rowsPerTask = size/(numTasks-1);

    count = rowsPerTask * (size + 2);

    for (i = 1; i<numTasks; i++){
        if(i != (numTasks - 1)){
            MPI_Recv(&nextW[(i-1)*rowsPerTask+1][0], count, MPI_CHAR, i, (curIteration * numTasks * 2) + (2*numTasks) +i, MPI_COMM_WORLD, &Stat);
        }else{
            MPI_Recv(&nextW[(i-1)*rowsPerTask+1][0], (size - (rowsPerTask * (i-1))) * (size + 2), MPI_CHAR, i, (curIteration * numTasks * 2) + (2*numTasks) +i, MPI_COMM_WORLD, &Stat);
        }
    }
}

void returnWorld(char** nextSmallWorld, int size, int numTasks, int rowsReceived, int curIteration, int rank){
    int tag, count;
    tag = (curIteration * numTasks * 2) + (2 * numTasks) + rank;
    count = (rowsReceived - 2) * (size +2);

    MPI_Send(&nextSmallWorld[1][0], count, MPI_CHAR, 0, tag, MPI_COMM_WORLD);
}

int worldRowsToReceive(int size, int numTasks, int rank){
    int rowsPerTask, rowsReceived;

    rowsPerTask = size/(numTasks-1);

    if(rank == (numTasks -1)){
        rowsReceived = size - (rowsPerTask * (rank-1)) + 2;    
    }else{
        rowsReceived = rowsPerTask + 2;
    }

    return rowsReceived;    
}

int searchRowsToReceive(int size, int pSize, int numTasks, int rank){
    int rowsPerTask, rowsToReceive;

    rowsPerTask = (size - (pSize - 1)) / (numTasks-1);
    if(rank == (numTasks -1)){
        rowsToReceive = size - ((numTasks -1) * rowsPerTask) + rowsPerTask + 2;
    }else{
        rowsToReceive = rowsPerTask + pSize;
    }
    

    return rowsToReceive;    
}

int receiveWorld(char **curSmallWorld, int size, int numTasks, int curIteration, int rank){
    char* recvWorldBuffer;
    MPI_Status Stat;
    int rowsPerTask, elementCount, rowsReceived, i;

    rowsPerTask = size/(numTasks-1);

    if(rank == (numTasks -1)){
        rowsReceived = size - (rowsPerTask * (rank-1)) + 2;    
    }else{
        rowsReceived = rowsPerTask + 2;
    }
    
    elementCount = rowsReceived*(size +2) ;    
    recvWorldBuffer = (char*) malloc(sizeof(char*) * elementCount);
   // printf("Rank: %d - Waiting to receive %d rows from root....\n", rank, rowsReceived);
    MPI_Recv(&recvWorldBuffer[0], elementCount, MPI_CHAR, 0, (curIteration*numTasks*2)+ numTasks + rank, MPI_COMM_WORLD, &Stat);
   // printf("Rank: %d - Received %d rows from root - SUCCESS\n", rank, rowsReceived);

    curSmallWorld[0] = recvWorldBuffer;
    for (i = 1; i < rowsReceived; i++){
        curSmallWorld[i] = &recvWorldBuffer[i*(size+2)];
    }    

    return rowsReceived;
}

int receiveSearch(char **curSmallWorld, int size, int numTasks, int rowsReceived, int iterations, int curIteration, int rank){
    char* recvWorldBuffer;
    MPI_Status Stat;
    int elementCount, i;
    
    elementCount = rowsReceived*(size +2) ;    
    recvWorldBuffer = (char*) malloc(sizeof(char*) * elementCount);
   // printf("Rank: %d - Waiting to receive %d SEARCH rows from root....\n", rank, rowsReceived);
    MPI_Recv(&recvWorldBuffer[0], elementCount, MPI_CHAR, 0, (iterations*numTasks*2)+ (curIteration * numTasks) +rank, MPI_COMM_WORLD, &Stat);
   // printf("Rank: %d - Received %d SEARCH rows from root - SUCCESS\n", rank, rowsReceived);

    curSmallWorld[0] = recvWorldBuffer;
    for (i = 1; i < rowsReceived; i++){
        curSmallWorld[i] = &recvWorldBuffer[i*(size+2)];
    }    

    return rowsReceived;
}

void receivePattern(char **pattern, int patternSize, int iterations, int numTasks, int rank){
    char* patternRecvBuffer = (char*)malloc(sizeof(char) * patternSize*patternSize);
    MPI_Status Stat;
    int i;
   // printf("Rank: %d - Waiting to receive pattern...\n", rank);

    MPI_Recv(&patternRecvBuffer[0], patternSize*patternSize, MPI_CHAR, 0, (iterations*numTasks*2) + (iterations*numTasks) + (2*numTasks) + rank, MPI_COMM_WORLD, &Stat);
   // printf("Rank: %d - Pattern received - SUCCESS\n", rank);

    pattern[0] = patternRecvBuffer;
    for(i = 1; i<patternSize; i++){
        pattern[i] = &patternRecvBuffer[i*patternSize];
    }

   // printf("Rank: %d - Pattern received in function:\n", rank);
  //  printSquareMatrix(pattern, patternSize);
}

void distributePattern(char **pattern, int patternSize, int numTasks, int iterations){
    int task; 

   // printf("Rank: 0 - Inside distribute pattern\n");

    for(task = 1; task<numTasks; task++){
   //     printf("Rank: 0 - distributing pattern to task %d....\n", task);
        MPI_Send(&pattern[0][0], patternSize*patternSize, MPI_CHAR, 
            task, (iterations*numTasks*2) + (iterations*numTasks) + (2*numTasks) + task, MPI_COMM_WORLD);
   //     printf("Rank: 0 - distribute pattern to task %d - SUCCESS\n", task);        
    }
}

void distributeWorld(char** curWorld, int size, int numTasks, int curIteration){
    int rowsPerTask, rowsGiven, elementCount, i;

   // printf("Rank 0 - Inside distribute world\n");
    rowsPerTask = size/(numTasks-1);
    rowsGiven = rowsPerTask + 2;
    elementCount = rowsGiven * (size + 2);

    for (i = 1; i < numTasks; i++){
  //      printf("Rank 0 - Distributing world to task: %d....\n", i);
        if(i != (numTasks -1)){
            MPI_Send(&curWorld[rowsPerTask*(i-1)][0], elementCount, MPI_CHAR, i, (curIteration*numTasks*2)+ numTasks +i, MPI_COMM_WORLD);
        }else{
            MPI_Send(&curWorld[rowsPerTask*(i-1)][0], (size - (rowsPerTask * (i-1)) + 2) * (size + 2), MPI_CHAR, i, (curIteration*numTasks*2) + numTasks +i, MPI_COMM_WORLD);
        }
    //    printf("Rank 0 - Distributed world to task: %d - SUCCESS\n", i);
    }    
}

void distributeSearch(char** curWorld, int size, int pSize, int numTasks, int iterations, int curIteration){
    int rowsPerTask, rowsGiven, elementCount, i;

   // printf("Rank 0 - Inside distribute search\n");
    rowsPerTask = (size - (pSize - 1)) / (numTasks-1);
    rowsGiven = rowsPerTask + pSize;
    elementCount = rowsGiven * (size + 2);

    for (i = 1; i < numTasks; i++){
    //    printf("Rank 0 - Distributing search to task: %d....\n", i);
        if(i != (numTasks -1)){
            MPI_Send(&curWorld[rowsPerTask*(i-1)][0], elementCount, MPI_CHAR, i, (iterations*numTasks*2)+ (curIteration * numTasks) +i, MPI_COMM_WORLD);
        }else{
            rowsGiven = size - ((numTasks-1) * rowsPerTask) + rowsPerTask + 2;
            elementCount = rowsGiven * (size + 2);        
            MPI_Send(&curWorld[rowsPerTask*(i-1)][0], elementCount, MPI_CHAR, i, (iterations*numTasks*2)+ (curIteration * numTasks) +i, MPI_COMM_WORLD);
        }
      //  printf("Rank 0 - Distributed search to task: %d - SUCCESS\n", i);
    }    
}

void die(int lineNo)
{
    fprintf(stderr, "Error at line %d. Exiting\n", lineNo);
    exit(1);
}

long long wallClockTime( )
{
#ifdef __linux__
    struct timespec tp;
    clock_gettime(CLOCK_REALTIME, &tp);
    return (long long)(tp.tv_nsec + (long long)tp.tv_sec * 1000000000ll);
#else
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (long long)(tv.tv_usec * 1000 + (long long)tv.tv_sec * 1000000000ll);
#endif
}

/***********************************************************
  Square matrix related functions, used by both world and pattern
***********************************************************/


char** allocateRektMatrix( int size, int rows, char defaultValue )
{
    //printf("Rank: %d - Just entered rekt\n", rank);
    char* contiguous;
    char** matrix;
    int i;

    //Using a least compiler version dependent approach here
    //C99, C11 have a nicer syntax.    
    contiguous = (char*) malloc(sizeof(char) * size * rows);
    if (contiguous == NULL) 
        die(__LINE__);


    memset(contiguous, defaultValue, size * rows );

    //Point the row array to the right place
    matrix = (char**) malloc(sizeof(char*) * rows );
    if (matrix == NULL) 
        die(__LINE__);

    matrix[0] = contiguous;
    for (i = 1; i < rows; i++){
        matrix[i] = &contiguous[i*size];
    }

    //printf("Rank: %d - Leaving rekt\n", rank);
    return matrix;
}

char** allocateSquareMatrix( int size, char defaultValue )
{

    char* contiguous;
    char** matrix;
    int i;

    //Using a least compiler version dependent approach here
    //C99, C11 have a nicer syntax.    
    contiguous = (char*) malloc(sizeof(char) * size * size);
    if (contiguous == NULL) 
        die(__LINE__);


    memset(contiguous, defaultValue, size * size );

    //Point the row array to the right place
    matrix = (char**) malloc(sizeof(char*) * size );
    if (matrix == NULL) 
        die(__LINE__);

    matrix[0] = contiguous;
    for (i = 1; i < size; i++){
        matrix[i] = &contiguous[i*size];
    }

    return matrix;
}

void printRektMatrix(char** matrix, int size, int row){
    int i, j;
    for (i = 0; i < row; i++){
        for (j = 0; j < size; j++){
            printf("%c", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void printSquareMatrix( char** matrix, int size )
{
    int i,j;
    
    for (i = 0; i < size; i++){
        for (j = 0; j < size; j++){
            printf("%c", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void freeSquareMatrix( char** matrix )
{
    if (matrix == NULL) return;

    free( matrix[0] );
}

/***********************************************************
   World  related functions
***********************************************************/

char** readWorldFromFile( char* fname, int* sizePtr )
{
    FILE* inf;
    
    char temp, **world;
    int i, j;
    int size;

    inf = fopen(fname,"r");
    if (inf == NULL)
        die(__LINE__);


    fscanf(inf, "%d", &size);
    fscanf(inf, "%c", &temp);
    
    //Using the "halo" approach
    // allocated additional top + bottom rows
    // and leftmost and rightmost rows to form a boundary
    // to simplify computation of cell along edges
    world = allocateSquareMatrix( size + 2, DEAD );

    for (i = 1; i <= size; i++){
        for (j = 1; j <= size; j++){
            fscanf(inf, "%c", &world[i][j]);
        }
        fscanf(inf, "%c", &temp);
    }

    *sizePtr = size;    //return size
    return world;
    
}

int countNeighbours(char** world, int row, int col)
//Assume 1 <= row, col <= size, no check 
{
    int i, j, count;

    count = 0;
    for(i = row-1; i <= row+1; i++){
        for(j = col-1; j <= col+1; j++){
            count += (world[i][j] == ALIVE );
        }
    }

    //discount the center
    count -= (world[row][col] == ALIVE);

    return count;

}

void evolveWorld(char** curWorld, char** nextWorld, int size, int rows)
{
    int i, j, liveNeighbours;

    for (i = 1; i <= rows; i++){
        for (j = 1; j <= size; j++){
            liveNeighbours = countNeighbours(curWorld, i, j);
            nextWorld[i][j] = DEAD;

            //Only take care of alive cases
            if (curWorld[i][j] == ALIVE) {

                if (liveNeighbours == 2 || liveNeighbours == 3)
                    nextWorld[i][j] = ALIVE;

            } else if (liveNeighbours == 3)
                    nextWorld[i][j] = ALIVE;
        } 
    }
}

void evolveWorldSlave(){

}
/***********************************************************
   Search related functions
***********************************************************/

char** readPatternFromFile( char* fname, int* sizePtr )
{
    FILE* inf;
    
    char temp, **pattern;
    int i, j;
    int size;

    inf = fopen(fname,"r");
    if (inf == NULL)
        die(__LINE__);


    fscanf(inf, "%d", &size);
    fscanf(inf, "%c", &temp);
    
    pattern = allocateSquareMatrix( size, DEAD );

    for (i = 0; i < size; i++){
        for (j = 0; j < size; j++){
            fscanf(inf, "%c", &pattern[i][j]);
        }
        fscanf(inf, "%c", &temp);
    }
    
    *sizePtr = size;    //return size
    return pattern;
}


void rotate90(char** current, char** rotated, int size)
{
    int i, j;

    for (i = 0; i < size; i++){
        for (j = 0; j < size; j++){
            rotated[j][size-i-1] = current[i][j];
        }
    }
}

void searchPatterns(char** world, int wSize, int rowsReceived, int iteration, 
        char** patterns[4], int pSize, int* countBuffer, int* foundCount, int rank, int rowsPerTask)
{
    int dir;

    for (dir = N; dir <= W; dir++){
        searchSinglePattern(world, wSize, rowsReceived, iteration, 
                patterns[dir], pSize, dir, countBuffer, foundCount, rank, rowsPerTask);
    }

}

void searchSinglePattern(char** world, int wSize, int numRows, int iteration,
        char** pattern, int pSize, int rotation, int* countBuffer, int* foundCount, int rank, int rowsPerTask)
{
    int wRow, wCol, pRow, pCol, match;

    int offset = (rank-1)*rowsPerTask;

   // printf("INSIDE SINGLE PATTERN FOUNDCOUNT: %d\n", *foundCount);
    for (wRow = 1; wRow <= ((numRows-1)-pSize+1); wRow++){
        for (wCol = 1; wCol <= (wSize-pSize+1); wCol++){
            match = 1;
#ifdef DEBUGMORE
            printf("S:(%d, %d)\n", wRow-1, wCol-1);
#endif
            for (pRow = 0; match && pRow < pSize; pRow++){
                for (pCol = 0; match && pCol < pSize; pCol++){
                    if(world[wRow+pRow][wCol+pCol] != pattern[pRow][pCol]){
#ifdef DEBUGMORE
                        printf("\tF:(%d, %d) %c != %c\n", pRow, pCol,
                            world[wRow+pRow][wCol+pCol], pattern[pRow][pCol]);
#endif
                        match = 0;    
                    }
                }
            }
            if (match){
                //insertEnd(list, iteration, wRow-1, wCol-1, rotation);
                countBuffer[*foundCount]=iteration;
                *foundCount +=1;
                countBuffer[*foundCount]=wRow-1 + offset;
                *foundCount +=1;
                countBuffer[*foundCount]=wCol-1;
                *foundCount +=1;
                countBuffer[*foundCount]=rotation;
                *foundCount +=1;
#ifdef DEBUGMORE
printf("*** Row = %d, Col = %d\n", wRow-1, wCol-1);
#endif
            }
        }
    }
}

/***********************************************************
   Simple circular linked list for match records
***********************************************************/

MATCHLIST* newList()
{
    MATCHLIST* list;

    list = (MATCHLIST*) malloc(sizeof(MATCHLIST));
    if (list == NULL)
        die(__LINE__);

    list->nItem = 0;
    list->tail = NULL;

    return list;
}

void deleteList( MATCHLIST* list)
{
    MATCH *cur, *next;
    int i;
    //delete items first

    if (list->nItem != 0 ){
        cur = list->tail->next;
        next = cur->next;
        for( i = 0; i < list->nItem; i++, cur = next, next = next->next ) {
            free(cur); 
        }

    }
    free( list );
}

void insertEnd(MATCHLIST* list, 
        int iteration, int row, int col, int rotation)
{
    MATCH* newItem;

    newItem = (MATCH*) malloc(sizeof(MATCH));
    if (newItem == NULL)
        die(__LINE__);

    newItem->iteration = iteration;
    newItem->row = row;
    newItem->col = col;
    newItem->rotation = rotation;

    if (list->nItem == 0){
        newItem->next = newItem;
        list->tail = newItem;
    } else {
        newItem->next = list->tail->next;
        list->tail->next = newItem;
        list->tail = newItem;
    }

    (list->nItem)++;

}

MATCHLIST* joinLists(MATCHLIST* list0, MATCHLIST* list1){
    MATCHLIST* joined;
    if((list0->nItem ==0)&&(list1->nItem ==0)){
        return joined;
    }else if(list0->nItem ==0){
        joined = list1;
        return joined;
    }else if(list1->nItem ==0){
        joined = list0;
        return joined;
    }else{
        MATCH* tempItem;

        tempItem = (MATCH*) malloc(sizeof(MATCH));
        if (tempItem == NULL)
            die(__LINE__);      

        tempItem->next = list0->tail->next;
        list0->tail->next = list1->tail->next;
        list1->tail->next = tempItem->next;
        joined = list1;
        joined->nItem = list0->nItem + list1->nItem;
        return joined; 
    }
}

void printList(MATCHLIST* list)
{
    int i;
    MATCH* cur;

    printf("List size = %d\n", list->nItem);    


    if (list->nItem == 0) return;

    cur = list->tail->next;
    for( i = 0; i < list->nItem; i++, cur=cur->next){
        printf("%d:%d:%d:%d\n", 
                cur->iteration, cur->row, cur->col, cur->rotation);
    }
}

