// This file requires at least C99 to compile

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h> // memset()


// unsigned int type used
typedef unsigned long long bi_block;

// Structure of bigint
#define BIGINT_LENGTH 5
typedef struct {
    size_t length;
    //bi_block blocks[BIGINT_LENGTH];
    bi_block* blocks;
} BigInt;

// Number of LSB used: (4 bits/byte).
#define bi_block_used (4 * sizeof(bi_block))  

// Base used (must be power of 2).
bi_block const bi_base = ((bi_block) 1) << bi_block_used;

// Results of comparisonss
enum bi_cmp {
    bi_below = -1,
    bi_equal = 0,
    bi_above = 1
};

/** Initialize a number with the given blocks
 * @param bi Number to initialize
 * @param val Blocs to assign, NULL for 0
 * @param len Number of elements in the table val
 * @return True if success
**/
bool bi_init(BigInt* bi, bi_block const* val, size_t len);

/** Destructor
 * @param Number to destroy, must have been initialized
**/
void bi_cleanup(BigInt*);

/** Comparison of two BigInts
 * @param Number A, must have been initialized
 * @param Number B, must have been initialized
 * @return 'bi_below' if A < B, 'bi_above' if A > B, 'bi_equal' if A == B
**/
enum bi_cmp bi_compare(BigInt const* A, BigInt const* B);

/** Modulo of BigInt by one bloc.
 * @param Number A, must have been initialized
 * @param Modulo m, positive
 * @return A % m
**/
bi_block bi_mod(BigInt const* A, bi_block m);

/** Modular exponentiation.
 * @param Base B, must have been initialized
 * @param Exponant E, must have been initialized
 * @param Modulo m, positive (return 0 if m is null)
 * @return B^E % m
**/
bi_block bi_expmod(BigInt const* B, BigInt const* E, bi_block m);

/** Addition of two BigInts, initialize a third as result.
 * @param Number A, must have been initialized
 * @param Number B, must have been initialized
 * @param Number R, initialized with value A + B
 * @return True if success
**/
bool bi_sum_init(BigInt const* A, BigInt const* B, BigInt* R);

/** Addition of two BigInts, place result in A
 * @param Nombre A, will be replaced by A + B, must have been initialized
 * @param Number B, must have been initialized
 * @return True if success
**/
bool bi_sum_over(BigInt* A, BigInt const* B);

/** Multiplication of two BigInts, initialize a third as result.
 * @param Number A, must have been initialized
 * @param Number B, must have been initialized
 * @param Number R, initialized with value A*B
 * @return True if success
**/
bool bi_mul_init(BigInt const* A, BigInt const* B, BigInt* R);


/** Multiplication of two BigInts, place result in A
 * @param Nombre A, will be replaced by A*B, must have been initialized
 * @param Number B, must have been initialized
 * @return True if success
**/
bool bi_mul_over(BigInt* A, BigInt const* B);

// -------------------------------------------------------------------------- //
// Fundamental operations (operations over 'bi_block')

/** Compute x modulo base
 * @param x Value of bi_block
 * @return x % base
**/
bi_block bi_lower(bi_block x) {
    return x % bi_base;
}

/** Compute x/base
 * @param x Value of bi_block
 * @return x % base
**/
bi_block bi_upper(bi_block x) {
    return x / bi_base;
}

// -------------------------------------------------------------------------- //

/** Print a BigInt.
 * @param p Prefix string
 * @param x Number given, must have been initialized
**/
void bi_println(const char* p, BigInt const* x) {
    printf("%s{", p);
    for (size_t i = 0; i < x->length; ++i)
        printf(" 0x%08llx", x->blocks[i]);
    puts(" }");
}

/* ===========================================================================
 * Definition of functions
 * =========================================================================== */
#define ERR_INVALID_ARGUMENT false;
#define ERR_NONE true;


bool bi_init(BigInt* bi, bi_block const* val, size_t len){

    // if non valid arguments
    if (val == NULL || bi == NULL){
        bi_cleanup(bi);
        return ERR_INVALID_ARGUMENT;
    }

    // if valid arguments

    //case len = 0
    if(len == 0){
        bi->length = 0;
        bi->blocks = NULL;
    }
    // otherwise
    else {
        bi->length = len;
        //allocate memory
        bi->blocks = calloc(len, sizeof(bi_block));
        //fill in values
        for(size_t i = 0; i < len; i++){
            bi->blocks[i] = val[i];
        }
    }
    
    return ERR_NONE;
}

void bi_cleanup(BigInt* bi){
    if(bi != NULL){
        memset(bi->blocks, 0, sizeof(bi_block)*bi->length);
        bi->length = 0;    
    }
}

enum bi_cmp bi_compare(BigInt const* A, BigInt const* B){

    // check if A and B have been initialized
    if(A == NULL || B == NULL){
        return bi_equal;
    }
    else if(A->blocks == NULL || B->blocks == NULL){
        return bi_equal;
    }


    // if A and B have different length, check if the bigger doesnt have 0 on higher power
    if(A->length > B->length){
        for(size_t i = B->length; i < A->length; i++) {
            if(A->blocks[i] > 0) return bi_above;
        }
    }
    if(A->length < B->length) {
        for(size_t i = A->length; i < B->length; i++) {
            if(B->blocks[i] > 0) return bi_below;
        }
    }
    

    // A and B same size or different sizes and the bigger has 0 on higher power
    // check each bi_block starting from highest power
    int counter = (A->length < B->length) ? (A->length - 1) : (B->length - 1); //min of lengths

    while(counter >= 0){
        if(A->blocks[counter] > B->blocks[counter]) return bi_above;
        else if(A->blocks[counter] < B->blocks[counter]) return bi_below;
        counter -= 1;
    }

    return bi_equal;    
}


bi_block bi_mod(BigInt const* A, bi_block m){
    
    //check if A has been initialized
    if(A == NULL){
        return 0;
    }
    else if(A->blocks == NULL){
        return 0;
    }

    bi_block p = 1;
    bi_block s = 0;

    for(size_t i = 0; i < A->length - 1; i++){
        s = (s + (A->blocks[i] * p) % m) % m;
        p = (bi_base * p) % m; 
    }
    return s;
}


bi_block bi_expmod(BigInt const* B, BigInt const* E, bi_block m){
    // check if B and E have been initialized
    if(B == NULL || E == NULL){
        return 0;
    }
    else if(B->blocks == NULL || E->blocks == NULL){
        return 0;
    }

    bi_block p = 1;
    bi_block b = bi_mod(B, m);
    for(size_t i = 0; i < E->length; i++){

        bi_block e = E->blocks[i];
        
        for(size_t j = 0; j < bi_block_used; j++){
            if(e % 2 != 0){
                p = (p * b) % m;
            }
            e /= 2;
            b = (b*b) % m;
        }  
    }
    return p;
    
}

bool bi_sum_init(BigInt const* A, BigInt const* B, BigInt* R){
    
    if(bi_init(R, A->blocks, A->length)){
        return bi_sum_over(R, B);
    }
    return ERR_INVALID_ARGUMENT;
}

bool bi_sum_over(BigInt* A, BigInt const* B){
    // Check if A and B have been initialized
    if(A == NULL || B == NULL){
        return ERR_INVALID_ARGUMENT;
    }
    else if(A->blocks == NULL || B->blocks == NULL){
        return ERR_NONE;
    }

    //Compute minimum length of A and B and maximum length that A can have
    size_t minLength = (A->length > B->length) ? B->length : A->length;
    size_t maxLength = ((A->length < B->length) ? B->length : A->length) + 1;
    
    A->blocks = realloc(A->blocks, sizeof(bi_block)*maxLength);

    //Add the blocks of A and B until we reach one of their end
    bi_block carry = 0;
    bi_block oldCarry = 0;
    for(size_t i = 0; i < minLength; i++){
        carry = bi_upper(A->blocks[i] + B->blocks[i] + oldCarry);
        A->blocks[i] = bi_lower(A->blocks[i] + B->blocks[i] + oldCarry);
        oldCarry = carry;
    }
    
    

    //Add remaining blocks of the biggest bigint (keep carry value and transfer it)
    if(A->length > B->length){
         for(size_t i = minLength; i < maxLength-1; i++){
            carry = bi_upper(A->blocks[i] + oldCarry);
            A->blocks[i] = bi_lower(A->blocks[i] + oldCarry);
            oldCarry = carry;
        }
    }
    if(A->length < B->length){
        for(size_t i = minLength; i < maxLength-1; i++){
            carry = bi_upper(B->blocks[i] + oldCarry);
            A->blocks[i] = bi_lower(B->blocks[i] + oldCarry);
            oldCarry = carry;
        }
    } 

    A->blocks[maxLength-1] = carry;
    
    //if there was no carry in the end
    if(A->blocks[maxLength-1] == 0){
        maxLength -= 1;
        A->blocks = realloc(A->blocks, sizeof(bi_block)*maxLength);
    }

    //setting (new) length
    A->length = maxLength;


    return ERR_NONE; 
}


bool bi_mul_over(BigInt* A, BigInt const* B){

    BigInt R;
    R.blocks = calloc(0, sizeof(bi_block));
    R.length = 0;

    bool v1 = bi_mul_init(A, B, &R);
    bool v2 = bi_init(A, R.blocks, R.length);

    return v1 && v2;
}
    
    

bool bi_mul_init(BigInt const* A, BigInt const* B, BigInt* R){
    
    // Check if A and B have been initialized
    if(A == NULL || B == NULL){
        return ERR_INVALID_ARGUMENT;
    }
    else if(A->blocks == NULL || B->blocks == NULL){
        return ERR_NONE;
    }

    
    //Number of elements in C
    size_t sizeOfC = (A->length > B->length) ? B->length : A->length;
    
    //Size of biggest bigint
    size_t bigSize = (A->length > B->length) ? A->length+1 : B->length+1;
    
    BigInt C[sizeOfC];

    // Init C
    for(size_t i = 0; i < sizeOfC; i++){  
        bi_block* vals = calloc(bigSize + i, sizeof(bi_block));
        bi_init(&C[i], vals, bigSize + i);
    }

        
    bi_block carry = 0;
    bi_block oldCarry = 0;

    for(size_t i = 0; i < sizeOfC; i++){

        //Start loop from i (keep a 0 for the power)

        if(A->length > B->length){
            for(size_t j = 0; j < bigSize; j++){
                carry = bi_upper(A->blocks[j]*B->blocks[i] + oldCarry);
                C[i].blocks[j+i] = bi_lower(A->blocks[j]*B->blocks[i] + oldCarry);
                oldCarry = carry;
            }
        }
        else {
            for(size_t j = 0; j < bigSize; j++){
                carry = bi_upper(A->blocks[i]*B->blocks[j] + oldCarry);
                C[i].blocks[j+i] = bi_lower(A->blocks[i]*B->blocks[j] + oldCarry);
                oldCarry = carry;
            }
        }

        if(carry != 0){
            C[i].blocks[C[i].length] = carry;
        }
        else{
            C[i].length = C[i].length - 1;
            C[i].blocks = realloc(C[i].blocks, sizeof(bi_block)*C[i].length);
        }
        
        carry = 0;
        oldCarry = 0;
    }
    
    //Sum up the blocks from multiplication
    
    // Sum over the last (biggest) BigInt
    for(size_t i = 0; i < sizeOfC-1; i++){
        bi_sum_over(&C[sizeOfC-1], &C[i]);
    }

    //delete extra 0
    int counter = C[sizeOfC-1].length - 1;
    while(C[sizeOfC-1].blocks[counter] == 0 && counter != 0){
        C[sizeOfC-1].length = C[sizeOfC-1].length - 1;
        C[sizeOfC-1].blocks = realloc(C[sizeOfC-1].blocks, sizeof(bi_block)*C[sizeOfC-1].length);
        counter -= 1;
    }
    
    
    bi_init(R, C[sizeOfC-1].blocks, C[sizeOfC-1].length);
    
    return ERR_NONE;
}



// -------------------------------------------------------------------------- //
// Tests

/** Tests sommaires.
 * @return 0 if success, 1 otherwise
**/
int main(void) {
    { // Tests comparison and modular exponentiation
        bi_block m = 373;
        BigInt x, y, z;
        memset(&x, 0, sizeof(x));
        memset(&y, 0, sizeof(y));
        memset(&z, 0, sizeof(z));
        bi_block const vx[] = {0xffffffffull, 0xffffffffull, 0xffffffffull, 0xffffffffull};
        bi_block const vy[] = {0x0000ffffull, 0x0000ffffull, 0x0000ffffull, 0x0000ffffull, 0x0000ffffull, 0x0000ffffull};
        bi_block const vz[] = {0xffffffffull, 0xffffffffull, 0xffffffffull, 0xffffffffull, 0x00000000ull};
        if (!bi_init(&x, vx, sizeof(vx) / sizeof(*vx))) {
            puts("Error: bi_init(&x, ...)");
            return 1;
        }
        if (!bi_init(&y, vy, sizeof(vy) / sizeof(*vy))) {
            bi_cleanup(&x);
            puts("Error: bi_init(&y, ...)");
            return 1;
        }
        if (!bi_init(&z, vz, sizeof(vz) / sizeof(*vz))) {
            bi_cleanup(&y);
            bi_cleanup(&x);
            puts("Error: bi_init(&z, ...)");
            return 1;
        }
        bi_println("x = ", &x);
        bi_println("y = ", &y);
        bi_println("z = ", &z);
        printf("bi_compare(&x, &y): %s\n", bi_compare(&x, &y) == bi_below ? "ok" : "error");
        printf("bi_compare(&x, &z): %s\n", bi_compare(&x, &z) == bi_equal ? "ok" : "error");
        printf("bi_compare(&y, &z): %s\n", bi_compare(&y, &z) == bi_above ? "ok" : "error");
        printf("bi_mod(&y, %llu): %s\n", m, bi_mod(&y, m) == 182ull ? "ok" : "error");
        printf("bi_expmod(&y, &z, %llu): %s\n", m, bi_expmod(&y, &z, m) == 240ull ? "ok" : "error");
        bi_cleanup(&z);
        bi_cleanup(&y);
        bi_cleanup(&x);
    
    /* 
     * Tests 'bi_sum_init', 'bi_sum_over', 'bi_mul_init'
     *      and'bi_mul_over'.
     */

        BigInt a, b, c;
        memset(&a, 0, sizeof(a));
        memset(&b, 0, sizeof(b));
        memset(&c, 0, sizeof(c));
        
        //bi_block const va[] = {0xffffffffull, 0x000000001ull};
        //bi_block const vb[] = {0x00000002ull, 0x00000000ull};
        //bi_block const vc[] = {0x11110000ull, 0xffffffffull, 0xffffffffull, 0xffffffffull, 0x00000000ull};

        bi_block const va[] = {0x00001111ull, 0x00001111ull, 0x00001111ull, 0x00000000ull, 0x00000002ull, 0xffffffffull, 0x22220000ull};
        bi_block const vb[] = {0x00000000ull, 0xffffffffull, 0x11110000ull, 0x22220000ull, 0xffffffffull};
        bi_block const vc[] = {0x11110000ull, 0xffffffffull, 0xffffffffull, 0xffffffffull, 0x00000000ull};


        if (!bi_init(&a, va, sizeof(va) / sizeof(*va))) {
            puts("Error: bi_init(&x, ...)");
            return 1;
        }
        if (!bi_init(&b, vb, sizeof(vb) / sizeof(*vb))) {
            puts("Error: bi_init(&x, ...)");
            return 1;
        }
        if (!bi_init(&c, vc, sizeof(vc) / sizeof(*vc))) {
            puts("Error: bi_init(&x, ...)");
            return 1;
        }
        
        bi_println("a = ", &a);
        bi_println("b = ", &b);
        bi_println("c = ", &c);
        
        BigInt aPlusb;
        bi_mul_init(&a, &b, &aPlusb);
        bi_mul_over(&a, &b);
        bi_println("a * b = ", &aPlusb);
        bi_println("a over = ", &a);

    }

    return 0;
}

