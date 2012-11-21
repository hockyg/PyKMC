/* *
 * The C Protein Folding Library.
 * Copyright (C) 2009 Andres Colubri.
 * Contact: andres.colubri 'AT' gmail.com
 *
 * This library was written at the Institute of Biophysical Dynamics at the University of Chicago.
 * Gordon Center for Integrated Science, W101. 929 East 57th Street, Chicago, IL 60637, USA.
 * Homepage: http://ibd.uchicago.edu/
 * 
 * Permission to use, copy, modify, and distribute this software and its
 * documentation with or without modifications and for any purpose and
 * without fee is hereby granted, provided that any copyright notices
 * appear in all copies and that both those copyright notices and this
 * permission notice appear in supporting documentation, and that the
 * names of the contributors or copyright holders not be used in
 * advertising or publicity pertaining to distribution of the software
 * without specific prior permission.
 * 
 * THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
 * WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
 * CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
 * OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
 * OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
 * OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
 * OR PERFORMANCE OF THIS SOFTWARE.
 * 
 */

/**
 * Implementation of a Mersenne Twister Random Number Generator. 
 * This random number generator is described in the article by
 * M. Matsumoto & T. Nishimura, in:
 * ACM Transactions on Modeling and Computer Simulation,
 * vol. 8, no. 1, 1998, pp. 3-30.
 * I has excellent statistical properties and it is very well suited for Monte Carlo simulations
 * This code has been ported to C from the original C++ code by Agner Fog. This code can be
 * found at http://www.agner.org.
 * For more information on the Mersenne Twister, check prof. M. Matsumoto's website:
 * http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
 *
 */

#ifndef __RANDOM_H__
#define __RANDOM_H__

#define MAXINT 4294967296 // 65536*65536 = 2^32
#define INVMAXINT 0.00000000023283064365386962890625 // 2^-32
#define DOUBLESHIFT 0.000000000116415321826934814453125 // 2^-33

// Define integer types with known size: int32_t, uint32_t, int64_t, uint64_t.
// If this doesn't work then insert compiler-specific definitions here:
#if defined(__GNUC__)
    // Compilers supporting C99 or C++0x have inttypes.h defining these integer types
    #include <inttypes.h>
    #define INT64_SUPPORTED // Remove this if the compiler doesn't support 64-bit integers
#elif defined(_WIN16) || defined(__MSDOS__) || defined(_MSDOS) 
     // 16 bit systems use long int for 32 bit integer
    typedef   signed long int int32_t;
    typedef unsigned long int uint32_t;
#elif defined(_MSC_VER)
    // Microsoft have their own definition
    typedef   signed __int32  int32_t;
    typedef unsigned __int32 uint32_t;
    typedef   signed __int64  int64_t;
    typedef unsigned __int64 uint64_t;
    #define INT64_SUPPORTED // Remove this if the compiler doesn't support 64-bit integers
#else
    // This works with most compilers
    typedef signed int          int32_t;
    typedef unsigned int       uint32_t;
    typedef long long           int64_t;
    typedef unsigned long long uint64_t;
    #define INT64_SUPPORTED // Remove this if the compiler doesn't support 64-bit integers
#endif

// #define MT11213A // Uses the Mersenne Twister 11213A (otherwise 19937 is used).

// Choose which version of Mersenne Twister you want:
#ifdef MT11213A
    // Define constants for type MT11213A:
    #define MERS_N   351
    #define MERS_M   175
    #define MERS_R   19
    #define MERS_U   11
    #define MERS_S   7
    #define MERS_T   15
    #define MERS_L   17
    #define MERS_A   0xE4BD75F5
    #define MERS_B   0x655E5280
    #define MERS_C   0xFFD58000
#else    
    // or constants for type MT19937:
    #define MERS_N   624
    #define MERS_M   397
    #define MERS_R   31
    #define MERS_U   11
    #define MERS_S   7
    #define MERS_T   15
    #define MERS_L   18
    #define MERS_A   0x9908B0DF
    #define MERS_B   0x9D2C5680
    #define MERS_C   0xEFC60000
#endif

/**
 * Basic initialization procedure.
 * @param seed int
 */
void init_mersenne(int seed);

/**
 * Initialize and seed generator.
 * @param seed int
 */

void set_seed(int seed);

/**
 * Seed by more than 32 bits,
 * @param seeds[] int const
 * @param nseeds int
 */
void set_seed_array(int const seeds[], int nseeds);

/**
 * Generates 32 random bits,
 */
uint32_t get_random_bits();

/**
 * Generates random float number in the interval 0 <= x < 1
 */
double get_frandom();

/**
 * Generates random double number in the interval 0 < x < 1
 *  by shifting value from get_frandom by 2^-33
 */
double get_frandom_2();

/**
 * Generates random integer in the interval min <= x <= max
 * Relative error on frequencies < 2^-32
 * @param min int
 * @param max int
 */
int get_irandom(int min, int max);

/**
 * Output random integer in the interval min <= x <= max
 * Each output value has exactly the same probability.
 * This is obtained by rejecting certain bit values so that the number
 * of possible bit values is divisible by the interval length
 * @param min int
 * @param max int
 */
int get_irandomx(int min, int max);     // Output random integer, exact

/**
 * Output two random floats with mean 0 and std dev 1 given
 * two floats from a uniform distribution in [0,1]
 * uses the Box-Muller transformation
 * 
 * @param y1 *float
 * @param y2 *float
 * @param x1 float
 * @param x2 float
 */
int random_g(double *y1,double *y2 );




#endif
