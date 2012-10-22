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

#include "random.h"
#include "math.h"

uint32_t mt[MERS_N];                    // State vector
int mti;                                // Index into mt
uint32_t LastInterval = 0;              // Last interval length for IRandomX
uint32_t RLimit;                        // Rejection limit used by IRandomX

void init_mersenne(int seed) 
{
    const uint32_t factor = 1812433253UL;
    mt[0]= seed;
    for (mti = 1; mti < MERS_N; mti++) 
    {
        mt[mti] = (factor * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
    }
}

void set_seed(int seed) 
{
    init_mersenne(seed);

    // Randomize some more
    int i;
    for (i = 0; i < 37; i++) get_random_bits();
}

void set_seed_array(int const seeds[], int nseeds) 
{
    int i, j, k;

    // Initialize
    init_mersenne(19650218);

    if (nseeds <= 0) return;

    // Randomize mt[] using whole seeds[] array
    i = 1;  j = 0;
    k = (MERS_N > nseeds ? MERS_N : nseeds);
    for (; k; k--) 
    {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL)) + (uint32_t)seeds[j] + j;
        i++; j++;
        if (i >= MERS_N) {mt[0] = mt[MERS_N-1]; i=1;}
        if (j >= nseeds) j=0;
    }
    for (k = MERS_N-1; k; k--) 
    {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL)) - i;
        if (++i >= MERS_N) { mt[0] = mt[MERS_N-1]; i=1; }
    }
    mt[0] = 0x80000000UL;  // MSB is 1; assuring non-zero initial array

    // Randomize some more
    mti = 0;
    for (i = 0; i <= MERS_N; i++) get_random_bits();
}

uint32_t get_random_bits() 
{
    // Generate 32 random bits
    uint32_t y;

    if (mti >= MERS_N) 
    {
        // Generate MERS_N words at one time
        const uint32_t LOWER_MASK = (1LU << MERS_R) - 1;       // Lower MERS_R bits
        const uint32_t UPPER_MASK = 0xFFFFFFFF << MERS_R;      // Upper (32 - MERS_R) bits
        static const uint32_t mag01[2] = {0, MERS_A};

        int kk;
        for (kk=0; kk < MERS_N-MERS_M; kk++) 
        {    
            y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
            mt[kk] = mt[kk+MERS_M] ^ (y >> 1) ^ mag01[y & 1];
        }

        for (; kk < MERS_N-1; kk++) 
        {    
            y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
            mt[kk] = mt[kk+(MERS_M-MERS_N)] ^ (y >> 1) ^ mag01[y & 1];
        }

        y = (mt[MERS_N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
        mt[MERS_N-1] = mt[MERS_M-1] ^ (y >> 1) ^ mag01[y & 1];
        mti = 0;
    }
    y = mt[mti++];

    // Tempering (May be omitted):
    y  ^=  y >> MERS_U;
    y ^= (y << MERS_S) & MERS_B;
    y ^= (y << MERS_T) & MERS_C;
    y ^=  y >> MERS_L;

    return y;
}

double get_frandom() 
{
    // Multiply by 2^(-32)
    return (double)get_random_bits() * (1./(65536.*65536.));
}

int get_irandom(int min, int max)
{
    if (max <= min) 
    {
        if (max == min) return min; else return 0x80000000;
    }
    // Multiply interval with random and truncate
    int r = (int)((double)(uint32_t)(max - min + 1) * get_frandom()  + min); 
    if (r > max) r = max;
    return r;
}

int random_g(double *y1,double *y2){
    double w,x1,x2;
    
    do {
        x1=2.0*get_frandom() -1.0;
        x2=2.0*get_frandom() -1.0;
        w=x1*x1+x2*x2;
    } while (w>= 1.0);
    w=sqrt( ( -2.0*log(w))/ w );
    *y1=x1*w;
    *y2=x2*w;

    return 0;
}

int get_irandomx(int min, int max) 
{
    if (max <= min) 
    {
        if (max == min) return min; else return 0x80000000;
    }

#ifdef  INT64_SUPPORTED
    // 64 bit integers available. Use multiply and shift method
    uint32_t interval;                    // Length of interval
    uint64_t longran;                     // Random bits * interval
    uint32_t iran;                        // Longran / 2^32
    uint32_t remainder;                   // Longran % 2^32

    interval = (uint32_t)(max - min + 1);
    if (interval != LastInterval) 
    {
        // Interval length has changed. Must calculate rejection limit
        // Reject when remainder >= 2^32 / interval * interval
        // RLimit will be 0 if interval is a power of 2. No rejection then
        RLimit = (uint32_t)(((uint64_t)1 << 32) / interval) * interval - 1;
        LastInterval = interval;
    }
    do 
    { // Rejection loop
        longran  = (uint64_t)get_random_bits() * interval;
        iran = (uint32_t)(longran >> 32);
        remainder = (uint32_t)longran;
    } while (remainder > RLimit);
    // Convert back to signed and return result
    return (int32_t)iran + min;

#else
    // 64 bit integers not available. Use modulo method
    uint32_t interval;                    // Length of interval
    uint32_t bran;                        // Random bits
    uint32_t iran;                        // bran / interval
    uint32_t remainder;                   // bran % interval

    interval = uint32_t(max - min + 1);
    if (interval != LastInterval) 
    {
        // Interval length has changed. Must calculate rejection limit
        // Reject when iran = 2^32 / interval
        // We can't make 2^32 so we use 2^32-1 and correct afterwards
        RLimit = (uint32_t)0xFFFFFFFF / interval;
        if ((uint32_t)0xFFFFFFFF % interval == interval - 1) RLimit++;
    }
    do 
    { // Rejection loop
        bran = get_random_bits();
        iran = bran / interval;
        remainder = bran % interval;
    } while (iran >= RLimit);
    // Convert back to signed and return result
    return (int32_t)remainder + min;

#endif
}
