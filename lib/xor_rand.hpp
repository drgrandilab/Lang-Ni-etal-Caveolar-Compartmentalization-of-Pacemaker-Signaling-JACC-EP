

#ifndef XOR_RAND_HPP
#define XOR_RAND_HPP

#define _USE_MATH_DEFINES
#include <cmath>
// #include <fstream>
// #include <cfloat>
#include <climits>
// #include <cstdlib>
#include <stdint.h>

class xor_rand
{
public:
	xor_rand(unsigned int seed = 0, unsigned int id = 0);
	~xor_rand() {};

	double gen_rand();
	void reset(unsigned int seed, unsigned int id);

	unsigned int gen_rand_uint();

	unsigned int xx, yy, zz, ww;

};






/*inline unsigned int xorshift(unsigned int *xx, unsigned int *yy, unsigned int *zz, unsigned int *ww)
{
	unsigned int t = (*xx ^ (*xx << 11)); *xx = *yy; *yy = *zz; *zz = *ww;
	return ( *ww = (*ww ^ (*ww >> 19)) ^ (t ^ (t >> 8)) ) / (double)(UINT_MAX);
}*/

#endif



// class xor_rand_256
// {
// public:
// 	xor_rand_256() {};
// 	~xor_rand_256() {};
// 	uint64_t s[4];
// 	uint64_t gen_rand(void) {
// 		const uint64_t result_starstar = rotl(s[1] * 5, 7) * 9;

// 		const uint64_t t = s[1] << 17;

// 		s[2] ^= s[0];
// 		s[3] ^= s[1];
// 		s[1] ^= s[2];
// 		s[0] ^= s[3];

// 		s[2] ^= t;

// 		s[3] = rotl(s[3], 45);

// 		return result_starstar;
// 	}


// 	/* This is the jump function for the generator. It is equivalent
// 	   to 2^128 calls to gen_rand(); it can be used to generate 2^128
// 	   non-overlapping subsequences for parallel computations. */

// 	void jump(void) {
// 		static const uint64_t JUMP[] = { 0x180ec6d33cfd0aba, 0xd5a61266f0c9392c, 0xa9582618e03fc9aa, 0x39abdc4529b1661c };

// 		uint64_t s0 = 0;
// 		uint64_t s1 = 0;
// 		uint64_t s2 = 0;
// 		uint64_t s3 = 0;
// 		for (int i = 0; i < sizeof JUMP / sizeof * JUMP; i++)
// 			for (int b = 0; b < 64; b++) {
// 				if (JUMP[i] & UINT64_C(1) << b) {
// 					s0 ^= s[0];
// 					s1 ^= s[1];
// 					s2 ^= s[2];
// 					s3 ^= s[3];
// 				}
// 				gen_rand();
// 			}

// 		s[0] = s0;
// 		s[1] = s1;
// 		s[2] = s2;
// 		s[3] = s3;
// 	}



// 	/* This is the long-jump function for the generator. It is equivalent to
// 	   2^192 calls to gen_rand(); it can be used to generate 2^64 starting points,
// 	   from each of which jump() will generate 2^64 non-overlapping
// 	   subsequences for parallel distributed computations. */

// 	void long_jump(void) {
// 		static const uint64_t LONG_JUMP[] = { 0x76e15d3efefdcbbf, 0xc5004e441c522fb3, 0x77710069854ee241, 0x39109bb02acbe635 };

// 		uint64_t s0 = 0;
// 		uint64_t s1 = 0;
// 		uint64_t s2 = 0;
// 		uint64_t s3 = 0;
// 		for (int i = 0; i < sizeof LONG_JUMP / sizeof * LONG_JUMP; i++)
// 			for (int b = 0; b < 64; b++) {
// 				if (LONG_JUMP[i] & UINT64_C(1) << b) {
// 					s0 ^= s[0];
// 					s1 ^= s[1];
// 					s2 ^= s[2];
// 					s3 ^= s[3];
// 				}
// 				gen_rand();
// 			}

// 		s[0] = s0;
// 		s[1] = s1;
// 		s[2] = s2;
// 		s[3] = s3;
// 	}

// };