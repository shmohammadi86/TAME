#ifndef _MIC_UTILS_H
#define _MIC_UTILS_H

#include <omp.h>
#include <cmath>
#include <iostream>
#include <immintrin.h>
using namespace std;

#define _NEW_OPT_ 
#define _ACTIVATE_MIC_CODE_

#define _ACTIVATE_TEST_CODE_

#ifdef __MIC__

#define SIMDMASKTYPE __mmask16
#define SIMDFPTYPE __m512
#define SIMDINTTYPE __m512i
#define _MM_LOAD(A)  _mm512_extload_ps((__m512*)(A), _MM_UPCONV_PS_NONE, _MM_BROADCAST32_NONE, _MM_HINT_NONE)
#define _MM_LOAD_INT(A)  _mm512_extload_epi32((__m512i*)(A), _MM_UPCONV_EPI32_NONE, _MM_BROADCAST32_NONE, _MM_HINT_NONE)
#define _MM_LOAD_MASK_INT(v_old, mask, A)  _mm512_mask_extload_epi32(v_old, mask, (__m512i*)(A), _MM_UPCONV_EPI32_NONE, _MM_BROADCAST32_NONE, _MM_HINT_NONE)
#define _MM_STORE(A, B)  _mm512_extstore_ps((__m512*)(A),   B, _MM_DOWNCONV_PS_NONE, _MM_HINT_NONE);
#define _MM_STORE_MASK(A, B, mask)  _mm512_mask_extstore_ps((__m512*)(A), mask,  B, _MM_DOWNCONV_PS_NONE, _MM_HINT_NONE);
#define _MM_SET(x) _mm512_set1_ps(x)
#define _MM_SET_INT(x) _mm512_set1_epi32(x)

#define _MM_STORE_INT(A, B)  _mm512_extstore_epi32((__m512i*)(A),   B, _MM_DOWNCONV_EPI32_NONE, _MM_HINT_NONE);
#define _MM_STORE_INT_MASK(A, B, mask)  _mm512_mask_extstore_epi32((__m512i*)(A), mask,  B, _MM_DOWNCONV_EPI32_NONE, _MM_HINT_NONE);
#define _MM_CMP_LE(x,y) _mm512_cmple_ps_mask(x,y)
#define _MM_CMP_LE_MASK(m,x,y) _mm512_mask_cmple_ps_mask(m,x,y)
#define _MM_CMP_LT(x,y) _mm512_cmplt_ps_mask(x,y)
#define _MM_CMP_LT_MASK(m,x,y) _mm512_mask_cmplt_ps_mask(m,x,y)
#define _MM_CMP_EQ_MASK(m,x,y) _mm512_mask_cmpeq_ps_mask(m,x,y)
#define _MM_CMP_LE_MASK_INT(m,x,y) _mm512_mask_cmple_epi32_mask(m,x,y)
#define _MM_CMP_LT_MASK_INT(m,x,y) _mm512_mask_cmplt_epi32_mask(m,x,y)
inline SIMDFPTYPE _MM_LOADU(float *X)
{
	SIMDFPTYPE vdst = _mm512_undefined_ps();
	vdst = _mm512_extloadunpacklo_ps(vdst, (float *)(X), _MM_UPCONV_PS_NONE, _MM_HINT_NONE);
	vdst = _mm512_extloadunpackhi_ps(vdst, (float *)(X+16), _MM_UPCONV_PS_NONE, _MM_HINT_NONE);
	return vdst;
}

inline SIMDFPTYPE _MM_LOADU_MASK(float *X, SIMDMASKTYPE mask, SIMDFPTYPE xmm0)
{
	SIMDFPTYPE vdst = _mm512_mask_extloadunpacklo_ps(xmm0, mask, (float *)(X), _MM_UPCONV_PS_NONE, _MM_HINT_NONE);
	vdst = _mm512_mask_extloadunpackhi_ps(vdst, mask, (float *)(X+16), _MM_UPCONV_PS_NONE, _MM_HINT_NONE);
	return vdst;
}
inline SIMDINTTYPE _MM_LOADU_INT(int *X)
{
	SIMDINTTYPE vdst = _mm512_undefined_epi32();
	vdst = _mm512_extloadunpacklo_epi32(vdst, (int *)(X), _MM_UPCONV_EPI32_NONE, _MM_HINT_NONE);
	vdst = _mm512_extloadunpackhi_epi32(vdst, (int *)(X+16), _MM_UPCONV_EPI32_NONE, _MM_HINT_NONE);
	return vdst;
}

inline SIMDINTTYPE _MM_LOADU_MASK_INT(int *X, SIMDMASKTYPE mask, SIMDINTTYPE xmm0)
{
	SIMDINTTYPE vdst = _mm512_mask_extloadunpacklo_epi32(xmm0, mask, (int *)(X), _MM_UPCONV_EPI32_NONE, _MM_HINT_NONE);
	vdst = _mm512_mask_extloadunpackhi_epi32(vdst, mask, (int *)(X+16), _MM_UPCONV_EPI32_NONE, _MM_HINT_NONE);
	return vdst;
}

inline void _MM_STOREU_INT(int *X, SIMDINTTYPE xmm0)
{
	_mm512_extpackstorelo_epi32((X), xmm0, _MM_DOWNCONV_EPI32_NONE, _MM_HINT_NONE);
	_mm512_extpackstorehi_epi32((X+16), xmm0, _MM_DOWNCONV_EPI32_NONE, _MM_HINT_NONE);
}

inline void _MM_STOREU_MASK_INT(int *X, SIMDMASKTYPE mask, SIMDINTTYPE xmm0)
{
	_mm512_mask_extpackstorelo_epi32((X), mask, xmm0, _MM_DOWNCONV_EPI32_NONE, _MM_HINT_NONE);
	_mm512_mask_extpackstorehi_epi32((X+16), mask, xmm0, _MM_DOWNCONV_EPI32_NONE, _MM_HINT_NONE);
}

#endif 





#endif //BMATCHING_H
