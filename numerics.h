//  numerics.h
//  xcssecovid
//
//  Created by Dr. Rolf Jansen on 2020-03-26.
//  Copyright © 2020 Dr. Rolf Jansen. All rights reserved.
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice,
//     this list of conditions and the following disclaimer.
//
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
//  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
//  INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
//  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
//  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
//  OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
//  OF THE POSSIBILITY OF SUCH DAMAGE.


typedef long double ldouble;


#pragma mark ••• Bulirsch–Stoer Differential Equation Solver •••

typedef enum
{
   deq_noError,
   deq_tolError,
   deq_stiffError,
   deq_stepsError
} deq_error;

typedef enum
{
   undefined = -1,
   inactive  =  0,
   active    =  1
} tris;

static inline ldouble sqrl(ldouble x)
{
   return x*x;
}

typedef void (*odeset)(ldouble t, ldouble *Y, ldouble *dY, ldouble *A);
deq_error ODEInt(int m, ldouble t0, ldouble t, ldouble *Y, ldouble *A, odeset equations);


#pragma mark ••• Curve Fitting by Levenberg–Marquardt least squares minimization •••

#define mpar 10

typedef int (*initvals)(ldouble t1, ldouble min, ldouble max, ldouble *A, int *f);
typedef int (*function)(ldouble t, ldouble *Y, ldouble *A, bool init);
ldouble curveFit(int n, ldouble *T, ldouble *Y,
                 int k, ldouble *A, ldouble *dA, int *f, function curve);


#pragma mark ••• LU Decomposition •••

ldouble LUdecomposition(int m, ldouble **A, ldouble **LU, int *idx);
void LUbacksubstitution(int m, ldouble **LU, int *idx, ldouble *B, ldouble *X);
void LUrefinment(int m, ldouble **A, ldouble **LU, int *idx, ldouble *B, ldouble *X);
void LUinversion(int m, ldouble **A, ldouble **LU, int *idx);
