//  models.c
//  xcssecovid
//
//  Created by Dr. Rolf Jansen on 2020-03-27.
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


#include <stdbool.h>
#include <math.h>

#include "numerics.h"
#include "models.h"


#pragma mark ••• Logistic Function •••

// https://en.wikipedia.org/wiki/Logistic_function
char *modelDescription_LF =
"# Model: Logistic Function\n"\
"#   y = a1/(1 + exp(-a0·(t - a2)))";

int initialValues_LF(ldouble t1, ldouble min, ldouble max, ldouble A[mpar], int f[mpar])
{
   if (isnan(A[0])) A[0] = 0.2L;
   if (isnan(A[1])) A[1] = 2.0L*max;
   if (isnan(A[2])) A[2] = t1 + 20.0L;

   int i, k = 0;
   for (i = 0; i < mpar; i++)
      if (f[i] != undefined) k++;
   if (k == 0)
      f[k++] = 0, f[k++] = 1, f[k++] = 2;
   return k;
}

void modelFunction_LF(ldouble t, ldouble *Y, ldouble A[mpar], bool init)
{
   *Y = A[1]/(1 + expl(-A[0]*(t - A[2])));
}


#pragma mark ••• Logistic Differential Equation •••

char *modelDescription_LDE =
"# Model: Logistic Differential Equation\n"\
"#   dy/dt = a0·y·(1 - y/a1) || y(a3) = a2";

int initialValues_LDE(ldouble t1, ldouble min, ldouble max, ldouble A[mpar], int f[mpar])
{
   if (isnan(A[0])) A[0] = 0.2L;
   if (isnan(A[1])) A[1] = 2.0L*max;
   if (isnan(A[2])) A[2] = min;
   if (isnan(A[3])) A[3] = t1;

   int i, k = 0;
   for (i = 0; i < mpar; i++)
      if (f[i] != undefined) k++;
   if (k == 0)
      f[k++] = 0, f[k++] = 1, f[k++] = 2;
   return k;
}

static void lde(ldouble t, ldouble *Y, ldouble *dY, ldouble A[mpar])
{
   *dY = A[0]**Y*(1 - *Y/A[1]);
}

void modelFunction_LDE(ldouble t, ldouble *Y, ldouble A[mpar], bool init)
{
   static ldouble t0;
   static ldouble Y0;

   if (init)
   {
      Y0 = A[2];
      ODEInt(1, A[3], t, &Y0, A, lde);
   }
   else
      ODEInt(1, t0, t, &Y0, A, lde);

   t0 = t;
   *Y = Y0;
}


#pragma mark ••• SI Differential Equations •••

char *modelDescription_SI =
"# Model: SI Differential Equations\n"\
"#   dy0/dt = -a0/a1·y0·y1 || y0(a3) = a1\n"\
"#   dy1/dt =  a0/a1·y0·y1 || y1(a3) = a2";

int initialValues_SI(ldouble t1, ldouble min, ldouble max, ldouble A[mpar], int f[mpar])
{
   if (isnan(A[0])) A[0] = 0.2;
   if (isnan(A[1])) A[1] = 2.0L*max;
   if (isnan(A[2])) A[2] = min;
   if (isnan(A[3])) A[3] = t1;

   int i, k = 0;
   for (i = 0; i < mpar; i++)
      if (f[i] != undefined) k++;
   if (k == 0)
      f[k++] = 0, f[k++] = 1, f[k++] = 2;
   return k;
}

static void sides(ldouble t, ldouble *Y, ldouble *dY, ldouble A[mpar])
{
   dY[0] = -A[0]/A[1]*Y[0]*Y[1];    // dS/dt
   dY[1] =  A[0]/A[1]*Y[0]*Y[1];    // dI/dt
}

void modelFunction_SI(ldouble t, ldouble *Y, ldouble A[mpar], bool init)
{
   static ldouble t0;
   static ldouble Y0[2];

   if (init)
   {
      Y0[0] = A[1];
      Y0[1] = A[2];
      ODEInt(2, A[3], t, Y0, A, sides);
   }
   else
      ODEInt(2, t0, t, Y0, A, sides);

   t0 = t;
   *Y = Y0[1];  // I - do curve fit of I(t)
}


#pragma mark ••• SIR Differential Equations •••

char *modelDescription_SIR =
"# Model: SIR Differential Equations\n"\
"#   dy0/dt = -a0/a1·y0·y1         || y0(a5) = a1\n"\
"#   dy1/dt =  a0/a1·y0·y1 - a3·y1 || y1(a5) = a2\n"\
"#   dy2/dt =  a3·y1               || y2(a5) = a4";

int initialValues_SIR(ldouble t1, ldouble min, ldouble max, ldouble A[mpar], int f[mpar])
{
   if (isnan(A[0])) A[0] = 0.6;
   if (isnan(A[1])) A[1] = 2.0L*max;
   if (isnan(A[2])) A[2] = 1.5L*min;
   if (isnan(A[3])) A[3] = 0.4L;
   if (isnan(A[4])) A[4] = 0.5L*min;
   if (isnan(A[5])) A[5] = t1;

   int i, k = 0;
   for (i = 0; i < mpar; i++)
      if (f[i] != undefined) k++;
   if (k == 0)
      f[k++] = 0, f[k++] = 1, f[k++] = 2;
   return k;
}

static void sirdes(ldouble t, ldouble *Y, ldouble *dY, ldouble A[mpar])
{
   dY[0] = -A[0]/A[1]*Y[0]*Y[1];                // dS/dt
   dY[1] =  A[0]/A[1]*Y[0]*Y[1] - A[3]*Y[1];    // dI/dt
   dY[2] =  A[3]*Y[1];                          // dR/dt
}

void modelFunction_SIR(ldouble t, ldouble *Y, ldouble A[mpar], bool init)
{
   static ldouble t0;
   static ldouble Y0[3];

   if (init)
   {
      Y0[0] = A[1];
      Y0[1] = A[2];
      Y0[2] = A[4];
      ODEInt(3, A[5], t, Y0, A, sirdes);
   }
   else
      ODEInt(3, t0, t, Y0, A, sirdes);

   t0 = t;
   *Y = Y0[2];  // R - do curve fit of R(t)
}


#pragma mark ••• Shifted Error Function •••

// generic form: https://en.wikipedia.org/wiki/Error_function
// the error function is rotational symmetrical at its inflection point at [0, 0].
// we need to shift it up by one, in order to make it all positive.
char *modelDescription_ERF =
"# Model: Shifted Error Function\n"\
"#   y = 0.5·a1·(1 + erf(a0·(t - a2)))";

int initialValues_ERF(ldouble t1, ldouble min, ldouble max, ldouble A[mpar], int f[mpar])
{
   if (isnan(A[0])) A[0] = 0.1L;
   if (isnan(A[1])) A[1] = 2.0L*max;
   if (isnan(A[2])) A[2] = t1 + 20.0L;

   int i, k = 0;
   for (i = 0; i < mpar; i++)
      if (f[i] != undefined) k++;
   if (k == 0)
      f[k++] = 0, f[k++] = 1, f[k++] = 2;
   return k;
}

void modelFunction_ERF(ldouble t, ldouble *Y, ldouble A[mpar], bool init)
{
   *Y = 0.5*A[1]*(1 + erfl(A[0]*(t - A[2])));
}


#pragma mark ••• Generalized Logistic Function •••

// https://en.wikipedia.org/wiki/Generalised_logistic_function
char *modelDescription_GLF =
"# Model: Generalised Logistic Function\n"\
"#   y = a1/(1 + exp(-a0·a3·(t - a2)))^(1/a3)";

int initialValues_GLF(ldouble t1, ldouble min, ldouble max, ldouble A[mpar], int f[mpar])
{
   if (isnan(A[0])) A[0] = 0.2L;
   if (isnan(A[1])) A[1] = 2.0L*max;
   if (isnan(A[2])) A[2] = t1 + 20.0L;
   if (isnan(A[3])) A[3] = 0.25L;

   int i, k = 0;
   for (i = 0; i < mpar; i++)
      if (f[i] != undefined) k++;
   if (k == 0)
      f[k++] = 0, f[k++] = 1, f[k++] = 2;
   return k;
}

void modelFunction_GLF(ldouble t, ldouble *Y, ldouble A[mpar], bool init)
{
   *Y = A[1]/powl(1 + expl(-A[0]*A[3]*(t - A[2])), 1/A[3]);
}
