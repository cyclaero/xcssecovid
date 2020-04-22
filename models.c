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

int modelFunction_LF(ldouble t, ldouble *Y, ldouble A[mpar], bool init)
{
   *Y = A[1]/(1 + expl(-A[0]*(t - A[2])));
   return 0;
}


#pragma mark ••• Logistic Differential Equation •••

// https://en.wikipedia.org/wiki/Logistic_function#Logistic_differential_equation
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

int modelFunction_LDE(ldouble t, ldouble *Y, ldouble A[mpar], bool init)
{
   int rc;

   static ldouble t0;
   static ldouble Y0;

   if (init)
   {
      Y0 = A[2];
      rc = ODEInt(1, A[3], t, &Y0, A, lde);
   }
   else
      rc = ODEInt(1, t0, t, &Y0, A, lde);

   t0 = t;
   *Y = Y0;

   return rc;
}


#pragma mark ••• SI Differential Equations •••

// https://de.wikipedia.org/wiki/SI-Modell
char *modelDescription_SI =
"# Model: SI Differential Equations\n"\
"#   dy0/dt = -a0/a1·y0·y1 || y0(a3) = a1-a2\n"\
"#   dy1/dt =  a0/a1·y0·y1 || y1(a3) = a2";

int initialValues_SI(ldouble t1, ldouble min, ldouble max, ldouble A[mpar], int f[mpar])
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

static void sides(ldouble t, ldouble *Y, ldouble *dY, ldouble A[mpar])
{
   dY[0] = -A[0]/A[1]*Y[0]*Y[1];    // dS/dt
   dY[1] =  A[0]/A[1]*Y[0]*Y[1];    // dI/dt
}

int modelFunction_SI(ldouble t, ldouble *Y, ldouble A[mpar], bool init)
{
   int rc;

   static ldouble t0;
   static ldouble Y0[2];

   if (init)
   {
      Y0[0] = A[1] - A[2];
      Y0[1] = A[2];
      rc = ODEInt(2, A[3], t, Y0, A, sides);
   }
   else
      rc = ODEInt(2, t0, t, Y0, A, sides);

   t0 = t;
   *Y = Y0[1];  // I - curve fit of I(t)

   return rc;
}


#pragma mark ••• SIR Differential Equations •••

// https://en.wikipedia.org/wiki/Mathematical_modelling_of_infectious_disease#The_SIR_model
char *modelDescription_SIR =
"# Model: SIR Differential Equations\n"\
"# S dy0/dt = -a0/a1·y0·y1         || y0(a5) = a1-a2-a4\n"\
"# I dy1/dt =  a0/a1·y0·y1 - a3·y1 || y1(a5) = a2\n"\
"# R dy2/dt =  a3·y1               || y2(a5) = a4 <- a3·a2";

int initialValues_SIR(ldouble t1, ldouble min, ldouble max, ldouble A[mpar], int f[mpar])
{
   if (isnan(A[0])) A[0] = 0.6L;
   if (isnan(A[1])) A[1] = 2.0L*max;
   if (isnan(A[2])) A[2] = 2.5L*min;
   if (isnan(A[3])) A[3] = 0.4L;
   if (isnan(A[4])) A[4] = A[3]*A[2];
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

int modelFunction_SIR(ldouble t, ldouble *Y, ldouble A[mpar], bool init)
{
   int rc;

   static ldouble t0;
   static ldouble Y0[3];

   if (init)
   {
      Y0[0] = A[1] - A[2] - A[4];
      Y0[1] = A[2];
      Y0[2] = A[4];
      rc = ODEInt(3, A[5], t, Y0, A, sirdes);
   }
   else
      rc = ODEInt(3, t0, t, Y0, A, sirdes);

   t0 = t;
   *Y = Y0[2];  // R - curve fit of R(t)

   return rc;
}


#pragma mark ••• SEIR Differential Equations •••

// https://www.idmod.org/docs/hiv/model-seir.html
char *modelDescription_SEIR =
"# Model: SEIR Differential Equations\n"\
"# S  dy0/dt = -a0/a1·y0·y1                  || y0(a7) = a1-a2-a5-a6\n"\
"# E  dy1/dt =  a0/a1·y0·y1 - a3·y1          || y1(a7) = a2\n"\
"# I  dy2/dt =  a3·y1 - a4·y2                || y2(a7) = a5 <- a2·a3\n"\
"# R  dy3/dt =  a4·y2                        || y3(a7) = a6 <- a2·a3·a4";

int initialValues_SEIR(ldouble t1, ldouble min, ldouble max, ldouble A[mpar], int f[mpar])
{
   if (isnan(A[0])) A[0] =  0.6L;            // beta  - infection rate
   if (isnan(A[1])) A[1] =  2.0L*max;        // population
   if (isnan(A[2])) A[2] = 10.0L*min;        // total number of exposed individuals at t1
   if (isnan(A[3])) A[3] =  0.4L;            // sigma - incubation rate (2.5 d until an infected individual becomes infectuous)
   if (isnan(A[4])) A[4] =  0.2L;            // gamma - removal rate (more 5 d until the infectuous individual can be removed from the chain of infection)
   if (isnan(A[5])) A[5] = A[2]*A[3];        // I(t1) boundary value at t1
   if (isnan(A[6])) A[6] = A[2]*A[3]*A[4];   // R(t1) boundary value at t1
   if (isnan(A[7])) A[7] = t1;

   int i, k = 0;
   for (i = 0; i < mpar; i++)
      if (f[i] != undefined) k++;
   if (k == 0)
      f[k++] = 0, f[k++] = 1, f[k++] = 3;
   return k;
}

static void seirdes(ldouble t, ldouble *Y, ldouble *dY, ldouble A[mpar])
{
   dY[0] = -A[0]/A[1]*Y[0]*Y[1];             // dS/dt
   dY[1] =  A[0]/A[1]*Y[0]*Y[1] - A[3]*Y[1]; // dE/dt
   dY[2] =  A[3]*Y[1] - A[4]*Y[2];           // dI/dt
   dY[3] =  A[4]*Y[2];                       // dR/dt
}

int modelFunction_SEIR(ldouble t, ldouble *Y, ldouble A[mpar], bool init)
{
   int rc;

   static ldouble t0;
   static ldouble Y0[4];

   if (init)
   {
      Y0[0] = A[1] - A[2] - A[5] - A[6];
      Y0[1] = A[2];
      Y0[2] = A[5];
      Y0[3] = A[6];
      rc = ODEInt(4, A[7], t, Y0, A, seirdes);
   }
   else
      rc = ODEInt(4, t0, t, Y0, A, seirdes);

   t0 = t;
   *Y = Y0[3];  // R - curve fit of R(t)

   return rc;
}


#pragma mark ••• SIRX Differential Equations •••

// http://rocs.hu-berlin.de/corona/docs/forecast/model/
char *modelDescription_SIRX =
"# Model: SIRX Differential Equations\n"\
"# S  dy0/dt = -a0/a1·y0·y1 - a4·y0                || y0(a8) = a1-a2-a6-a7\n"\
"# I  dy1/dt =  a0/a1·y0·y1 - (a3 + a4 + a5)·y1    || y1(a8) = a2\n"\
"# R  dy2/dt =  a3·y1 + a4·y0                      || y2(a8) = a6 <- a3·a2 + a4·a1\n"\
"# X  dy3/dt =  (a4 + a5)·y1                       || y3(a8) = a7 <- (a4 + a5)·a2";

int initialValues_SIRX(ldouble t1, ldouble min, ldouble max, ldouble A[mpar], int f[mpar])
{
   if (isnan(A[0])) A[0] =  0.6L;                  // alpha  - infection rate
   if (isnan(A[1])) A[1] = 80.0e6L;                // population
   if (isnan(A[2])) A[2] = 10.0L*min;              // total number of infectious individuals at t1
   if (isnan(A[3])) A[3] =  0.38L;                 // beta   - removal rate
   if (isnan(A[4])) A[4] =  0.001L;                // kappa0 - general containment rate (all)
   if (isnan(A[5])) A[5] =  0.0005L;               // kappa  - quarantine rate (infected)
   if (isnan(A[6])) A[6] = A[3]*A[2] + A[4]*A[1];  // R(t1) boundary value at t1
   if (isnan(A[7])) A[7] = (A[4] + A[5])*A[2];     // X(t1) boundary value at t1
   if (isnan(A[8])) A[8] = t1;

   int i, k = 0;
   for (i = 0; i < mpar; i++)
      if (f[i] != undefined) k++;
   if (k == 0)
      f[k++] = 2, f[k++] = 4, f[k++] = 5;
   return k;
}

static void sirxdes(ldouble t, ldouble *Y, ldouble *dY, ldouble A[mpar])
{
   dY[0] = -A[0]/A[1]*Y[0]*Y[1] - A[4]*Y[0];                   // dS/dt
   dY[1] =  A[0]/A[1]*Y[0]*Y[1] - (A[3] + A[4] + A[5])*Y[1];   // dI/dt
   dY[2] =  A[3]*Y[1] + A[4]*Y[0];                             // dR/dt
   dY[3] =  (A[4] + A[5])*Y[1];                                // dX/dt
}

int modelFunction_SIRX(ldouble t, ldouble *Y, ldouble A[mpar], bool init)
{
   int rc;

   static ldouble t0;
   static ldouble Y0[4];

   if (init)
   {
      Y0[0] = A[1] - A[2] - A[6] - A[7];
      Y0[1] = A[2];
      Y0[2] = A[6];
      Y0[3] = A[7];
      rc = ODEInt(4, A[8], t, Y0, A, sirxdes);
   }
   else
      rc = ODEInt(4, t0, t, Y0, A, sirxdes);

   t0 = t;
   *Y = Y0[3];  // X - curve fit of X(t)

   return rc;
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

int modelFunction_ERF(ldouble t, ldouble *Y, ldouble A[mpar], bool init)
{
   *Y = 0.5*A[1]*(1 + erfl(A[0]*(t - A[2])));
   return 0;
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
   if (isnan(A[3])) A[3] = 0.333333L;

   int i, k = 0;
   for (i = 0; i < mpar; i++)
      if (f[i] != undefined) k++;
   if (k == 0)
      f[k++] = 0, f[k++] = 1, f[k++] = 2;
   return k;
}

int modelFunction_GLF(ldouble t, ldouble *Y, ldouble A[mpar], bool init)
{
   *Y = A[1]/powl(1 + expl(-A[0]*A[3]*(t - A[2])), 1/A[3]);
   return 0;
}
