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

int modelFunction_LF(ldouble t, ldouble *Y, ldouble A[mpar], int f[mpar], bool fit, bool init)
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

int modelFunction_LDE(ldouble t, ldouble *Y, ldouble A[mpar], int f[mpar], bool fit, bool init)
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

int modelFunction_SI(ldouble t, ldouble *Y, ldouble A[mpar], int f[mpar], bool fit, bool init)
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

   if (fit)
      *Y = Y0[1];  // I - curve fit of I(t)
   else
      Y[0] = Y0[0],
      Y[1] = Y0[1],
      Y[2] = NAN;

   return rc;
}


#pragma mark ••• SIR Differential Equations •••

// https://en.wikipedia.org/wiki/Mathematical_modelling_of_infectious_disease#The_SIR_model
char *modelDescription_SIR =
"# Model: SIR Differential Equations\n"\
"# S  dy0/dt = -a0/a1·y0·y1            || y0(a5) = a1-a2-a4\n"\
"# I  dy1/dt =  a0/a1·y0·y1 - a3·y1    || y1(a5) = a2 <- a4/a3\n"\
"# R  dy2/dt =  a3·y1                  || y2(a5) = a4";

int initialValues_SIR(ldouble t1, ldouble min, ldouble max, ldouble A[mpar], int f[mpar])
{
   if (isnan(A[0])) A[0] = 0.4L;
   if (isnan(A[1])) A[1] = 2.0L*max;
   if (isnan(A[3])) A[3] = 0.08L;
   if (isnan(A[4])) A[4] = min;
   if (isnan(A[5])) A[5] = t1;

   if (isnan(A[2])) A[2] = A[4]/A[3];

   int i, k = 0;
   for (i = 0; i < mpar; i++)
      if (f[i] != undefined) k++;
   if (k == 0)
      f[k++] = 0, f[k++] = 1, f[k++] = 3;

   return k;
}

static void sirdes(ldouble t, ldouble *Y, ldouble *dY, ldouble A[mpar])
{
   dY[0] = -A[0]/A[1]*Y[0]*Y[1];                // dS/dt
   dY[1] =  A[0]/A[1]*Y[0]*Y[1] - A[3]*Y[1];    // dI/dt
   dY[2] =  A[3]*Y[1];                          // dR/dt
}

int modelFunction_SIR(ldouble t, ldouble *Y, ldouble A[mpar], int f[mpar], bool fit, bool init)
{
   int rc;

   static ldouble t0;
   static ldouble Y0[3];

   if (init)
   {
      bool isvar2 = false;

      for (int i = 0; i < mpar && f[i] != undefined; i++)
         if (f[i] == 2)
            isvar2 = true;

      if (!isvar2) A[2] = A[4]/A[3];

      Y0[0] = A[1] - A[2] - A[4];
      Y0[1] = A[2];
      Y0[2] = A[4];
      rc = ODEInt(3, A[5], t, Y0, A, sirdes);
   }
   else
      rc = ODEInt(3, t0, t, Y0, A, sirdes);

   t0 = t;

   if (fit)
      *Y = Y0[2];  // R - curve fit of R(t)
   else
      Y[0] = Y0[0],
      Y[1] = Y0[1],
      Y[2] = Y0[2],
      Y[3] = NAN;

   return rc;
}


#pragma mark ••• SEIR Differential Equations •••

// https://www.idmod.org/docs/hiv/model-seir.html
char *modelDescription_SEIR =
"# Model: SEIR Differential Equations\n"\
"# S  dy0/dt = -a0/a1·y0·y2 + a8/y0    || y0(a7) = a1-a2-a5-a6\n"\
"# E  dy1/dt =  a0/a1·y0·y2 - a3·y1    || y1(a7) = a2 <- a6/a4/a3\n"\
"# I  dy2/dt =  a3·y1 - a4·y2          || y2(a7) = a5 <- (1 - a4)·a3·a2\n"\
"# R  dy3/dt =  a4·y2                  || y3(a7) = a6";

int initialValues_SEIR(ldouble t1, ldouble min, ldouble max, ldouble A[mpar], int f[mpar])
{
   if (isnan(A[0])) A[0] =  0.5L;                     // beta  - infection rate
   if (isnan(A[1])) A[1] =  max;                      // VSP   - Virtual Susceptible Population
   if (isnan(A[3])) A[3] =  0.4L;                     // sigma - incubation rate  (2.5 d (latency) until an infected individual becomes infectious)
   if (isnan(A[4])) A[4] =  0.0625L;                  // gamma - removal rate (more 16 d until the infectious individual can be removed from the chain of infection)
   if (isnan(A[6])) A[6] =  min;                      // R(t1) boundary value at t1
   if (isnan(A[7])) A[7] =  t1;
   if (isnan(A[8])) A[8] =  0.0L;                     // d·TP  - diffusion constant x Total Population, for example 0.05*80e6 = 4000000

   if (isnan(A[2])) A[2] = A[6]/A[4]/A[3];            // total number of exposed individuals at t1
   if (isnan(A[5])) A[5] = (1.0L - A[4])*A[3]*A[2];   // I(t1) boundary value at t1

   int i, k = 0;
   for (i = 0; i < mpar; i++)
      if (f[i] != undefined) k++;
   if (k == 0)
      f[k++] = 0, f[k++] = 1, f[k++] = 2;

   return k;
}

static void seirdes(ldouble t, ldouble *Y, ldouble *dY, ldouble A[mpar])
{
   dY[0] = -A[0]/A[1]*Y[0]*Y[2] + A[8]/Y[0];          // dS/dt
   dY[1] =  A[0]/A[1]*Y[0]*Y[2] - A[3]*Y[1];          // dE/dt
   dY[2] =  A[3]*Y[1] - A[4]*Y[2];                    // dI/dt
   dY[3] =  A[4]*Y[2];                                // dR/dt
}

int modelFunction_SEIR(ldouble t, ldouble *Y, ldouble A[mpar], int f[mpar], bool fit, bool init)
{
   int rc;

   static ldouble t0;
   static ldouble Y0[4];

   if (init)
   {
      bool isvar2 = false,
           isvar5 = false;

      for (int i = 0; i < mpar && f[i] != undefined; i++)
         if (f[i] == 2)
            isvar2 = true;
         else if (f[i] == 5)
            isvar5 = true;

      if (!isvar2) A[2] = A[6]/A[4]/A[3];
      if (!isvar5) A[5] = (1.0L - A[4])*A[3]*A[2];

      Y0[0] = A[1] - A[2] - A[5] - A[6];
      Y0[1] = A[2];
      Y0[2] = A[5];
      Y0[3] = A[6];
      rc = ODEInt(4, A[7], t, Y0, A, seirdes);
   }
   else
      rc = ODEInt(4, t0, t, Y0, A, seirdes);

   t0 = t;

   if (fit)
      *Y = Y0[3];  // R - curve fit of R(t)
   else
      Y[0] = Y0[0],
      Y[1] = Y0[1],
      Y[2] = Y0[2],
      Y[3] = Y0[3],
      Y[4] = NAN;

   return rc;
}


char *modelDescription_SEIR_de =
"# Model: SEIR Differential Equations\n"\
"# || f = 1, c = 0\n"\
"# || f =     9 if 143 < t and t <= 147   -- Gütersloh/Göttingen\n"\
"# || f =     6 if 172 < t and t <= 177   -- Vechta, Mettmann, ...\n"\
"# || f =     5 if 188 < t and t <= 208   -- ... Mamming, ...\n"\
"# || f =     3 if 208 < t and t <= 251   -- ... Parties/Coronades, ...\n"\
"# || c =   250 if 160 < t and t <= 208   -- Vacation & Outdoor in July 2020\n"
"# || c =   500 if 209 < t and t <= 230   -- Back to school in August 2020\n"
"# || c =  2000 if 230 < t and t <= 251   -- Back to school in September 2020\n"
"# || c = 10500 if 251 < t and t <= 271   -- Begin of a damp and cold season\n"
"# || c = 23000 if 271 < t and t <= 289   -- Manifestation of a damp and cold season\n"
"# || c = 13000 if 289 < t and t <= 307   -- Lock down light\n"
"# || c = 27000 if 307 < t                -- Winter time\n"
"# S  dy0/dt = -f·a0/a1·y0·y2 + a8/y0 + c || y0(a7) = a1-a2-a5-a6\n"\
"# E  dy1/dt =  f·a0/a1·y0·y2 - a3·y1     || y1(a7) = a2 <- a6/a4/a3\n"\
"# I  dy2/dt =  a3·y1 - a4·y2             || y2(a7) = a5 <- (1 - a4)·a3·a2\n"\
"# R  dy3/dt =  a4·y2                     || y3(a7) = a6";

int initialValues_SEIR_de(ldouble t1, ldouble min, ldouble max, ldouble A[mpar], int f[mpar])
{
   if (isnan(A[0])) A[0] =  0.7L;                     // beta  - infection rate
   if (isnan(A[1])) A[1] =  160000.0L;                // VSP   - Virtual Susceptible Population
   if (isnan(A[3])) A[3] =  0.25L;                    // sigma - incubation rate  (4 d (latency) until an infected individual becomes infectious)
   if (isnan(A[4])) A[4] =  0.09L;                    // gamma - removal rate (more 11 d until the infectious individual can be removed from the chain of infection)
   if (isnan(A[6])) A[6] =  min;                      // R(t1) boundary value at t1
   if (isnan(A[7])) A[7] =  t1;
   if (isnan(A[8])) A[8] =  0.0L;                     // d·TP  - diffusion constant x Total Population, for example 0.05*80e6 = 4000000

   if (isnan(A[2])) A[2] = A[6]/A[4]/A[3];            // total number of exposed individuals at t1
   if (isnan(A[5])) A[5] = (1.0L - A[4])*A[3]*A[2];   // I(t1) boundary value at t1

   int i, k = 0;
   for (i = 0; i < mpar; i++)
      if (f[i] != undefined) k++;
   if (k == 0)
      f[k++] = 0, f[k++] = 1, f[k++] = 2;

   return k;
}

static void seirdes_de(ldouble t, ldouble *Y, ldouble *dY, ldouble A[mpar])
{
   ldouble f = 1.0L, c = 0.0L;                        // f is an acceleration factor which may serve to model local outbreaks
                                                      // 1 = no acceleration
   if (143.0L < t && t <= 147.0L)
      f = 9.0L;                                       // Gütersloh/Göttingen
   else if (172.0L < t && t <= 177.0L)
      f = 6.0L;                                       // Vechta, Mettmann, ...
   else if (188.0L < t && t <= 208.0L)
      f = 5.0L;                                       // ... Mamming, ...
   else if (208.0L < t && t <= 251.0L)
      f = 3.0L;                                       // ... Parties/Coronades, ...

   if (160.0L < t && t <= 208.0L)                     // c is a constant summand to the virtual susceptibles and may serve for modeling behavioural and seasonal changes
      c =   250.0L;                                   // Vacation & Outdoor in July 2020
   else if (208.0L < t && t <= 230.0L)
      c =   500.0L;                                   // Back to school in August 2020
   else if (230.0L < t && t <= 251.0L)
      c =  2000.0L;                                   // Back to school in September 2020
   else if (251.0L < t && t <= 271.0L)
      c = 10500.0L;                                   // Begin of a damp and cold season
   else if (271.0L < t && t <= 289.0L)
      c = 23000.0L;                                   // Manifestation of a damp and cold season
   else if (289.0L < t && t <= 307.0L)
      c = 13000.0L;                                   // Lock down light
   else if (307.0L < t)
      c = 27000.0L;                                   // Winter time

   dY[0] = -f*A[0]/A[1]*Y[0]*Y[2] + A[8]/Y[0] + c;    // dS/dt
   dY[1] =  f*A[0]/A[1]*Y[0]*Y[2] - A[3]*Y[1];        // dE/dt
   dY[2] =    A[3]*Y[1] - A[4]*Y[2];                  // dI/dt
   dY[3] =    A[4]*Y[2];                              // dR/dt
}

int modelFunction_SEIR_de(ldouble t, ldouble *Y, ldouble A[mpar], int f[mpar], bool fit, bool init)
{
   int rc;

   static ldouble t0;
   static ldouble Y0[4];

   if (init)
   {
      bool isvar2 = false,
           isvar5 = false;

      for (int i = 0; i < mpar && f[i] != undefined; i++)
         if (f[i] == 2)
            isvar2 = true;
         else if (f[i] == 5)
            isvar5 = true;

      if (!isvar2) A[2] = A[6]/A[4]/A[3];
      if (!isvar5) A[5] = (1.0L - A[4])*A[3]*A[2];

      Y0[0] = A[1] - A[2] - A[5] - A[6];
      Y0[1] = A[2];
      Y0[2] = A[5];
      Y0[3] = A[6];
      rc = ODEInt(4, A[7], t, Y0, A, seirdes_de);
   }
   else
      rc = ODEInt(4, t0, t, Y0, A, seirdes_de);

   t0 = t;

   if (fit)
      *Y = Y0[3];  // R - curve fit of R(t)
   else
      Y[0] = Y0[0],
      Y[1] = Y0[1],
      Y[2] = Y0[2],
      Y[3] = Y0[3],
      Y[4] = NAN;

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
   if (isnan(A[0])) A[0] =  0.6L;                     // alpha  - infection rate
   if (isnan(A[1])) A[1] = 80.0e6L;                   // population
   if (isnan(A[2])) A[2] = 10.0L*min;                 // total number of infectious individuals at t1
   if (isnan(A[3])) A[3] =  0.38L;                    // beta   - removal rate
   if (isnan(A[4])) A[4] =  0.001L;                   // kappa0 - general containment rate (all)
   if (isnan(A[5])) A[5] =  0.0005L;                  // kappa  - quarantine rate (infected)
   if (isnan(A[6])) A[6] = A[3]*A[2] + A[4]*A[1];     // R(t1) boundary value at t1
   if (isnan(A[7])) A[7] = (A[4] + A[5])*A[2];        // X(t1) boundary value at t1
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

int modelFunction_SIRX(ldouble t, ldouble *Y, ldouble A[mpar], int f[mpar], bool fit, bool init)
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
   if (fit)
      *Y = Y0[3];  // X - curve fit of X(t)
   else
      Y[0] = Y0[0],
      Y[1] = Y0[1],
      Y[2] = Y0[2],
      Y[3] = Y0[3],
      Y[4] = NAN;

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

int modelFunction_ERF(ldouble t, ldouble *Y, ldouble A[mpar], int f[mpar], bool fit, bool init)
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

int modelFunction_GLF(ldouble t, ldouble *Y, ldouble A[mpar], int f[mpar], bool fit, bool init)
{
   *Y = A[1]/powl(1 + expl(-A[0]*A[3]*(t - A[2])), 1/A[3]);
   return 0;
}
