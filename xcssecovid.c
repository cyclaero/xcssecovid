//  xcssecovid.c
//  xcssecovid
//
//  Created by Dr. Rolf Jansen on 2020-03-22.
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
//
//  Usage:
//
//  1. Compile this file on either of FreeBSD, Linux or macOS:
//
//     cc -g0 -O3 -march=native xcssecovid.c -Wno-parentheses -lm -o xcssecovid
//
//  2. Download the daily updated time series of confirmed Covid-19 cases
//     from CSSE's (at Johns Hopkins University) GitHub site CSSE COVID-19 Dataset - https://github.com/CSSEGISandData/COVID-19/tree/master/csse_covid_19_data
//
//     curl -O https://github.com/CSSEGISandData/COVID-19/raw/master/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Confirmed.csv
//
//  3. Extract the time series row of a given country starting on 2020-01-21 as day 0
//     and write it out together with the country's cases and a simulated curve-fitted
//     logistic function into a three-column TSV output file:
//
//     ./xcssecovid Germany time_series_19-covid-Confirmed.csv Germany.tsv
//
//  4. Open the TSV file with your favorite graphing and/or data analysis application.


#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>


typedef long double ldouble;
typedef ldouble    *Vector;
typedef ldouble   **Matrix;
typedef int        *Index;


static ldouble sqr(ldouble x)
{
   return x*x;
}

ldouble LogisticFunctionGen(ldouble t, Vector A)
{  // generic form: https://en.wikipedia.org/wiki/Logistic_function - parameters A[3] to A[5] = 0
   return A[0]/(1 + exp(-A[1]*(t - A[2]))) - (A[3] + A[4]*t + A[5]*sqr(t));
}


ldouble LogisticFunctionAlt(ldouble t, Vector A)
{  // alternate form: https://de.wikipedia.org/wiki/Logistische_Funktion#Weitere_Darstellungen - parameters A[3] to A[5] = 0
   return A[0]*0.5/(0.5 + (A[0] - 0.5)*exp(-A[1]*A[0]*(t - A[2]))) - (A[3] + A[4]*t + A[5]*sqr(t));
}

ldouble curveFit(int n, Vector T, Vector C, Vector S,
                 int k, Vector A, Vector dA,
                 ldouble (*model)(ldouble t, Vector A));

static inline size_t collen(const char *col)
{
   if (!col || !*col)
      return 0;

   size_t l;
   for (l = 0; col[l] && col[l] != ','; l++)
      ;
   return l;
}

int main(int argc, char *const argv[])
{
   FILE *csv, *tsv;
   char *country = argv[1];

   if (*country)
      if (csv = (*(uint16_t *)argv[2] == *(uint16_t *)"-")
                ? stdin
                : fopen(argv[2], "r"))
      {
         if (tsv = (*(uint16_t *)argv[3] == *(uint16_t *)"-")
                   ? stdout
                   : fopen(argv[3], "w"))
         {
            int     i;
            char   *line;
            char    d[65536];
            ldouble t[366];  // 2020 is a leap year - t[0] = 2020-01-01; t[365] = 2020-12-31
            ldouble c[366];
            ldouble l[366];

            // day 1 of the CSSE time series is 2020-01-22, hence day 0 is 2020-01-21, i.e. t[20]
            for (i = 0; i < 366; i++)
               t[i] = i - 20, c[i] = l[i] = NAN;

            // find the line with the series of the specified country
            while (line = fgets(d, 65536, csv))
            {
               size_t len = strlen(line);
               size_t col = collen(line);
               if (len != col && (line = strcasestr(line+col+1, country)))
                  break;   // found!
            }


            if (line)
            {
               // skip the 3 more fields (Country/Region,Lat,Long)
               for (i = 0; i < 3; i++)
                  line += collen(line) + 1;

               // read the case numbers into the n array, starting at day 1 = index 20+1;
               int p, q = p = 21;
               char *chk;
               while (*line && *line != '\n' && *line != '\r')
               {
                  ldouble v = strtod(line, &chk);
                  if (chk > line)
                     c[q] = v, line = chk + 1;
                  else
                     line += collen(line) + 1;
                  q++;
               }

               if (q > p)
               {
                  ldouble (*model)(ldouble t, Vector A) = LogisticFunctionGen;

                  int j, k;
                  ldouble ChiSqr,
                          A[6] = {2.0L*c[q-1], (model == LogisticFunctionGen) ? 0.2L : 5.0e-6L, 1.0L, 0.0L, 0.0L, 0.0L},
                         dA[6] = {};

                  char *funcStr = (model == LogisticFunctionGen)
                                ? "a/(1 + exp(-b·(x - c)))"
                                : "a·0.5/(0.5 + (a - 0.5)·exp(-b·a·(x - c)))";

                  if (isfinite(ChiSqr = curveFit(q - p, &t[p], &c[p], &l[p], k = 1, A, dA, model)))
                  {
                     dA[0] = 0.0;
                     if (isfinite(ChiSqr = curveFit(q - p, &t[p], &c[p], &l[p], k = 3, A, dA, model)))
                     {
                        if (k <= 3)
                           fprintf(tsv, "# Model: %s\n", funcStr);
                        else if (k == 5)
                           fprintf(tsv, "# Model: %s\n", funcStr);
                        else if (k == 6)
                           fprintf(tsv, "# Model: %s\n", funcStr);

                        for (j = 0; j < k; j++)
                           fprintf(tsv, "#        %c = %9.6Lg ± %.5Lg %%\n", 'a'+j, A[j], dA[j]);
                        fprintf(tsv, "#   ChiSqr = %9.7Lg\n", ChiSqr);
                     }
                     else
                        fprintf(tsv, "# Curve fit failed\n");
                  }

                  // write the column header formular symbols and units.
                  // - the formular symbol of time is 't', the unit symbol of day is 'd'
                  // - the formular symbol of number of cases is C without a unit
                  // - the formular symbol of the siumulated loogistic function is L without a unit
                  fprintf(tsv, "t/d\tC\tL\n");
                  for (i = 0; i < 366; i++)
                     if (isfinite(c[i]))
                        fprintf(tsv, "%.0Lf\t%.0Lf\t%.6Lf\n", t[i], c[i], l[i] = model(t[i], A));
                     else
                        fprintf(tsv, "%.0Lf\t*\t%.6Lf\n", t[i], l[i] = model(t[i], A));
               }
               else
                  fprintf(tsv, "# No values for country %s encountered.\n", country);
            }
            else
               fprintf(tsv, "# No data for country %s found.\n", country);

            if (tsv != stdout)
               fclose(tsv);
         }

         if (csv != stdin)
            fclose(csv);
      }

   return 0;
}


// LU decomposition
ldouble LUdcmp(int m, Matrix A, Matrix LU, Index idx)
{
   int    i, j, k, imax = 0;
   ldouble max, sum, dum, d;
   Vector V = calloc(m, sizeof(ldouble));

   if (LU != A)
      for (i = 0; i < m; i++)
         for (j = 0; j < m; j++)
            LU[i][j] = A[i][j];

   for (i = 0; i < m; i++)
   {
      max = 0.0L;
      for (j = 0; j < m; j++)
         if ((dum = fabsl(LU[i][j])) > max)
            max = dum;

      if (max != 0.0L)
         V[i] = 1.0L/max;
      else
      {
         free(V);
         return NAN;
      }
   }

   d = 1.0L;
   for (j = 0; j < m; j++)
   {
      for (i = 0; i < j; i++)
      {
         sum = LU[i][j];
         for (k = 0; k < i; k++)
            sum -= LU[i][k]*LU[k][j];
         LU[i][j] = sum;
      }

      max = 0.0L;
      for (; i < m; i++)
      {
         sum = LU[i][j];
         for (k = 0; k < j; k++)
            sum -= LU[i][k]*LU[k][j];
         LU[i][j] = sum;

         if ((dum = V[i]*fabsl(sum)) >= max)
         {
            max = dum;
            imax = i;
         }
      }

      if (j != imax)
      {
         for (k = 0; k < m; k++)
         {
            dum = LU[imax][k];
            LU[imax][k] = LU[j][k];
            LU[j][k] = dum;
         }
         V[imax] = V[j];
         d = -d;
      }
      idx[j] = imax;

      if (LU[j][j] == 0.0L)
         LU[j][j] = __LDBL_EPSILON__;

      if (j < m-1)
      {
         dum = 1.0L/LU[j][j];
         for (i = j+1; i < m; i++)
            LU[i][j] *= dum;
      }
   }

   free(V);
   return d;
}


void LUbksb(int m, Matrix LU, Index idx, Vector B, Vector X)
{
   int     i, j, k, l = -1;
   ldouble sum;

   if (X != B)
      for (i = 0; i < m; i++)
         X[i] = B[i];

   for (i = 0; i < m; i++)
   {
      k = idx[i];
      sum = X[k];
      X[k] = X[i];
      if (l >= 0)
         for (j = l; j <= i-1; j++)
            sum -= LU[i][j]*X[j];
      else if (sum != 0.0L)
         l = i;
      X[i] = sum;
   }

   for (i = m-1; i >= 0; i--)
   {
      sum = X[i];
      for (j = i+1; j < m; j++)
         sum -= LU[i][j]*X[j];
      X[i] = sum/LU[i][i];
   }
}


void LUimpr(int m, Matrix A, Matrix LU, Index idx, Vector B, Vector X)
{
   int    i, j;
   ldouble sdp;
   Vector R = calloc(m, sizeof(ldouble));

   for (i = 0; i < m; i++)
   {
      sdp = -B[i];
      for (j = 0; j < m; j++)
         sdp += A[i][j]*X[j];
      R[i] = sdp;
   }
   LUbksb(m, LU, idx, R, R);

   for (i = 0; i < m; i++)
      X[i] -= R[i];

   free(R);
}


void LUInvr(int m, Matrix A, Matrix LU, Index idx)
{
   int  i, j;
   Matrix  Ai = malloc(m*sizeof(ldouble*));
   for (j = 0; j < m; j++)
      Ai[j] = calloc(m, sizeof(ldouble));

   for (i = 0; i < m; i++)
      Ai[i][i] = 1.0L;

   for (j = 0; j < m; j++)
      LUbksb(m, LU, idx, Ai[j], Ai[j]);

   for (i = 0; i < m; i++)
      for (j = 0; j < m; j++)
         A[i][j] = Ai[i][j];

   for (j = 0; j < m; j++)
      free(Ai[j]);
   free(Ai);
}


ldouble ERelNorm(int m, Vector B, Vector R)
{
   int    i;
   ldouble sqSum = 0.0L;

   for (i = 0; i < m; i++)
      if (R[i] != 0.0L && isfinite(R[i]))
         sqSum += sqr(B[i]/R[i]);
      else
         sqSum += sqr(B[i]);
   return sqrt(sqSum);
}


ldouble calcChiSqr(int n, Vector T, Vector C, Vector S,
                   int k, Vector A,
                   ldouble (*model)(ldouble t, Vector A))
{
   int     i, cnt = -1;
   ldouble chiSqr = 0.0L;

   for (i = 0; i < n; i++)
   {
      S[i] = model(T[i], A);
      chiSqr += sqr(S[i] - C[i]);
      cnt++;
   }

   if (cnt - k > 0)
      return chiSqr/(cnt - k);
   else if (cnt > 0)
      return chiSqr/cnt;
   else
      return NAN;
}


ldouble calcGradientCurvature(int n, Vector T, Vector  C,
                              int k, Vector A, Vector dA,
                              Vector beta, Matrix alpha, ldouble lambda,
                              ldouble (*model)(ldouble t, Vector A))
{
   int      i, p, q, cnt = 0;
   ldouble  ap, app, dYdA, yi;
   ldouble **V = malloc(k*sizeof(ldouble*));
   for (p = 0; p < k; p++)
      V[p] = calloc(n, sizeof(ldouble));

   for (p = 0; p < k; p++)
   {
      beta[p] = 0.0L;
      for (q = 0; q < k; q++)
         alpha[p][q] = 0.0L;

      ap = A[p];
      A[p] += dA[p];

      for (i = 0; i < n; i++)
         V[p][i] = model(T[i], A);

      A[p] = ap;
   }

   for (i = 0; i < n; i++)
   {
      cnt++;
      yi = model(T[i], A);
      for (p = 0; p < k; p++)
      {
         dYdA = (V[p][i] - yi)/dA[p];
         beta[p] += dYdA*(yi - C[i]);
         for (q = 0; q <= p; q++)
            alpha[p][q] += dYdA*(V[q][i] - yi)/dA[q];
      }
   }

   for (p = 1; p < k; p++)
      for (q = 0; q < p; q++)
         alpha[q][p] = alpha[p][q];

   for (p = 0; p < k; p++)
      for (q = 0; q < k; q++)
         alpha[p][q] /= beta[p];

   if (lambda == 0.0L)
   {
      for (p = 0; p < k; p++)
         if ((app = fabsl(alpha[p][p])) > lambda)
            lambda = app;
      lambda = sqrt(lambda)/cnt;
   }

   for (p = 0; p < k; p++)
   {
      alpha[p][p] *= (1.0L + lambda);
      free(V[p]);
   }
   free(V);
   return lambda;
}



ldouble curveFit(int n, Vector T, Vector C, Vector S,
                int k, Vector A, Vector dA,
                ldouble (*model)(ldouble t, Vector A))
{
   #define rtol 1.0e-9L
   #define atol 1.0e-12L

   static ldouble unitvector[10] = {1.0L, 1.0L, 1.0L, 1.0L, 1.0L, 1.0L, 1.0L, 1.0L, 1.0L, 1.0L};

   bool      revoke;
   int       i, j, iter = 0;
   int       badIter = 0;
   int       idx[k];
   ldouble   ChiSqr, ChiSqr0, dChiSqr, dRelA = NAN, lambda0 = 0.0L, lambda = 0.0L;
   ldouble   beta[k], delta[k], B[k];
   ldouble **alpha, **alphaLU;

   alpha   = malloc(k*sizeof(ldouble *));
   alphaLU = malloc(k*sizeof(ldouble *));
   for (j = 0; j < k; j++)
   {
      dA[j]      = (A[j] != 0.0L) ? fabsl(A[j])*rtol : rtol;
      alpha[j]   = calloc(k, sizeof(ldouble));
      alphaLU[j] = calloc(k, sizeof(ldouble));
   }

   ChiSqr = calcChiSqr(n, T, C, S, k, A, model);
   do
   {
      iter++;
      lambda = calcGradientCurvature(n, T, C,
                                     k, A, dA, beta, alpha,
                                     lambda, model);
      if (!isfinite(LUdcmp(k, alpha, alphaLU, idx)))
         break;
      LUbksb(k, alphaLU, idx, unitvector, delta);
      LUimpr(k, alpha, alphaLU, idx, unitvector, delta);

      dRelA = ERelNorm(k, delta, A);

      for (j = 0; j < k; j++)
         B[j] = A[j], A[j] -= delta[j];

      ChiSqr0 = ChiSqr;
      ChiSqr  = calcChiSqr(n, T, C, S, k, A, model);
      dChiSqr = ChiSqr - ChiSqr0;
      revoke  = dChiSqr >= 0;

      lambda0 = lambda;
      if (!revoke)
         lambda *= 0.2L;
      else
      {
         lambda = (lambda < 100.0L) ? lambda*3.0L : 0.0L;
         for (j = 0; j < k; j++)
            A[j] = B[j];
         ChiSqr = ChiSqr0;

         if (dChiSqr < atol)
            badIter++;
         else
            badIter = 0;
      }
   }
   while (iter < 1000 && badIter < 10 && isfinite(dRelA) && (dRelA > atol && fabsl(dChiSqr) > atol || dChiSqr > 0));

   if (isfinite(dRelA))
   {
      // extract lambda for preparation of the covarince matrix
      for (i = 0; i < k; i++)
      {
         for (j = 0; j < k; j++)
            alpha[i][j] *= beta[i];
         alpha[i][i] /= (1.0 + lambda0);
      }

      ldouble det = LUdcmp(k, alpha, alphaLU, idx);
      if (det != 0.0L && isfinite(det))
      {
         LUInvr(k, alpha, alphaLU, idx);
         for (j = 0; j < k; j++)
            dA[j] = sqrtl(ChiSqr*fabsl(alpha[j][j]))*100.0L/fabsl(A[j]);
      }
   }
   else
      ChiSqr = NAN;

   for (j = 0; j < k; j++)
   {
      free(alpha[j]);
      free(alphaLU[j]);
   }
   free(alpha);
   free(alphaLU);

   #undef rtol
   #undef atol

   return ChiSqr;
}
