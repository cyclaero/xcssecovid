//  numerics.c
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


#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

#include "numerics.h"


#pragma mark ••• Bulirsch–Stoer Differential Equation Solver •••

#define tol (ldexp(__LDBL_EPSILON__, 21))

#define nuse 7
#define imax 13
static int nseq[imax] = { 2, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 256 };

void mmidpt(int m, int steps, ldouble t0, ldouble htot, ldouble *Y0, ldouble *dY, ldouble *Y, ldouble *A, odeset model)
{
   int      i, j;
   ldouble  h, h2, t;
   ldouble  Ym[m], Yn[m], Ys[m];
   ldouble *z, *Zm = Ym, *Zn = Yn, *Zs = Ys;

   memcpy(Ym, Y0, m*sizeof(ldouble));
   for (h = htot/steps, j = 0; j < m; j++)
      Yn[j] = Y0[j] + h*dY[j];

   model(t = t0 + h, Yn, Y, A);
   for (h2 = 2.0L*h, i = 1; i < steps; i++)
   {
      for (j = 0; j < m; j++)
         Zs[j] = Zm[j] + h2*Y[j];
      z = Zm; Zm = Zn; Zn = Zs; Zs = z;
      model(t += h, Zn, Y, A);
   }

   for (j = 0; j < m; j++)
      Y[j] = 0.5L*(Zm[j] + Zn[j] + h*Y[j]);
}

void rzextr(int i, int m, ldouble *Yest, ldouble *Y, ldouble *dY, ldouble rzt[imax], ldouble rzdY[nuse][m])
{
   int    j, k, m1 = (i < nuse-1) ? i : nuse-1;
   ldouble v, yj, dyj = 0.0L;
   ldouble b, b1, c;
   ldouble fx[nuse];

   for (k = 0; k < m1; k++)
      fx[k+1] = rzt[i-k-1]/rzt[i];

   for (j = 0; j < m; j++)
   {
      v = rzdY[0][j];
      c = yj = rzdY[0][j] = Yest[j];
      for (k = 1; k <= m1; k++)
      {
         b1 = v*fx[k];
         b = b1 - c;
         if (fabsl(b) > __LDBL_EPSILON__)
         {
            b = (c - v)/b;
            dyj = c*b;
            c = b1*b;
         }
         else
            dyj = v;

         if (k != m1)
            v = rzdY[k][j];

         yj += dyj;
         rzdY[k][j] = dyj;
      }

       Y[j] = yj;
      dY[j] = dyj;
   }
}

deq_error bsstep(int m, ldouble h, ldouble *t, ldouble *hdid, ldouble *hnext, ldouble *Y, ldouble *dY, ldouble *fY, ldouble *A, odeset model)
{
   int      i, j;
   ldouble  err, errmax;
   ldouble  Y0[m], Yseq[m], Yerr[m];

   ldouble  rzt[imax];
   ldouble rzdY[nuse][m];

   memcpy(Y0, Y, m*sizeof(ldouble));
   for (;;)
   {
      for (i = 0; i < imax; i++)
      {
         mmidpt(m, nseq[i], *t, h, Y0, dY, Yseq, A, model);
         rzt[i] = sqrl(h/nseq[i]);
         if (i == 0)
         {
            memcpy(Y,       Yseq, m*sizeof(ldouble));
            memcpy(Yerr,    Yseq, m*sizeof(ldouble));
            memcpy(rzdY[0], Yseq, m*sizeof(ldouble));
         }
         else
            rzextr(i, m, Yseq, Y, Yerr, rzt, rzdY);

         if (i > 2)
         {
            for (errmax = 0.0L, j = 0; j < m; j++)
               if (!isfinite(Yerr[j]))
                  goto tolErrCheck;
               else if ((err = fabsl(Yerr[j]/fY[j])) > errmax)
                  errmax = err;

            if (errmax < tol)
            {
               *t += h;
               *hdid = h;
               if (i == nuse-1)
                  *hnext = 0.95L*h;
               else if (i == nuse-2)
                  *hnext = 1.2L*h;
               else
                  *hnext = h*nseq[nuse-2]/nseq[i];

               return deq_noError;
            }
         }
      }

   tolErrCheck:
      for (h *= 0.25L, i = 0; i < (imax-nuse)/2; i++)
         h *= 0.5L;
      if (*t + h == *t)
         return deq_tolError;
   }
}

deq_error ODEInt(int m, ldouble t1, ldouble t2, ldouble *Y, ldouble *A, odeset model)
{
   if (t1 == t2)
      return deq_noError;

   int     i, j, err;
   ldouble t = t1;
   ldouble h = t2 - t1;
   ldouble hnext, hdid;
   ldouble dY[m], fY[m];

   for (i = 0; i < 10000; i++)
   {
      model(t, Y, dY, A);
      for (j = 0; j < m; j++)
         fY[j] = fabsl(Y[j]) + fabsl(h*dY[j]) + __LDBL_EPSILON__;

      if ((t + h - t2)*(t + h - t1) > __LDBL_EPSILON__)
         h = t2 - t;

      if (err = bsstep(m, h, &t, &hdid, &hnext, Y, dY, fY, A, model))
         return err;

      if ((t - t2)*(t2 - t1) >= 0.0L)
         return deq_noError;

      if (fabsl(hnext) < __LDBL_EPSILON__)
         return deq_stiffError;

      h = hnext;
   }

   return deq_stepsError;
}


#pragma mark ••• Curve Fitting by Levenberg–Marquardt least squares minimization •••
// https://en.wikipedia.org/wiki/Levenberg–Marquardt_algorithm

ldouble calcChiSqr(int n, ldouble *T, ldouble *Y,
                   int k, ldouble *A, function curve)
{
   ldouble yi, chiSqr = 0.0L;

   for (int i = 0; i < n; i++)
   {
      curve(T[i], &yi, A, i == 0);
      chiSqr += sqrl(yi - Y[i]);
   }

   if (--n - k > 0)
      return chiSqr/(n - k);
   else if (n > 0)
      return chiSqr/n;
   else
      return NAN;
}

ldouble calcGradientCurvature(int n, ldouble *T, ldouble *Y,
                              int k, ldouble *A, ldouble *dA, int *f,
                              ldouble *beta, ldouble **alpha,
                              ldouble lambda, function curve)
{
   int       i, p, q;
   ldouble   ap, app, dYdA, yi;
   ldouble **V = malloc(k*sizeof(ldouble*));
   for (p = 0; p < k; p++)
      V[p] = calloc(n, sizeof(ldouble));

   for (p = 0; p < k; p++)
   {
      beta[p] = 0.0L;
      for (q = 0; q < k; q++)
         alpha[p][q] = 0.0L;

      ap = A[f[p]], A[f[p]] += dA[p];

      for (i = 0; i < n; i++)
         curve(T[i], &V[p][i], A, i == 0);

      A[f[p]] = ap;
   }

   for (i = 0; i < n; i++)
   {
      curve(T[i], &yi, A, i == 0);
      for (p = 0; p < k; p++)
      {
         dYdA = (V[p][i] - yi)/dA[p];
         beta[p] += dYdA*(yi - Y[i]);
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

   if (fabsl(lambda) < __LDBL_EPSILON__)
   {
      for (p = 0; p < k; p++)
         if ((app = fabsl(alpha[p][p])) > lambda)
            lambda = app;
      lambda = sqrtl(lambda)/n;
   }

   for (p = 0; p < k; p++)
   {
      alpha[p][p] *= (1.0L + lambda);
      free(V[p]);
   }
   free(V);

   return lambda;
}

ldouble ERelNorm(int m, ldouble *B, ldouble *R)
{
   ldouble sqSum = 0.0L;

   for (int i = 0; i < m; i++)
      sqSum += (fabsl(R[i]) > __LDBL_EPSILON__ && isfinite(R[i]))
             ? sqrl(B[i]/R[i])
             : sqrl(B[i]);

   return sqrtl(sqSum);
}

ldouble curveFit(int n, ldouble *T, ldouble *Y,
                 int k, ldouble *A, ldouble *dA, int *f, function curve)
{
   #define rtol 1.0e-9L
   #define atol 1.0e-12L

   static ldouble unitv[mpar] = {1.0L, 1.0L, 1.0L, 1.0L, 1.0L, 1.0L, 1.0L, 1.0L, 1.0L, 1.0L};

   bool      revoke;
   int       i, j, iter = 0;
   int       badIter = 0;
   int       idx[k];
   ldouble   chiSqr, chiSqr0, dChiSqr, dRelB = NAN, lambda0 = 0.0L, lambda = 0.0L;
   ldouble   beta[k], delta[k], B[k], dB[k], C[k];
   ldouble **alpha, **alphaLU;

   alpha   = malloc(k*sizeof(ldouble *));
   alphaLU = malloc(k*sizeof(ldouble *));
   for (j = 0; j < k; j++)
   {
       B[j]      = A[f[j]];
      dB[j]      = (fabsl(B[j]) > __LDBL_EPSILON__) ? fabsl(B[j])*rtol : rtol;
      alpha[j]   = calloc(k, sizeof(ldouble));
      alphaLU[j] = calloc(k, sizeof(ldouble));
   }

   chiSqr = calcChiSqr(n, T, Y, k, A, curve);
   do
   {
      iter++;
      lambda = calcGradientCurvature(n, T, Y,
                                     k, A, dB, f,
                                     beta, alpha, lambda, curve);
      if (!isfinite(dRelB = LUdecomposition(k, alpha, alphaLU, idx)))
         break;
      LUbacksubstitution(k, alphaLU, idx, unitv, delta);
      LUrefinment(k, alpha, alphaLU, idx, unitv, delta);

      dRelB = ERelNorm(k, delta, B);

      for (j = 0; j < k; j++)
         C[j] = B[j], B[j] -= delta[j], A[f[j]] = B[j];

      chiSqr0 = chiSqr;
      chiSqr  = calcChiSqr(n, T, Y, k, A, curve);
      dChiSqr = chiSqr - chiSqr0;
      revoke  = dChiSqr > 0.0L;

      lambda0 = lambda;
      if (!revoke)
         lambda *= 0.2L;
      else
      {
         lambda = (lambda < 100.0L) ? lambda*3.0L : 0.0L;
         for (j = 0; j < k; j++)
            A[f[j]] = B[j] = C[j];
         chiSqr = chiSqr0;

         if (dChiSqr < atol)
            badIter++;
         else
            badIter = 0;
      }

      for (j = 0; j < k; j++)
         A[f[j]] = B[j];
   }
   while (iter < 1000 && badIter < 10 && isfinite(dRelB) && (dRelB > atol && fabsl(dChiSqr) > atol || dChiSqr > 0));

   if (isfinite(dRelB))
   {
      // extract lambda for preparation of the covarince matrix
      for (i = 0; i < k; i++)
      {
         for (j = 0; j < k; j++)
            alpha[i][j] *= beta[i];
         alpha[i][i] /= (1.0L + lambda0);
      }

      ldouble det = LUdecomposition(k, alpha, alphaLU, idx);
      if (fabsl(det) > __LDBL_EPSILON__ && isfinite(det))
      {
         LUinversion(k, alpha, alphaLU, idx);
         for (j = 0; j < k; j++)
            dA[f[j]] = sqrtl(chiSqr*fabsl(alpha[j][j]))*100.0L/fabsl(B[j]);
      }
   }
   else
      chiSqr = NAN;

   for (j = 0; j < k; j++)
   {
      free(alpha[j]);
      free(alphaLU[j]);
   }
   free(alpha);
   free(alphaLU);

   #undef rtol
   #undef atol

   return chiSqr;
}


#pragma mark ••• LU Decomposition •••
// https://en.wikipedia.org/wiki/LU_decomposition

ldouble LUdecomposition(int m, ldouble **A, ldouble **LU, int *idx)
{
   int      i, j, k, maxi = 0;
   ldouble  max, sum, dum, d;
   ldouble *V = calloc(m, sizeof(ldouble));

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

      if (fabsl(max) > __LDBL_EPSILON__)
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
            maxi = i;
         }
      }

      if (j != maxi)
      {
         for (k = 0; k < m; k++)
         {
            dum = LU[maxi][k];
            LU[maxi][k] = LU[j][k];
            LU[j][k] = dum;
         }
         V[maxi] = V[j];
         d = -d;
      }
      idx[j] = maxi;

      if (fabsl(LU[j][j]) < __LDBL_EPSILON__)
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

void LUbacksubstitution(int m, ldouble **LU, int *idx, ldouble *B, ldouble *X)
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
      else if (fabsl(sum) > __LDBL_EPSILON__)
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

void LUrefinment(int m, ldouble **A, ldouble **LU, int *idx, ldouble *B, ldouble *X)
{
   int      i, j;
   ldouble  sdp;
   ldouble *R = calloc(m, sizeof(ldouble));

   for (i = 0; i < m; i++)
   {
      sdp = -B[i];
      for (j = 0; j < m; j++)
         sdp += A[i][j]*X[j];
      R[i] = sdp;
   }
   LUbacksubstitution(m, LU, idx, R, R);

   for (i = 0; i < m; i++)
      X[i] -= R[i];

   free(R);
}

void LUinversion(int m, ldouble **A, ldouble **LU, int *idx)
{
   int       i, j;
   ldouble **Ai = malloc(m*sizeof(ldouble*));
   for (j = 0; j < m; j++)
      Ai[j] = calloc(m, sizeof(ldouble));

   for (i = 0; i < m; i++)
      Ai[i][i] = 1.0L;

   for (j = 0; j < m; j++)
      LUbacksubstitution(m, LU, idx, Ai[j], Ai[j]);

   for (i = 0; i < m; i++)
      for (j = 0; j < m; j++)
         A[i][j] = Ai[i][j];

   for (j = 0; j < m; j++)
      free(Ai[j]);
   free(Ai);
}
