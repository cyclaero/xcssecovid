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
//  1. Compile numerics.c, models.c and xcssecovid.c on either of FreeBSD or macOS:
//
//     clang -g0 -O3 -march=native numerics.c models.c xcssecovid.c -Wno-parentheses -lm -o xcssecovid
//
//  2. Download the daily updated time series of confirmed Covid-19 cases
//     from CSSE's (at Johns Hopkins University) GitHub site CSSE COVID-19 Dataset
//     https://github.com/CSSEGISandData/COVID-19/tree/master/csse_covid_19_data
//
//     fetch https://github.com/CSSEGISandData/COVID-19/raw/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv
//
//  3. Extract the time series row of a given country starting on 2020-01-21 as day 0
//     and write it out together with the country's cases and a simulated curve-fitted
//     logistic function into a three-column TSV output file:
//
//     ./xcssecovid Germany time_series_covid19_confirmed_global.csv Germany.tsv
//
//  4. Open the TSV file with your favorite graphing and/or data analysis application.


#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/stat.h>

#include "numerics.h"
#include "models.h"


static const char *version = "Version 1.0.2";

static inline const char *execname(const char *cmd)
{
   const char *r = cmd + strlen(cmd);
   while (--r >= cmd && *r != '/');
   return r + 1;
}

int usage(const char *exe)
{
   printf("\nExtract, curve fit an epidemiological model and transpose CSSE@JHU's Covid-19 cases data per country - Copyright Dr. Rolf Jansen (c) 2020 - %s\n", version);
   printf("Usage: %s [-a<0-9> value] [-f (0-9)+] [-m model] [-e] [-r] [-s] [-o day#] [-z day#] [-h|-?|?] <Country> <CSV Input file> <TSV Output File>\n\n\
       -a<0..9> value      Optionally set initial values for the model's parameters the Differential Equation Solver and Curve Fitting.\n\
                           The models deduce its initial parameters from the boundaries of the imported time series and by common\n\
                           knowledge/best educated guesses. Example: -a0 0.57 -a1 125000\n\n\
       -f (0-9)+           Overrides the default selection of a model's parameters which take part in curve fitting, and which usually is a0, a1, a2.\n\
                           Example: -f 1245 would lead to curve fitting against the paramters a1, a2, a4 and a5, while a0 and a would be left untouched.\n\
                           Different models got different number of parameters, which is currently 3 to 6.\n\n\
       -m model            Select the model for curve fitting and simulation:\n\
                           - LF   Logistic Function -- https://en.wikipedia.org/wiki/Logistic_function\n\
                           - LDE  Logistic Differential Equation -- https://en.wikipedia.org/wiki/Logistic_function#Logistic_differential_equation\n\
                           - SI   Epidemiological SI-Model (basically another form of the LDE) -- https://de.wikipedia.org/wiki/SI-Modell\n\
                           - SIR  Epidemiological SIR-Model [default] -- https://en.wikipedia.org/wiki/Mathematical_modelling_of_infectious_disease#The_SIR_model\n\
                           - SIRX Extension of the SIR-Model by the HU/RKI -- http://rocs.hu-berlin.de/corona/docs/forecast/model/\n\
                           - ERF  Shifted Error Function -- (generic form) https://en.wikipedia.org/wiki/Error_function\n\
                           - GLF  Generalised Logistic Function -- https://en.wikipedia.org/wiki/Generalised_logistic_function\n\n\
       -e                  Only export the extracted and transposed time series without curve fitting and simulation of the model.\n\n\
       -r                  Only report the model description and the values of its parameters with or without fitting, but without exporting any curve data.\n\n\
       -s                  Simulate the model without curve fitting before.\n\n\
       -t thresh           Threshold of minimal number of cases to include into curve fitting and simulation [default: 17].\n\n\
       -o day#             The day# of the first data point in the imported time series to be included for curve fitting.\n\
                           [default: first day with more than <thresh> cases].\n\n\
       -z day#             The day# of the last data point in the imported time series to be included for curve fitting.\n\
                           [default: last day of the imported series].\n\n\
       -h|-?|?             Show these usage instructions.\n\n\
       <Country>           Select the country for which the time series shall be processed.\n\n\
       <CSV Input file>    Path to the CSSE@JHU's Covid-19 time series CSV file.\n\n\
       <TSV Output File>   Path to the TSV output file containing the extracted and transposed time series for the given <Country>, including\n\
                           a column for a simulated time series by the given model, using the parameter as resulted from curve fitting.\n\n", exe);
   return 1;
}


static inline size_t collen(const char *col)
{
   if (!col || !*col)
      return 0;

   size_t l;
   for (l = 0; col[l] && col[l] != ','; l++)
      ;
   return l;
}

#define ndays 366    // 2020 is a leap year - t[0] = 2020-01-01; t[365] = 2020-12-31

int main(int argc, char *const argv[])
{
   const char *exe = execname(argv[0]);

   char *modelDescription = modelDescription_SIR;
   initvals initialValues = initialValues_SIR;
   function modelFunction = modelFunction_SIR;

   bool  do_simulation = true,
         do_curve_fit  = true,
         export_series = true;

   ldouble thr = 17.0L,
             o = NAN,
             z = NAN;

   ldouble A[mpar] = {NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN},
          dA[mpar] = {};

   int     f[mpar] = {undefined, undefined, undefined, undefined, undefined,
                      undefined, undefined, undefined, undefined, undefined};

   bool aopt = false;
   int  opc;
   while ((opc = getopt(argc, argv, "af:0:1:2:3:4:5:6:7:8:9:m:erst:o:z:h?")) != -1)
   {
      char   *chk, c;
      int     i;
      ldouble v;

      switch (opc)
      {
         case 'a':
            if (aopt)
               return usage(exe);
            aopt = true;
            break;

         case 'f':
            if (optarg && *optarg)
            {
               for (i = 0; c = *optarg; optarg++)
                  if ('0' <= c && c <= '9' && f[i] == -1)
                     f[i++] = c - '0';
                  else
                     return usage(exe);
            }
            else
               return usage(exe);
            break;

         case '0' ... '9':
            if (aopt && optarg && *optarg)
            {
               if ('0' <= opc && opc <= '9' && isnan(A[opc -= '0'])
                && ((v = strtold(optarg, &chk)) != 0.0L || chk != optarg)
                && isfinite(v))
               {
                  A[opc] = v;
                  aopt = false;
               }
               else
                  return usage(exe);
            }
            else
               return usage(exe);
            break;

         case 'm':
            if (optarg && *optarg)
               if (*(uint16_t *)optarg == *(uint16_t *)"LF" && optarg[2] == '\0')
               {
                  modelDescription = modelDescription_LF;
                  initialValues = initialValues_LF;
                  modelFunction = modelFunction_LF;
               }

               else if (*(uint32_t *)optarg == *(uint32_t *)"LDE")
               {
                  modelDescription = modelDescription_LDE;
                  initialValues = initialValues_LDE;
                  modelFunction = modelFunction_LDE;
               }

               else if (*(uint16_t *)optarg == *(uint16_t *)"SI" && optarg[2] == '\0')
               {
                  modelDescription = modelDescription_SI;
                  initialValues = initialValues_SI;
                  modelFunction = modelFunction_SI;
               }

               else if (*(uint32_t *)optarg == *(uint32_t *)"SIR")
                  ; // default, do nothing

               else if (*(uint32_t *)optarg == *(uint32_t *)"SIRX" && optarg[4] == '\0')
               {
                  modelDescription = modelDescription_SIRX;
                  initialValues = initialValues_SIRX;
                  modelFunction = modelFunction_SIRX;
               }

               else if (*(uint16_t *)optarg == *(uint16_t *)"ERF")
               {
                  modelDescription = modelDescription_ERF;
                  initialValues = initialValues_ERF;
                  modelFunction = modelFunction_ERF;
               }

               else if (*(uint16_t *)optarg == *(uint16_t *)"GLF")
               {
                  modelDescription = modelDescription_GLF;
                  initialValues = initialValues_GLF;
                  modelFunction = modelFunction_GLF;
               }

               else
                  return usage(exe);
            else
               return usage(exe);
            break;

         case 'e':
            do_simulation = false;
            do_curve_fit  = false;
            export_series = true;
            break;

         case 'r':
            do_simulation = true;
            export_series = false;
            break;

         case 's':
            do_simulation = true;
            do_curve_fit  = false;
            break;

         case 't':
            if (optarg && *optarg && (v = strtold(optarg, NULL)) > 0.0L)
               thr = v;
            else
               return usage(exe);
            break;

         case 'o':
            if (optarg && *optarg
             && ((v = strtold(optarg, &chk)) != 0.0L || chk != optarg)
             && (v < z || isnan(z)))
               o = v;
            else
               return usage(exe);
            break;

         case 'z':
            if (optarg && *optarg
             && ((v = strtold(optarg, &chk)) != 0.0L || chk != optarg)
             && (v >= o || isnan(o)))
               z = v;
            else
               return usage(exe);
            break;

         case 'h':
         case '?':
         default:
            return usage(exe);
            break;
      }
   }

   if (argc-optind != 3 || argc && argv[optind][0] == '?' && argv[optind][1] == '\0')
      return usage(exe);

   FILE  *csv, *tsv;
   char  *country = argv[optind];
   size_t cl = strlen(country);

   if (*country)
      if (csv = (*(uint16_t *)argv[optind+1] == *(uint16_t *)"-")
                ? stdin
                : fopen(argv[optind+1], "r"))
      {
         if (tsv = (*(uint16_t *)argv[optind+2] == *(uint16_t *)"-")
                   ? stdout
                   : fopen(argv[optind+2], "w"))
         {
            int     i, m = 0, n = 0;
            char   *line;
            char    d[65536];
            ldouble t[ndays], c[ndays], l[ndays];

            // day 1 of the CSSE time series is 2020-01-22, hence day 0 is 2020-01-21, i.e. t[20]
            for (i = 0; i < ndays; i++)
               t[i] = i - 20, c[i] = l[i] = NAN;

            // find the lines with the series of the specified country
            while (line = fgets(d, 65536, csv))
            {
               size_t len = strlen(line);
               size_t col = collen(line);
               if (len != col && (line = strstr(line+col+1, country)) && line[cl] == ',')
               {
                  // skip 3 more fields (Country/Region,Lat,Long)
                  for (i = 0; i < 3; i++)
                     line += collen(line) + 1;

                  // read the case numbers and sum it up into the c array, starting at day 1 = index 20+1;
                  int p, q = p = 21;
                  while (*line && *line != '\n' && *line != '\r')
                  {
                     char   *chk;
                     ldouble v = strtod(line, &chk);
                     if (chk > line && isfinite(v))
                     {
                        if (isnan(c[q]))
                           c[q]  = v;
                        else
                           c[q] += v;
                        line = chk + 1;
                     }
                     else
                        line += collen(line) + 1;
                     q++;
                  }

                  if (q > p)
                  {
                     if (p > m) m = p;
                     if (q > n) n = q;
                  }
               }
            }

            if (isfinite(o))
            {
               for (i = 0; t[i] != o && i < ndays; i++);
               if (i < ndays && isfinite(c[i]))
                  m = i;
               else
                  return usage(exe);
            }
            else
               for (; isfinite(c[m]) && c[m] < thr && m < n; m++);

            if (isfinite(z))
            {
               for (i = ndays-1; t[i] != z && i >= 0; i--);
               if (i >= 0 && isfinite(c[i]))
                  n = i + 1;
               else
                  return usage(exe);
            }

            if (n > m)
            {
               fprintf(tsv, "# %s", exe);
               for (i = 1; i < argc; i++)
                  fprintf(tsv, " %s", argv[i]);
               fprintf(tsv, "\n# %s\n", country);

               if (do_simulation)
               {
                  int j, k = initialValues(t[m], c[m], c[n-1], A, f);
                  ldouble chiSqr = NAN;

                  if (!do_curve_fit
                   || isfinite(chiSqr = curveFit(n - m, &t[m], &c[m], k, A, dA, f, modelFunction)))
                  {
                     fprintf(tsv, "%s\n#\n", modelDescription);
                     for (j = 0; j < mpar && isfinite(A[j]); j++)
                        if (dA[j] != 0.0L)
                           fprintf(tsv, "#      %c%i = %11.6Lg ± %.5Lg %%\n", 'a', j, A[j], dA[j]);
                        else
                           fprintf(tsv, "#      %c%i = %11.6Lg\n", 'a', j, A[j]);

                     if (isfinite(chiSqr))
                        fprintf(tsv, "#\n#  ChiSqr = %11.1Lf\n", chiSqr);
                  }
                  else if (do_curve_fit)
                     fprintf(tsv, "# Curve fit failed\n");

                  if (export_series)
                  {
                     // write the column header with formular symbols and units.
                     // - the formular symbol of time is 't', the unit symbol of day is 'd'
                     // - the formular symbol of the number of cases is C without a unit
                     // - the formular symbol of the siumulated model is L without a unit
                     fprintf(tsv, "t/d\tC\tL\n");
                     for (i = 0; i < ndays; i++)
                     {
                        if (i >= m)
                           modelFunction(t[i], &l[i], A, i == m);
                        else
                           l[i] = 0.0L;

                        if (isfinite(c[i]) && isfinite(l[i]))
                           fprintf(tsv, "%.0Lf\t%.0Lf\t%.6Lf\n", t[i], c[i], l[i]);
                        else if (isfinite(c[i]))
                           fprintf(tsv, "%.0Lf\t%.0Lf\t*\n",     t[i], c[i]);
                        else if (isfinite(l[i]))
                           fprintf(tsv, "%.0Lf\t*\t%.6Lf\n",     t[i],       l[i]);
                        else
                           fprintf(tsv, "%.0Lf\t*\t*\n",         t[i]);
                     }
                  }
               }

               else
               {
                  // write the column header with formular symbols and units.
                  // - the formular symbol of time is 't', the unit symbol of day is 'd'
                  // - the formular symbol of the number of cases is C without a unit
                  fprintf(tsv, "t/d\tC\n");
                  for (i = 0; i < ndays; i++)
                  {
                     if (isfinite(c[i]))
                        fprintf(tsv, "%.0Lf\t%.6Lf\n",t[i], c[i]);
                     else
                        fprintf(tsv, "%.0Lf\t*\n",    t[i]);
                  }
               }
            }
            else
               fprintf(tsv, "# No values for country %s encountered.\n", country);

            if (tsv != stdout)
               fclose(tsv);
         }

         if (csv != stdin)
            fclose(csv);
      }

   return 0;
}
