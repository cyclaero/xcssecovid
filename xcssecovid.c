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

enum
{
    no_error,
   usg_error,
    fs_error,
   val_error,
   fit_error,
   par_error,
   sim_error
};

int usage(const char *exe)
{
   printf("\nExtract, curve fit an epidemiological model and transpose CSSE@JHU's Covid-19 cases data per country - Copyright Dr. Rolf Jansen (c) 2020 - %s\n", version);
   printf("Usage: %s [-a<0-9> value] [-f (0-9)+] [-m model] [-e] [-i] [-s] [-o day#] [-z day#] [-r depth] [-h|-?|?] <Country> <CSV Input file> <TSV Output File>\n\n\
       -a<0..9> value      Optionally set initial values for the model's parameters the Differential Equation Solver and Curve Fitting.\n\
                           The models deduce its initial parameters from the boundaries of the imported time series and by common\n\
                           knowledge/best educated guesses. Example: -a0 0.57 -a1 125000\n\n\
       -f (0-9)+           Overrides the default selection of a model's parameters which take part in curve fitting, and which usually is a0, a1, a2.\n\
                           Example: -f 1245 would lead to curve fitting against the paramters a1, a2, a4 and a5, while a0 and a3 would be left untouched.\n\
                           Different models got different number of parameters, which is currently 3 to 9.\n\n\
       -m model            Select the model for curve fitting and simulation:\n\
                           - LF   Logistic Function -- https://en.wikipedia.org/wiki/Logistic_function\n\
                           - LDE  Logistic Differential Equation -- https://en.wikipedia.org/wiki/Logistic_function#Logistic_differential_equation\n\
                           - SI   Epidemiological SI-Model (basically another form of the LDE) -- https://de.wikipedia.org/wiki/SI-Modell\n\
                           - SIR  Epidemiological SIR-Model [default] -- https://en.wikipedia.org/wiki/Mathematical_modelling_of_infectious_disease#The_SIR_model\n\
                           - SEIR Epidemiological SEIR-Model -- https://www.idmod.org/docs/hiv/model-seir.html\n\
                           - SIRX Extension of the SIR-Model by the HU/RKI -- http://rocs.hu-berlin.de/corona/docs/forecast/model/\n\
                           - ERF  Shifted Error Function -- (generic form) https://en.wikipedia.org/wiki/Error_function\n\
                           - GLF  Generalised Logistic Function -- https://en.wikipedia.org/wiki/Generalised_logistic_function\n\n\
       -e                  Only export the extracted and transposed time series without curve fitting and simulation of the model.\n\n\
       -i                  Only inform the model description and the values of its parameters with or without fitting, but without exporting any curve data.\n\n\
       -s                  Simulate the model without curve fitting before.\n\n\
       -t thresh           Threshold of minimal number of cases to include into curve fitting and simulation [default: 17].\n\n\
       -o day#             The day# of the first data point in the imported time series to be included for curve fitting.\n\
                           [default: first day with more than <thresh> cases].\n\n\
       -z day#             The day# of the last data point in the imported time series to be included for curve fitting.\n\
                           [default: last day of the imported series].\n\n\
       -r depth            Retrospective day by day curve fitting and simulation of the model back for depth number of days.\n\n\
       -h|-?|?             Show these usage instructions.\n\n\
       <Country>           Select the country for which the time series shall be processed.\n\n\
       <CSV Input file>    Path to the CSSE@JHU's Covid-19 time series CSV file.\n\n\
       <TSV Output File>   Path to the TSV output file containing the extracted and transposed time series for the given <Country>, including\n\
                           one ore more columns for the simulated time series by the given model, using the parameters as resulted from curve fitting.\n\n", exe);
   return usg_error;
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
   int rc = no_error;

   const char *exe = execname(argv[0]);

   char *modelDescription = modelDescription_SIR;
   initvals initialValues = initialValues_SIR;
   function modelFunction = modelFunction_SIR;

   bool do_simulation = true,
        do_curve_fit  = true,
        export_series = true;

   int     r = 0;
   ldouble o = NAN,
           z = NAN,
         thr = 17.0L;

   ldouble A0[mpar] = {NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN},
         **A        = NULL;

   int    f[mpar] = {undefined, undefined, undefined, undefined, undefined,
                     undefined, undefined, undefined, undefined, undefined};

   bool aopt = false;
   int  opc;
   while ((opc = getopt(argc, argv, "af:0:1:2:3:4:5:6:7:8:9:m:eist:o:r:z:h?")) != -1)
   {
      char   *chk, c;
      int     i;
      long    l;
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
               if ('0' <= opc && opc <= '9' && isnan(A0[opc -= '0'])
                && ((v = strtold(optarg, &chk)) != 0.0L || chk != optarg)
                && isfinite(v))
               {
                  A0[opc] = v;
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

               else if (*(uint32_t *)optarg == *(uint32_t *)"SEIR" && optarg[4] == '\0')
               {
                  modelDescription = modelDescription_SEIR;
                  initialValues = initialValues_SEIR;
                  modelFunction = modelFunction_SEIR;
               }

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

         case 'i':
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

         case 'r':
            if (optarg && *optarg && (l = strtol(optarg, NULL, 10)) > 0)
               r = (int)l;
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
            int     h, i, j, m = 0, n = 0;
            char   *line;
            char    d[65536];
            ldouble t[ndays],
                    c[ndays],
                  **l = malloc((r+1)*sizeof(ldouble *));

            for (h = 0; h <= r; h++)
               l[h] = malloc(ndays*sizeof(ldouble));

            // day 1 of the CSSE time series is 2020-01-22, hence day 0 is 2020-01-21, i.e. t[20]
            for (i = 0; i < ndays; i++)
               t[i] = i - 20, c[i] = NAN;

            for (h = 0; h <= r; h++)
               for (i = 0; i < ndays; i++)
                  l[h][i] = NAN;

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
               if (r < n-m - 3)
               {
                  fprintf(tsv, "# %s", exe);
                  for (i = 1; i < argc; i++)
                     fprintf(tsv, " %s", argv[i]);
                  fprintf(tsv, "\n#\n# %s\n", country);

                  if (do_simulation)
                  {
                     A = malloc((r+1)*sizeof(ldouble *));
                     for (h = 0; h <= r; h++)
                        A[h] = malloc(mpar*sizeof(ldouble));

                     for (j = 0; j < mpar; j++)
                        A[0][j] = A0[j];
                     int k = initialValues(t[m], c[m], c[n-1], A[0], f);

                     if (do_curve_fit)
                     {
                        ldouble chiSqr = NAN;
                        ldouble dA[mpar] = {};
                        if (isfinite(chiSqr = curveFit(n-m, &t[m], &c[m], k, A[0], dA, f, modelFunction)))
                        {
                           fprintf(tsv, "%s\n#\n", modelDescription);
                           for (j = 0; j < mpar && isfinite(A[0][j]); j++)
                              if (dA[j] != 0.0L)
                                 fprintf(tsv, "#      %c%d = %11.6Lg ± %.5Lg %%\n", 'a', j, A[0][j], dA[j]);
                              else
                                 fprintf(tsv, "#      %c%d = %11.6Lg\n", 'a', j, A[0][j]);

                           if (isfinite(chiSqr))
                              fprintf(tsv, "#\n#  ChiSqr = %11.1Lf\n", chiSqr);
                        }

                        else
                        {
                           fprintf(tsv,
                                "# Curve fit failed\n");
                           printf("Curve fit failed\n");
                           rc = fit_error;
                           goto cleanup;
                        }

                        if (rc == no_error && r)
                        {
                           for (h = 1; h <= r; h++)
                           {
                              for (j = 0; j < mpar; j++)
                                 A[h][j] = A0[j], dA[j] = 0.0L;
                              initialValues(t[m], c[m], c[n-1 - h], A[h], f);

                              if (isfinite(chiSqr = curveFit(n-m - h, &t[m], &c[m], k, A[h], dA, f, modelFunction)))
                              {
                                 fprintf(tsv, "#\n# Retropective curve fit #%d:\n%s\n#\n", h, modelDescription);
                                 for (j = 0; j < mpar && isfinite(A[h][j]); j++)
                                    if (dA[j] != 0.0L)
                                       fprintf(tsv, "#      %c%d = %11.6Lg ± %.5Lg %%\n", 'a', j, A[h][j], dA[j]);
                                    else
                                       fprintf(tsv, "#      %c%d = %11.6Lg\n", 'a', j, A[h][j]);

                                 if (isfinite(chiSqr))
                                    fprintf(tsv, "#\n#  ChiSqr = %11.1Lf\n", chiSqr);
                              }

                              else
                              {
                                 fprintf(tsv,
                                      "# Retrospective curve fit #%d failed\n", h);
                                 printf("Retrospective curve fit #%d failed\n", h);
                                 rc = fit_error;
                                 break;
                              }
                           }

                           fprintf(tsv, "#\n# Std. Dev. of all parameters of the retrospection over %d days\n", r+1);
                           ldouble *Z = malloc((r+1)*sizeof(ldouble));
                           for (j = 0; j < mpar && isfinite(A[0][j]); j++)
                           {
                              for (h = 0; h <= r; h++)
                                 Z[h] = A[h][j];
                              if (dA[j] != 0.0L)
                                 fprintf(tsv, "#     s%c%d = %11.6Lg %%\n", 'a', j, sdev(r+1, Z)/aave(r+1, Z)*100.0L);
                           }
                           free(Z);
                        }
                     }

                     if (rc == no_error && export_series)
                     {
                        // write the column header with formular symbols and units.
                        // - the formular symbol of time is 't', the unit symbol of day is 'd'
                        // - the formular symbol of the number of cases is C without a unit
                        // - the formular symbol of the siumulated model is L without a unit
                        fprintf(tsv, "t/d\tC\tL");
                        if (!r)
                           fprintf(tsv, "\n");
                        else
                        {
                           for (h = 0; h < r; h++)
                              fprintf(tsv, "%d\tL", h);
                           fprintf(tsv, "%d\n", h);
                        }

                        for (h = 0; h <= r; h++)
                           for (i = 0; i < ndays; i++)
                           {
                              if (i >= m)
                              {
                                 if (modelFunction(t[i], &l[h][i], A[h], i == m) != no_error)
                                    goto sim_err_out;
                              }
                              else
                                 l[h][i] = 0.0L;
                           }

                        for (i = 0; i < ndays; i++)
                        {
                           if (isfinite(c[i]))
                              fprintf(tsv, "%.0Lf\t%.0Lf", t[i], c[i]);
                           else
                              fprintf(tsv, "%.0Lf\t*",     t[i]);

                           for (h = 0; h <= r; h++)
                              if (isfinite(l[h][i]))
                                 fprintf(tsv, "\t%.0Lf",   l[h][i]);
                              else
                                 fprintf(tsv, "\t*");
                           fprintf(tsv, "\n");
                        }

                        goto cleanup;

                     sim_err_out:
                        fprintf(tsv,
                             "# Simulation failed\n");
                        printf("Simulation failed\n");
                        rc = sim_error;
                     }

                  cleanup:
                     if (A)
                     {
                        for (h = 0; h <= r; h++)
                           free(A[h]);
                        free(A);
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
                           fprintf(tsv, "%.0Lf\t%.6Lf\n", t[i], c[i]);
                        else
                           fprintf(tsv, "%.0Lf\t*\n",     t[i]);
                     }
                  }

                  for (h = 0; h <= r; h++)
                     free(l[h]);
                  free(l);
               }

               else
               {
                  fprintf(tsv,
                       "# The retrospective depth #%d is too large\n", r);
                  printf("The retrospective depth #%d is too large\n", r);
                  rc = par_error;
               }

            else
            {
               fprintf(tsv,
                    "# No values for country %s encountered.\n", country);
               printf("No values for country %s encountered.\n", country);
               rc = val_error;
            }

            if (tsv != stdout)
               fclose(tsv);
         }

         else
         {
            printf("Output file %s could not be opened for writing.\n", argv[optind+2]);
            rc = fs_error;
         }

         if (csv != stdin)
            fclose(csv);
      }

      else
      {
         printf("Input file %s not found.\n", argv[optind+1]);
         rc = fs_error;
      }

   return rc;
}
