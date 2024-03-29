### [ACTION REQUIRED] Your GitHub account, cyclaero, will soon require 2FA

Here is the deal: https://obsigna.com/articles/1693258424.html

---
 
# xcssecovid
Extract, curve fit an epidemiological model and transpose CSSE@JHU's Covid-19 case data per country

## Usage

  1. Compile `numerics.c`, `models.c` and `xcssecovid.c` on either of FreeBSD or macOS:

     `clang -g0 -O3 -march=native numerics.c models.c xcssecovid.c -Wno-parentheses -lm -o xcssecovid`

  2. Download the daily updated time series of confirmed Covid-19 cases
     from CSSE's (at Johns Hopkins University) GitHub site CSSE COVID-19 Dataset
     https://github.com/CSSEGISandData/COVID-19/tree/master/csse_covid_19_data

     `curl -O https://github.com/CSSEGISandData/COVID-19/raw/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv`

  3. Extract the time series row of a given country starting on 2020-01-21 as day 0
     and write it out together with the country's cases and a simulated curve-fitted
     logistic function into a three-column TSV output file:

     `./xcssecovid Germany time_series_covid19_confirmed_global.csv Germany.tsv`

  4. Open the TSV file with your favorite graphing and/or data analysis application
  
### All the options:

`./xcssecovid -h`

    Extract, curve fit an epidemiological model and transpose CSSE@JHU's Covid-19 case data per country - Copyright Dr. Rolf Jansen (c) 2020 - Version 1.0.3
    
    Usage: xcssecovid [-a<0-9> value] [-f (0-9)+] [-m model] [-e] [-i] [-s] [-o day#] [-z day#] [-r depth] [-h|-?|?] <Country> <CSV Input file> <TSV Output File>

       -a<0..9> value      Optionally set initial values for the model's parameters the Differential Equation Solver and Curve Fitting.
                           The models deduce its initial parameters from the boundaries of the imported time series and by common
                           knowledge/best educated guesses. Example: -a0 0.57 -a1 125000

       -f (0-9)+           Overrides the default selection of a model's parameters which take part in curve fitting, and which usually is a0, a1, a2.
                           Example: -f 1245 would lead to curve fitting against the paramters a1, a2, a4 and a5, while a0 and a3 would be left untouched.
                           Different models got different number of parameters, which is currently 3 to 9.

       -m model            Select the model for curve fitting and simulation:
                           - LF   Logistic Function -- https://en.wikipedia.org/wiki/Logistic_function
                           - LDE  Logistic Differential Equation -- https://en.wikipedia.org/wiki/Logistic_function#Logistic_differential_equation
                           - SI   Epidemiological SI-Model (basically another form of the LDE) -- https://de.wikipedia.org/wiki/SI-Modell
                           - SIR  Epidemiological SIR-Model [default] -- https://en.wikipedia.org/wiki/Mathematical_modelling_of_infectious_disease#The_SIR_model
                           - SEIR Epidemiological SEIR-Model -- https://www.idmod.org/docs/hiv/model-seir.html
                           - SIRX Extended SIR-Model by the HU/RKI -- http://rocs.hu-berlin.de/corona/docs/forecast/model/
                           - ERF  Shifted Error Function -- (generic form) https://en.wikipedia.org/wiki/Error_function
                           - GLF  Generalised Logistic Function -- https://en.wikipedia.org/wiki/Generalised_logistic_function

       -e                  Only export the extracted and transposed time series without curve fitting and simulation of the model.

       -i                  Only inform the model description and the values of its parameters with or without fitting, but without exporting any curve data.

       -s                  Simulate the model without curve fitting before.

       -t thresh           Threshold of minimal number of cases to include into curve fitting and simulation [default: 17].

       -o day#             The day# of the first data point in the imported time series to be included for curve fitting.
                           [default: first day with more than 17 cases].

       -z day#             The day# of the last data point in the imported time series to be included for curve fitting.
                           [default: last day of the imported series].

       -q                  Output the time series of the whole set of the simulated differential equations of the given model along with the extracted CSSE@JHU's Covid-19 case data.\n\n\

       -r depth            Retrospective day by day curve fitting and simulation of the model back for depth number of days.\n\
                           (Note: the -q and -r options are mutually exclusive, the last option on the command line counts).\n\n\

       -h|-?|?             Show these usage instructions.

       <Country>           Select the country for which the time series shall be processed.

       <CSV Input file>    Path to the CSSE@JHU's Covid-19 time series CSV file.

       <TSV Output File>   Path to the TSV output file containing the extracted and transposed time series for the given <Country>, including
                           one or more columns for the simulated time series by the given model, using the parameters as resulted from curve fitting.
   
## Use cases on my BLog (German)
Modelle zur Beschreibung der Ausbreitung des n-Coronavirus-2019: https://obsigna.com/articles/1586103085.html
English translation by MS-Online-Translation: https://www.translatetheweb.com/?from=de&to=en&a=https://obsigna.com/articles/1586103085.html
