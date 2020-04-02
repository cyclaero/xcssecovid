# xcssecovid
Extract, curve fit an epidemiological model and transpose CSSE@JHU's Covid-19 cases data per country

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

    Extract, curve fit an epidemiological model and transpose CSSE@JHU's Covid-19 cases data per country - Copyright Dr. Rolf Jansen (c) 2020 - Version 1.0.1
    
    Usage: xcssecovid [-a<0-9> value] [-f (0-9)+] [-m model] [-e] [-r] [-s] [-o day#] [-z day#] [-h|-?|?] <Country> <CSV Input file> <TSV Output File>
   
       -a<0..9> value      Optionally set initial values for model's parameters the Differential Equation Solver and Curve Fitting.
                           The models deduce its initial parameters from the boundaries of the imported time series and by common
                           knowledge/best educated guesses. Example: -a0 0.57 -a1 125000

       -f (0-9)+           Overrides the default selection of a model's parameters which take part in curve fitting, and which usually is a0, a1, a2.
                           Example: -f 1245 would lead to curve fitting against the paramters a1, a2, a4 and a5, while a0 and a would be left untouched.
                           Different models got different number of parameters, which is currently 3 to 6.

       -m model            Select the model for curve fitting and simulation:
                           - LF   Logistic Function -- https://en.wikipedia.org/wiki/Logistic_function
                           - LDE  Logistic Differential Equation -- https://en.wikipedia.org/wiki/Logistic_function#Logistic_differential_equation
                           - SI   Epidemiological SI-Model (basically another form of the LDE) -- https://de.wikipedia.org/wiki/SI-Modell
                           - SIR  Epidemiological SIR-Model [default] -- https://en.wikipedia.org/wiki/Mathematical_modelling_of_infectious_disease#The_SIR_model

       -e                  Only export the extracted and transposed time series without curve fitting and simulation of the model.

       -r                  Only report the model description and the values of its parameters with or without fitting, but without exporting any curve data.

       -s                  Simulate the model without curve fitting before.

       -o day#             The day# of the first data point in the imported time series to be included for curve fitting.
                           [default: first day with more than 17 cases].

       -z day#             The day# of the last data point in the imported time series to be included for curve fitting.
                           [default: last day of the imported eries].

       -h|-?|?             Show these usage instructions.

       <Country>           Select the country for which the time series shall be processed.

       <CSV Input file>    Path to the CSSE@JHU's Covid-19 time series CSV file.

       <TSV Output File>   Path to the TSV output file containing the extracted and transposed time series for the given <Country>, including
                           a column for a simulated time series by the given model, using the parameter as resulted from curve fitting.
   
## Use cases on my BLog (German)
Zeitreihenanalyse zur Ausbreitung von Covid-19 ausserhalb Chinas: https://obsigna.com/articles/1584931539.html
English translation by MS-Online-Translation: https://www.translatetheweb.com/?from=de&to=en&a=https://obsigna.com/articles/1584931539.html
