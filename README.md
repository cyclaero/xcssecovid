# xcssecovid
Extract, curve-fit a logistic function and transpose JHU-CSSE's Covid-19 cases data per country/region

## Usage

  1. Compile `xcssecovid.c` on either of FreeBSD, macOS or Linux:

     `cc -g0 -O3 -march=native xcssecovid.c -Wno-parentheses -lm -o xcssecovid`

  2. Download the daily updated time series of confirmed Covid-19 cases
     from CSSE's (at Johns Hopkins University) GitHub site CSSE COVID-19 Dataset  
     https://github.com/CSSEGISandData/COVID-19/tree/master/csse_covid_19_data

     `curl -O https://github.com/CSSEGISandData/COVID-19/raw/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv`

  3. Extract the time series row of a given country starting on 2020-01-21 as day 0
     and write it out together with the country's cases and a simulated curve-fitted
     logistic function into a three-column TSV output file:

     `./xcssecovid Germany time_series_covid19_confirmed_global.csv Germany.tsv`

  4. Open the TSV file with your favorite graphing and/or data analysis application.  
     
## Use cases on my BLog (German)
Zeitreihenanalyse zur Ausbreitung von Covid-19 ausserhalb Chinas: https://obsigna.com/articles/1584931539.html
English translation by MS-Online-Translation: https://www.translatetheweb.com/?from=de&to=en&a=https://obsigna.com/articles/1584931539.html
