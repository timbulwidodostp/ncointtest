{smcl}
{* 21jan2021}{...}
{hline}
help for {hi:ncointtest}
{hline}

{title:ncointtest computes the non-cointegration test between two time series in the frequency domain, as in Souza et al (2018)}

{title:Syntax}

{p 8 16 2}{cmd:ncointtest}
{it:varlist}
{cmd:,} {cmd:d(}{it:0<d<=1}{cmd:)}
{cmd:b(}{it:0<b<1}{cmd:)}
{cmd:r()}{#)}

{pstd }{cmd:ncointest} is for use with time-series data. You must {cmd:tsset} your data before using {cmd:cointest}; 
see help {cmd:tsset}.

{title:Description}

{pstd}{cmd:ncointest} applies the non-cointegration test in the frequency domain on a bivariate series {it:varlist}. This test considers 
that a cointegration process stablishes once the order of cointegration of the residuals of the linear combination of the series is of a lower values
than the original series. In the way, suppose {it:d} is the order of integration of the bivariate system, then d_{1} is the order of integration of the residuals 
of one series over the other, with 0<d_{1}<{it:d}<=1, that is, the two series share commonalities. {it:d} is the parmeter {it:d} in the test, which has no default and can be estimated using the {cmd:gphudak} package. Furthermore, this command also requires that {cmd:more_mata} package to be installed in the machine. The parameter {it:b} is the ratio of the numbers of observations used in this test, with default 0.8. {it:r} is the weight applied to neghbouring frequencies in the test. ncointtest does not allow for missing values in the series. 


{title:Options}

{phang}{cmd:d} is the shared order of the cointegration of the bivariate system, which as Leschinski et al (2020) can be the simple average of the order of cointegration of each time series, given that statistically they are similar.  

{phang}{cmd:b} is the bandwidth in the test and it is a percentage of the number of observations, with n^{it:b} being the number of observations used in the test. 

{phang}{cmd:r} is the number of neighboring frequencies to be utilized in estimating the density spectral matrix. 

 
{title:Examples}

{pstd}{stata "import delimited https://raw.githubusercontent.com/alanleal-econ/ncointtest/main/ncointtest_data.csv, clear" :. import delimited "https://raw.githubusercontent.com/alanleal-econ/ncointtest/main/ncointtest_data.csv", clear}

{pstd}{stata "gen mydate=ym(1996,02)+_n" :. gen mydate=ym(1996,02)+_n}

{pstd}{stata "tsset mydate" :. tsset mydate}

{pstd}{stata "ncointtest ldj lftse, d(.923327) b(0.8) r(3)" :. ncointtest ldj lftse, d(.923327) b(0.8) r(3)}

{pstd}{stata "ncointtest ldj lftse, d(.923327) b(0.7) r(3)" :. ncointtest ldj lftse, d(.923327) b(0.7) r(3)}

{pstd}{stata "ncointtest ldj lftse, d(.923327) b(0.8) r(3)" :. ncointtest ldj lftse, d(.923327) b(0.8) r(4)}

{pstd}{stata "ncointtest ldj lftse, d(.923327) b(0.7) r(3)" :. ncointtest ldj lftse, d(.923327) b(0.7) r(4)}



{title:Author}

{pstd} Alan Leal, University of São Paulo, Brazil{break} 
       prof@alanleal-econ.com
       
{title:References}

{pstd} Igor Viveiros Melo Souza, Valderio Anselmo Reisen, Glaura da Conceição Franco & Pascal Bondon (2018) The Estimation and Testing of the Cointegration Order Based on the Frequency Domain, Journal of Business & Economic Statistics, 36:4, 695-704, DOI: 10.1080/07350015.2016.1251442 

{pstd} Leschinski, C., Voges, M. & Sibbertsen, P. (2020) A comparison of semiparametric tests for fractional cointegration. Stat Papers. https://doi.org/10.1007/s00362-020-01169-1

{title:Disclaimer}


{pstd} This program is provided without warranty of any kind. The author is not responsible for any cost derived by the usage
 of this program.
 
{title:Acknowledgements} 

{pstd} This command was translated from code belonging to Prof. Igor Viveiro Melo Souza, who graciouly allowed for its translation. Mata function for fractionary differentiation was based on freely distributed code in the fracdiff package in R. Special thanks to Jorge Perez, whose discret fourier transform allowed for the fast translation of this code and also gave his authorization for use. 


