*! ncointtest computes the non-cointegration test between two time series in the frequency domain, as in Souza et al (2018)
*! v.1.0 January 19, 2021 Alan Leal (prof@alanleal-econ.com)

* It's necessary installing the  mm_seq function from more_mata package!
cap which moremata.hlp
if _rc {
ssc install moremata
}

program define ncointtest, rclass
	version 16
	syntax varlist(ts min=2 max=2) [,d(real 0) b(real 0.7) r(integer 3)]
	* Checking whether the data is time-set or not
	_ts
	tokenize `varlist'
	if (`d'<=0 | `d'>1){
		di in red "d must be a value between 0 and 1!"
		exit
	}
	if (`b'<0 | `b'>1){
		di in red "Bandwidth must be a value betweeen 0 and 1!"
		exit
	}
	if (`r'<=1){
		di in red "r must be larger than 1!"
		exit
	}
	mata: x=st_data(., ("`1'"))
	mata: y=st_data(., ("`2'"))
	mata: mis=sum(rowmissing((x,y)) :== 1)
	mata: st_numscalar("mis",mis)
	if (mis>0){
		di in red "Sample contains missing values!"
		exit
	}
	mata: st_matrix("r(b)",ncointtest(x,y,`d',`b',`r'))
	matrix r_e=r(b)
	di  "Non-cointegration test in the frequency domain"
	di as text "{hline 72}"
	di  "H0: `1' and `2' are non-cointegrated"
	di  "Note: Asymptotic values are based on Z ~ N(0,1)"
	di as text "{hline 72}"
	di   "Set up: d = `d', Bandwidth = `b', r = `r' and n = " (_N-1)
	di as text "{hline 72}"
	di  "         Est.b" "          ass.d"  "        t-like stat" "      p-value"
	di as text "{hline 72}"
	di %15.4f r_e[1,1]  %15.4f r_e[1,2]  %15.4f r_e[1,3] %15.4f r_e[1,4]
	di as text "{hline 72}"
end


// *! fracdiff_my computes the fractional differenced series for a d, with 0<=d<=1 - necessary to the ncointtest
// * This is a Stata translation of the code used for the fracdiff package in R, as written by Valdério Heisen
mata
	version 16
	real vector fracdiff_my(real vector x_f,real scalar b){
	real scalar X_r
	real vector PI, y_diff,x
	X_r=rows(x_f)
	x=x_f:-(sum(x_f)/X_r)
	PI=J(X_r,1,0)
	PI[1,1]=-b
	for (k=2;k<=X_r;k=k+1){
		PI[k,1]=PI[k-1,1]*(k-1-b)/k
	}
	y_diff=x
	for (i=2;i<=X_r;i=i+1){
		y_diff[i,1]=x[i,1]+sum(PI[(1..(i-1)),1]:*x[((i-1)..1),1])
	}
	return(y_diff)
	}
end


// *! dft computes the discrete fourier transform
// * This was copied by a post in StataList from Jorge Perez (who graciously shared his code for this command), available at https://www.stata.com/statalist/archive/2011-01/msg00088.html
mata
	version 16
	complex matrix dft( transmorphic matrix x)
	{	
		real matrix yr, yi, dftr, dfti, range
		real scalar n, i, ae
		complex matrix y, df
		y=C(x)
		yr=Re(y)
		yi=Im(y)
		n=rows(x)
		dftr=J(n,1,.)
		dfti=J(n,1,.)
		range=range(1,n,1)
		ae=-2*pi()/n
		for (i=0; i<=n-1; i++) {	
			dftr[i+1]=colsum(yr:*cos((range:-1):*ae:*i)-yi:*sin((range:-1):*ae:*i))
			dfti[i+1]=colsum(yr:*sin((range:-1):*ae:*i)+yi:*cos((range:-1):*ae:*i))
		}
		df=conj(C(dftr,dfti))
		return(df)
	}
end

*! ncointtest computes the non-cointegration test between two time series in the frequency domain, as in Souza et al (2018)
mata
version 16
real vector ncointtest(real vector x, real vector y, real fracdiff_v ,real scalar band, real scalar r){
	real matrix A,x_reg
	real vector di11, di12, di21, di22, L_det, freq,vec_1,y_reg
	real scalar m, band_temp,b,sd_as, t_ldr,p_value
	if (fracdiff_v~=0){
		x=fracdiff_my(x,fracdiff_v)
		x=x[(2..rows(x)),1]
		y=fracdiff_my(y,fracdiff_v)
		y=y[(2..rows(y)),1]
	};
	di11=((1/(2*pi()*rows(x)))*((dft(x):*conj(dft(x)))))[(2..rows(x)),1]
	di12=((1/(2*pi()*rows(x)))*((dft(x):*conj(dft(y)))))[(2..rows(x)),1]
	di21=((1/(2*pi()*rows(x)))*((dft(y):*conj(dft(x)))))[(2..rows(x)),1]
	di22=((1/(2*pi()*rows(x)))*((dft(y):*conj(dft(y)))))[(2..rows(x)),1]
	m=band
	band_temp=floor(rows(x)^band)/r
	if (mod(band_temp,0.5)==0){
		band=floor(band_temp)*r
	}
	else{
		band=round(band_temp)*r
	}
	A=J(band,rows(di11),0)
	for (i=1;i<=(rows(A)/r);i=i+1){
		A[i,(r*i-r+1..r*i-r+1+r-1)]=J(1,r,1)
	}
	L_det=Re(((A*di11):*(A*di22))-((A*di12):*(A*di21)))
	L_det=log(L_det[(1..(rows(A)/r)),1])
	freq=(log(2:-2*cos(2*pi()*(mm_seq(2,band,r)')/rows(x))))'
	vec_1=J(rows(freq),1,1)
	x_reg=(vec_1,freq)
	y_reg=L_det
	b=((pinv(x_reg'*x_reg)*x_reg'*y_reg))[2,1]
	sd_as=sqrt((trigamma(3)+trigamma(2)):/sum((freq:-(sum(freq)/rows(freq))):^2))
	t_ldr=b/sd_as
	p_value=1-normal(abs(t_ldr))
	return ((b,sd_as,t_ldr,p_value))
}
end
