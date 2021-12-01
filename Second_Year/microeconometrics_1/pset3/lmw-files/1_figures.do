clear all
set more off
local rootDir = "/Users/Research/LMW2019"

local data "`rootDir'/data/stata"
local temp "`rootDir'/data/temp"
local figDir "`rootDir'/figures"

local rate=87.94

********************************************************************************
* Figures
* Last updated: 23 Apr 2019
********************************************************************************

********************************************************************************
* Figure 1: The electric utility as a natural monopoly
********************************************************************************

* Note that there is no underlying data for these illustrative figures.

********************************************************************************
* Figure 2: Experimental evidence on the demand for rural electrification
********************************************************************************

***************
*** Panel A ***
***************

use `data'/demand.dta, clear
bysort price: egen takeup3=mean(takeup)
duplicates drop price, force
drop if price==.
keep price takeup3
gen rp=takeup3*100
drop takeup3
rename price kprice
sort kprice
save `temp'/rp.dta, replace

clear
set obs 15
gen kprice=.
replace kprice=0 in 1
replace kprice=3800 in 2
replace kprice=4000 in 3
replace kprice=5000 in 4
replace kprice=6000 in 5
replace kprice=7000 in 6
replace kprice=8000 in 7
replace kprice=9000 in 8
replace kprice=10000 in 9
replace kprice=12000 in 10
replace kprice=15000 in 11
replace kprice=17000 in 12
replace kprice=20000 in 13
replace kprice=25000 in 14
replace kprice=35000 in 15

gen price=kprice/`rate'

gen nes=.
replace nes=100 in 1
replace nes=100 in 2
replace nes=90 in 3
replace nes=90 in 4
replace nes=90 in 5
replace nes=80 in 6
replace nes=80 in 7
replace nes=70 in 8
replace nes=70 in 9
replace nes=60 in 10
replace nes=40 in 11
replace nes=30 in 12
replace nes=30 in 13
replace nes=20 in 14
replace nes=10 in 15

gen pap=.
replace pap=45 if kprice==15000
replace pap=20 if kprice==25000
replace pap=2.5 if kprice==35000
sort kprice

merge kprice using `temp'/rp.dta
drop if kprice==.
drop _merge
gen priorpoints=1
sort kprice

twoway ///
		(connected price nes, lcolor(gs10) lwidth(medthin) lpattern(longdash) msize(small) mcolor(gs8) msymbol(square)) ///
		(connected price rp, lcolor(gs5) lwidth(thick) lpattern(solid) msize(small) mcolor(gs2) msymbol(square)), ///
		xtitle("Take-up (%)", size(medium) margin(medsmall)) ///
		ytitle("Connection price (USD)", size(medium) margin(medsmall)) ///
		xlabel(0(20)100, labsize(medsmall)) ///
		ylabel(0(50)400, labsize(medsmall)) ///
		xsize(9) ysize(11) ///
		scheme(s1mono) ///
		legend(size(medsmall) ///
		label(1 "Kenyan gov't report") ///
		label(2 "Demand curve") ///
		symx(9) order(2 1) region(lstyle(none)) position(1) ring(0) rows(3))
		graph export `figDir'/figure2a.png, replace

rm `temp'/rp.dta
		
***************
*** Panel B ***
***************

use `data'/demand.dta, clear
drop if price==.
sum earn, detail
local l=`r(p25)'
local h=`r(p75)'
gen quartile=0
replace quartile=1 if earn<=`l'
replace quartile=4 if earn>=`h'
sum earn if quartile==1
sum earn if quartile==4
keep price earn quartile takeup 
bysort quartile price: egen take1=mean(takeup)
bysort price: egen take0=mean(takeup)
duplicates drop quartile price, force
drop takeup
sort quartile price

gen rp_low=take1*100 if quartile==1
gen rp_high=take1*100 if quartile==4
gen rp=take0*100 if quartile==4

rename price kprice
gen price=kprice/`rate'

twoway ///
		(connected price rp_low if quartile==1, lcolor(gs2) lwidth(medthin) lpattern(longdash) msize(small) mcolor(gs1) msymbol(square)) ///
		(connected price rp_high if quartile==4, lcolor(gs10) lwidth(medthin) lpattern(longdash) msize(small) mcolor(gs8) msymbol(square)) ///
		(connected price rp if quartile==4, lcolor(gs5) lwidth(thick) lpattern(solid) msize(small) mcolor(gs2) msymbol(square)), ///
		xtitle("Take-up (%)", size(medium) margin(medsmall)) ///
		ytitle("Connection price (USD)", size(medium) margin(medsmall)) ///
		xlabel(0(20)100, labsize(medsmall)) ///
		ylabel(0(50)400, labsize(medsmall)) ///
		xsize(9) ysize(11) ///
		scheme(s1mono) ///
		legend(size(medsmall) label( 1 "Monthly earnings, lower quartile") ///
		label(2 "Monthly earnings, upper quartile") ///
		label(3 "Demand curve") ///
		symx(9) order(3 1 2) region(lstyle(none)) margin(zero) position(1) ring(0) rows(3))
		graph export `figDir'/figure2b.png, replace

***************
*** Panel C ***
***************

use `data'/demand.dta, clear
bysort price wall: egen takeup3=mean(takeup)
bysort price: egen takeup4=mean(takeup)
duplicates drop price wall, force
keep price takeup3 takeup4 wall
gen temp=takeup3 if wall==1
bysort price: egen takeup5=max(temp)
keep if wall==0
drop if price==.
drop wall temp
gen rp_low = takeup3*100
gen rp_high = takeup5*100
gen rp = takeup4*100

rename price kprice
gen price=kprice/`rate'

twoway ///
		(connected price rp_low, lcolor(gs2) lwidth(medthin) lpattern(longdash) msize(small) mcolor(gs1) msymbol(square)) ///
		(connected price rp_high, lcolor(gs10) lwidth(medthin) lpattern(longdash) msize(small) mcolor(gs8) msymbol(square)) ///
		(connected price rp, lcolor(gs5) lwidth(thick) lpattern(solid) msize(small) mcolor(gs2) msymbol(square)), ///
		xtitle("Take-up (%)", size(medium) margin(medsmall)) ///
		ytitle("Connection price (USD)", size(medium) margin(medsmall)) ///
		xlabel(0(20)100, labsize(medsmall)) ///
		ylabel(0(50)400, labsize(medsmall)) ///
		xsize(9) ysize(11) ///
		scheme(s1mono) ///
		legend(size(medsmall) label( 1 "Low-quality walls subsample") ///
		label(2 "High-quality walls subsample") ///
		label(3 "Demand curve") ///
		symx(9) order(3 1 2) region(lstyle(none)) margin(zero) position(1) ring(0) rows(3))
		graph export `figDir'/figure2c.png, replace

********************************************************************************
* Figure 3: Experimental evidence on the social surplus implications of RE
********************************************************************************

***************
*** Panel A ***
***************

* gen ATC curve (OLS - Predicted)
use `data'/costs.dta, clear
gen costscatter=1
replace connections=connect_d if connections==.
gen atc_usd=cost_b1
gen connections2=connections^2

reg atc_usd connections connections2, vce(robust)
mat b=e(b)
matlist b
local a=b[1,3]
local b1=b[1,1]
local b2=b[1,2]

* note: we plot the population-weighted ATC curve (i.e., larger communities are assigned 
* higher weights) corresponding to the predicted cost of connecting various shares 
* of baseline-unconnected households for each community. 

egen totc=sum(compounds)
gen wa=compounds/totc
drop totc

forval n = 1/100 {
	gen x_`n'=`n'/100*ucompounds
	gen co_`n'=`a'+`b1'*x_`n'+`b2'*x_`n'^2	
	gen w_co_`n'=co_`n'*wa
	egen cuco_`n'=sum(w_co_`n')
	drop x_`n' co_`n' w_co_`n'
}

duplicates drop cuco_1, force
keep cuco_*
gen id=1

reshape long cuco_, i(id)
drop id
rename _j proportion
rename cuco_ cost_b1
keep cost_b1 proportion
gen ols=1
save `temp'/temp.dta, replace

* gen ATC curve (NL - Predicted & Weighted)
use `data'/costs.dta, clear
gen costscatter=1
replace connections=connect_d if connections==.

* estimate nl coefficients using quantities
nl (cost_b1={b0}/connections+{b1}+{b2}*connections)
mat b=e(b)
matlist b
/*
             | b0        | b1        | b2        
             |     _cons |     _cons |     _cons 
-------------+-----------+-----------+-----------
          y1 |  2453.441 |  999.4473 |  -3.24163 
*/
local b0=b[1,1]
local b1=b[1,2]
local b2=b[1,3]

* note: we plot the population-weighted ATC curve (i.e., larger communities are assigned 
* higher weights) corresponding to the predicted cost of connecting various shares 
* of baseline-unconnected households for each community. 

egen totc=sum(compounds)
gen wa=compounds/totc
drop totc

* predict ATC of achieving different proportions connected, weighting each community by total compounds
forval n = 1/100 {
	gen x_`n'=`n'/100*ucompounds
	gen co_`n'=`b0'/x_`n'+`b1'+`b2'*x_`n'
	gen mc_`n'=`b1'+2*`b2'*x_`n'	
	gen w_co_`n'=co_`n'*wa
	gen w_mc_`n'=mc_`n'*wa
	egen cuco_`n'=sum(w_co_`n')
	egen cumc_`n'=sum(w_mc_`n')	
	drop x_`n' co_`n' mc_`n' w_co_`n' w_mc_`n'
}

duplicates drop cuco_1, force
keep cuco_* cumc_*
gen id=1

reshape long cuco_ cumc_, i(id)
drop id
rename _j proportion
rename cuco_ cost_b1
rename cumc_ mc_b1
preserve
keep cost_b1 proportion
gen atc=1
save `temp'/atc.dta, replace
restore
keep mc_b1 proportion
gen mc=1
append using `temp'/atc.dta
sort proportion

append using `temp'/temp.dta
append using `data'/costs.dta
gen costscatter=1 if designed!=.

keep proportion cost_b1 mc_b1 atc mc ols costscatter designed
save `temp'/cost_data.dta, replace

* Extract and save demand scatter plot data
use `data'/demand.dta, clear
drop if connected==1
gen temp=takeup*100
drop takeup
rename temp takeup
bysort siteno: egen proportion=mean(takeup)
duplicates drop siteno, force
keep price proportion
gen temp=price/`rate'
gen p=round(temp,1)
keep p proportion
rename p price
gen demandscatter=1
sort price
save `temp'/demand_data.dta, replace

* Extract and save demand curve
use `data'/demand.dta, clear
bysort price: egen takeup3=mean(takeup)
duplicates drop price, force
drop if price==.
keep price takeup3
gen proportion=takeup3*100
drop takeup3
rename price kprice
sort kprice
gen price=kprice/`rate'
keep price proportion
gen demandcurve=1
append using `temp'/demand_data.dta
save `temp'/demand_data.dta, replace

* Plot
use `temp'/cost_data.dta, clear
append using `temp'/demand_data.dta

replace price=cost_b1 if atc==1 | costscatter==1
replace price=mc_b1 if mc==1

twoway ///
	(scatter price proportion if demandscatter==1 & proportion <=100, msize(small) mlcolor(gs6) mfcolor(gs8) msymbol(O)) ///
	(scatter price proportion if costscatter==1 & price<=3000 & designed==0 & proportion <=100, msize(small) mlcolor(gs10) mfcolor(gs13) msymbol(O)) ///
	(scatter price proportion if costscatter==1 & price<=3000 & designed==1 & proportion <=100, msize(small) mlcolor(gs10) mfcolor(white) msymbol(O)) ///
	(line price proportion if atc==1 & price<=3000 & proportion <=100, lcolor(gs8) lwidth(thick) lpattern(solid)) ///
	(line price proportion if mc==1 & price<=3000 & proportion <=100, lcolor(gs4) lwidth(thin) lpattern(dash) msymbol(none)) ///
	(connected price proportion if demandcurve==1 & proportion <=100, lcolor(gs5) lwidth(thick) lpattern(solid) msize(medsmall) mcolor(gs2) msymbol(square)), ///
	xtitle("Take-up, community coverage (%)", size(medium) margin(medsmall)) ///
	ytitle("Connection price, ATC per connection (USD)", size(medium) margin(medsmall)) ///
	ylabel(0(500)3000, labsize(medsmall)) ///
	xlabel(0(20)100, labsize(medsmall)) ///
	xsize(9) ysize(11) ///
	plotr(margin(0.3 0.3 0.4 3)) ///
	legend(size(medsmall) label(1 "Community takeup") label(2 "Community ATC (actual)") ///
	label(3 "Community ATC (designedd)") label(4 "ATC curve") label(5 "MC curve") ///
	label(6 "Demand curve") order(6 4 5) symx(9) region(lstyle(none)) position(1) ring(0) rows(5)) ///
	scheme(s1mono)
	graph export `figDir'/figure3a.png, replace

***************
*** Panel B ***
***************

sum cost_b1 if proportion==100 & atc==1, detail
* ATC at 100% is $ 739.29
local line=`r(mean)'
gen masscost=`line'
di `line'*84.7

* [0, 100] area under the lowest cost is $62,618 at 84.7 average community compound density.

* Extend demand curve through [0, 1.3] with regular projection
* Slope of first section: -19.602752
* b = 397.9986-(-19.602752*1.304348) = 423.56741 (projected intercept)

display (95.52631-23.7467)*170.5708*0.5 + ///
	(23.7467-7.105263)*(170.5708) + ///
	(23.7467-7.105263)*(284.2847-170.5708)*0.5 + ///
	(7.105263-1.304348)*(397.9986-284.2847)*0.5 + ///
	(7.105263-1.304348)*(284.2847) + ///
	(1.304348)*(423.56741-397.9986)*0.5 + ///
	(1.304348)*(397.9986)
			
* Measured [1.3, 100] area under the demand curve is $11885.411
* Implied [0, 100] area is 12421.215
* for now extend to zero

gen shadedcost=demandcurve
set obs 534
replace proportion=100 in 532
replace masscost=`line' in 532
replace shadedcost=1 in 532
replace proportion=0 in 533
replace masscost=`line' in 533
replace shadedcost=1 in 533

gen shadeddemand=demandcurve
replace proportion=0 in 534
replace price=423.56741 in 534
replace shadeddemand=1 in 534

sort proportion shadedcost shadeddemand

preserve
keep if price<=3000 & proportion <=100
keep if atc==1 | mc==1 
keep price proportion atc mc
sort atc proportion
save `temp'/atc_curve.dta, replace
restore

twoway ///
	(area masscost proportion if shadedcost==1, lpattern(solid) lstyle(unextended) lwidth(vthin) lcolor(gs13) fcolor(gs13)) ///
	(area price proportion if shadeddemand==1, lpattern(solid)  lstyle(unextended) lwidth(vthin) lcolor(gs11) fcolor(gs11)) ///
	(line price proportion if atc==1 & price<=3000 & proportion <=100, lcolor(gs8) lwidth(thick) lpattern(solid)) ///
	(line price proportion if demandcurve==1 & proportion <=100, lcolor(gs5) lwidth(thick) lpattern(solid)), ///
	xtitle("Take-up, community coverage (%)", size(medium) margin(medsmall)) ///
	ytitle("Connection price, ATC per connection (USD)", size(medium) margin(medsmall)) ///
	ylabel(0(500)3000, labsize(medsmall)) ///
	xlabel(0(20)100, labsize(medsmall)) ///
	xsize(9) ysize(11) ///
	plotr(margin(0.3 0.3 0.4 3)) ///
	text(650 70 "TC = $62,618", just(left) place(e) size(medsmall)) ///
	text(80 2 "CS = $12,421", just(left) place(e) size(medsmall)) ///
	legend(size(medsmall) label(1 "Area 1") label(2 "Area 2") label(3 "ATC curve") ///
	label(4 "Demand curve") order(4 3) symx(9) region(lstyle(none)) position(1) ring(0) rows(5)) ///
	scheme(s1mono)
	graph export `figDir'/figure3b.png, replace

display (62617.973-12421.215)
display (62617.973-12421.215)/84.7

* $50,197 gap between cost and demand
* $593 per household

***************
*** Panel C ***
***************

local low=15
local high=25

* CV question 1 (16b): "Would you be willing to pay XX Ksh for an electricity connection?"
use `data'/demand.dta, clear
bysort WTP_amt: egen takeup1=mean(WTP_r1)
drop if WTP_amt==94 | WTP_amt==95 | WTP_amt==5000 | WTP_amt==30000
* not enough observations (we asked about a couple different data points in the first week of surveying (randomly) before settling on set of price points)

tab WTP_amt
/*
   Randomly |
   assigned |
   price of |
       grid |
 connection |
     for CV |
   question |      Freq.     Percent        Cum.
------------+-----------------------------------
          0 |        202        8.86        8.86
      10000 |        375       16.45       25.31
      15000 |        382       16.75       42.06
      20000 |        369       16.18       58.25
      25000 |        396       17.37       75.61
      35000 |        370       16.23       91.84
      75000 |        186        8.16      100.00
------------+-----------------------------------
      Total |      2,280      100.00
*/
duplicates drop WTP_amt, force
keep WTP_amt takeup1
gen cv1=takeup1*100
drop takeup1
rename WTP_amt kprice
sort kprice
save `temp'/cv1.dta, replace

* CV question 2 (16c): "16c. Imagine that you were offered an electricity connection at this price today, and you were given 6 weeks to complete the payment. Would you accept the offer?"
use `data'/demand.dta, clear
bysort WTP_amt: egen takeup2=mean(WTP_r2)
drop if WTP_amt==94 | WTP_amt==95 | WTP_amt==5000 | WTP_amt==30000
* not enough observations (we asked about a couple different data points in the first week of surveying (randomly) before settling on set of price points)
duplicates drop WTP_amt, force
keep WTP_amt takeup2
gen cv2=takeup2*100
drop takeup2
rename WTP_amt kprice
sort kprice
save `temp'/cv2.dta, replace

* Experimental results
use `data'/demand.dta, clear
bysort price: egen takeup3=mean(takeup)
duplicates drop price, force
drop if price==.
keep price takeup3
gen rp=takeup3*100
drop takeup3
rename price kprice
sort kprice
save `temp'/rp.dta, replace

* Financing WTP results - Low discount rate, CV question 1 (16f)
use `data'/demand.dta, clear
drop if fin_group==.
bysort fin_group: egen temp=mean(fin_WTP_r1)
gen fincv1_`low'=temp*100
duplicates drop fin_group, force
keep fin_npv`low' fincv1_`low'
rename fin_npv`low' kprice
sort kprice
save `temp'/fin_cv1_`low'.dta, replace

* Financing WTP results - High discount rate, CV question 1 (16f)
use `data'/demand.dta, clear
drop if fin_group==.
bysort fin_group: egen temp=mean(fin_WTP_r1)
gen fincv1_`high'=temp*100
duplicates drop fin_group, force
keep fin_npv`high' fincv1_`high'
rename fin_npv`high' kprice
sort kprice
save `temp'/fin_cv1_`high'.dta, replace

* Financing WTP results - Low discount rate, CV question 2 (16g)
use `data'/demand.dta, clear
drop if fin_group==.
bysort fin_group: egen temp=mean(fin_WTP_r2)
gen fincv2_`low'=temp*100
duplicates drop fin_group, force
keep fin_npv`low' fincv2_`low'
rename fin_npv`low' kprice
sort kprice
save `temp'/fin_cv2_`low'.dta, replace

* Financing WTP results - High discount rate, CV question 2 (16g)
use `data'/demand.dta, clear
drop if fin_group==.
bysort fin_group: egen temp=mean(fin_WTP_r2)
gen fincv2_`high'=temp*100
duplicates drop fin_group, force
keep fin_npv`high' fincv2_`high'
rename fin_npv`high' kprice
sort kprice
save `temp'/fin_cv2_`high'.dta, replace

* Merge files
use `temp'/cv1.dta, clear
merge kprice using `temp'/cv2.dta
drop _merge 
sort kprice
merge kprice using `temp'/rp.dta
drop _merge
sort kprice
merge kprice using `temp'/fin_cv1_`low'.dta
drop _merge
sort kprice
merge kprice using `temp'/fin_cv1_`high'.dta
drop _merge
sort kprice
merge kprice using `temp'/fin_cv2_`low'.dta
drop _merge
sort kprice
merge kprice using `temp'/fin_cv2_`high'.dta
drop _merge
gen price=kprice/`rate'
gen experimentpoint=1
sort kprice
save `temp'/experimental_takeup.dta, replace
rm `temp'/cv1.dta
rm `temp'/cv2.dta
rm `temp'/rp.dta
rm `temp'/fin_cv1_`low'.dta
rm `temp'/fin_cv1_`high'.dta
rm `temp'/fin_cv2_`low'.dta
rm `temp'/fin_cv2_`high'.dta

* We assume that take-up rates for CV and CV (with financing) are the same for price of 0 KES 
replace fincv1_`low'=cv1 if kprice==0
replace fincv1_`high'=cv1 if kprice==0
replace fincv2_`low'=cv1 if kprice==0
replace fincv2_`high'=cv1 if kprice==0

append using `temp'/atc_curve.dta

twoway ///
		(line price proportion if proportion!=. & atc==1, lcolor(gs8) lwidth(thick) lpattern(solid)) ///
		(line price proportion if proportion!=. & mc==1, lcolor(gs4) lwidth(thin) lpattern(dash) msymbol(none)) ///
		(connected price fincv1_`low', lcolor(gs2) lwidth(med) lpattern(shortdash) msize(medsmall) mcolor(gs1) msymbol(circle)) ///
		(connected price fincv2_`low', lcolor(gs10) lwidth(med) lpattern(shortdash) msize(medsmall) mcolor(gs8) msymbol(circle)) ///
		(connected price cv1, lcolor(gs2) lwidth(medthin) lpattern(longdash) msize(small) mcolor(gs1) msymbol(square)) ///
		(connected price cv2, lcolor(gs10) lwidth(medthin) lpattern(longdash) msize(small) mcolor(gs8) msymbol(square)) ///		
		(connected price rp if experimentpoint==1, lcolor(gs5) lwidth(thick) lpattern(solid) msize(small) mcolor(gs2) msymbol(square)), ///
		xtitle("Take-up, community coverage (%)", size(medium) margin(medsmall)) ///
		ytitle("Connection price, ATC per connection (USD)", size(medium) margin(medsmall)) ///
		ylabel(0(500)3000, labsize(medsmall)) ///
		xlabel(0(20)100, labsize(medsmall)) ///
		xsize(9) ysize(11) ///
		scheme(s1mono) ///
		legend(size(medsmall) ///
		label(1 "ATC curve") ///
		label(2 "MC curve") ///
		label(3 "Credit") ///
		label(4 "Credit (time-limited offer)") ///
		label(5 "WTP") ///
		label(6 "WTP (time-limited offer)") ///	
		label(7 "Demand curve") ///	
		symx(9) order(7 1 2 5 6 3 4) region(lstyle(none) fcolor(none) margin(0 0 0 1)) position(1) ring(0) rows(7))
		graph export `figDir'/figure3c.png, replace

rm `temp'/atc_curve.dta
rm `temp'/atc.dta
rm `temp'/temp.dta
rm `temp'/cost_data.dta
rm `temp'/demand_data.dta
rm `temp'/experimental_takeup.dta




