clear all
set more off
local rootDir = "/Users/Research/LMW2019"

local data "`rootDir'/data/stata"
local temp "`rootDir'/data/temp"
local appDir "`rootDir'/appendix"

local rate=87.94

********************************************************************************
* Appendix Figures
* Last updated: 23 Apr 2019
********************************************************************************

********************************************************************************
* Figure B7: Stated reasons why households remain unconnected to electricity
********************************************************************************

use `data'/demand.dta, clear
keep if connected==0
keep reasons
rename reasons f2

split f2
destring f21 f22 f23 f24, replace

gen reason_unavailable = (f21==1|f22==1|f23==1|f24==1)
gen reason_aff_connect = (f21==2|f22==2|f23==2|f24==2)
gen reason_aff_wiring = (f21==3|f22==3|f23==3|f24==3)
gen reason_aff_bill = (f21==4|f22==4|f23==4|f24==4)
gen reason_aff_appliance = (f21==5|f22==5|f23==5|f24==5)
gen reason_satisfied = (f21==6|f22==6|f23==6|f24==6)
gen reason_permission = (f21==7|f22==7|f23==7|f24==7)
gen reason_theft = (f21==8|f22==8|f23==8|f24==8)
gen reason_app_hassle = (f21==9|f22==9|f23==9|f24==9)
gen reason_useneighbors = (f21==10|f22==10|f23==10|f24==10)
gen reason_nobenefit = (f21==11|f22==11|f23==11|f24==11)
gen reason_notsafe = (f21==12|f22==12|f23==12|f24==12)
gen reason_unreliable = (f21==13|f22==13|f23==13|f24==13)
drop f2 f21 f22 f23 f24

keep reason_*

local reasonList unavailable aff_connect aff_wiring aff_bill aff_appliance satisfied permission ///
	app_hassle useneighbors nobenefit notsafe unreliable
local i=1
foreach r of local reasonList {
	egen temp=mean(reason_`r')
	gen v`i'=temp*100
	drop temp
	local i=`i'+1
}

duplicates drop v1, force
drop reason_*
gen id=1

reshape long v, i(id) j(variable)
drop id
rename v proportion
gsort -proportion
gen n=_n

egen other=sum(proportion) if n>=6
replace proportion=other if n==6
drop if n>=7
drop variable other

label define reasons 1 "Connection cost" 2 "Wiring cost" 3 "Monthly cost" ///
	4 "Availability" 5 "Hassle" 6 "Other"
label values n reasons

gen prop=round(proportion,0.1)

twoway	(bar proportion n, barw(0.5) lcolor(gs5) lwidth(thin) color(gs11)) ///
		(scatter proportion n, mlabel(prop) mlabpos(12) ms(none) mlabsize(small) legend(off)), ///
		ytitle("Proportion of households (%)", size(medsmall) margin(medsmall)) ///
		xtitle("", size(medsmall) margin(medsmall)) ///
		ylabel(0(10)100, valuelabel labsize(small)) ///
		scheme(s1mono) plotregion(margin(b = 0)) ///
		xlabel(1 "Connection cost" 2 "Wiring cost" 3 "Monthly cost" ///
		4 "Availability" 5 "Hassle" 6 "Other reason", labsize(small)) ///
		xsize(15) ysize(10)
		graph export `appDir'/figb7.png, replace

********************************************************************************
* Figure B8: Demand for rural electrification
********************************************************************************

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
		(connected price pap, lcolor(gs2) lwidth(medthin) lpattern(longdash) msize(small) mcolor(gs1) msymbol(square)) ///
		(connected price rp, lcolor(gs5) lwidth(thick) lpattern(solid) msize(small) mcolor(gs2) msymbol(square)), ///
		xtitle("Take-up (%)", size(medium) margin(medsmall)) ///
		ytitle("Connection price (USD)", size(medium) margin(medsmall)) ///
		xlabel(0(20)100, labsize(medsmall)) ///
		ylabel(0(50)400, labsize(medsmall)) ///
		xsize(9) ysize(11) ///
		scheme(s1mono) ///
		legend(size(medsmall) ///
		label(1 "Kenyan gov't report") ///
		label(2 "Pre-analysis plan") ///
		label(3 "Experiment") ///
		symx(9) order(3 1 2) region(lstyle(none)) position(1) ring(0) rows(3))
		graph export `appDir'/figb8.png, replace

rm `temp'/rp.dta

********************************************************************************
* Figure B9A: Costs of rural electrification
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

twoway ///
	(scatter cost_b1 proportion if designed==0 & costscatter==1, msize(medsmall) mlcolor(gs10) mfcolor(gs12) msymbol(O)) ///
	(scatter cost_b1 proportion if designed==1 & costscatter==1, msize(medsmall) mlcolor(gs10) mfcolor(white) msymbol(O)) ///
	(line cost_b1 proportion if ols==1, sort lcolor(gs11) lwidth(med) lpattern(solid)) ///
	(line cost_b1 proportion if atc==1, sort lcolor(gs8) lwidth(thick) lpattern(solid)), ///
	ytitle("ATC per connection (USD)", size(medium) margin(medsmall)) ///
	xtitle("Community coverage (%)", size(medium) margin(medsmall)) ///
	ylabel(0(1000)6000, labsize(medsmall)) ///
	xlabel(0(20)100, labsize(medsmall)) ///
	xsize(9) ysize(11) ///
	legend(size(medsmall) label(1 "Sample communities") label(2 "Designed communities") ///
	label(4 "ATC curve (Nonlinear)") label(3 "ATC curve (Quadratic)") ///
	order(3 4 1 2) symx(9) region(lstyle(none)) position(1) ring(0) rows(4)) ///
	scheme(s1mono)
	graph export `appDir'/figb9aa.png, replace

rm `temp'/temp.dta
rm `temp'/atc.dta

***************
*** Panel B ***
***************

use `data'/costs.dta, clear
gen costscatter=1
drop if connections==.
nl (cost_b1={b0}/connections+{b1}+{b2}*connections)
mat b=e(b)
matlist b
local b0=b[1,1]
local b1=b[1,2]
local b2=b[1,3]

egen totc=sum(compounds)
gen wa=compounds/totc
drop totc

forval n = 1/38 {
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
keep cost_b1 proportion
gen atc=1
gen s=1
save `temp'/atc_s.dta, replace

use `data'/costs.dta, clear
gen costscatter=1
replace connections=connect_d if connections==.
nl (cost_b1={b0}/connections+{b1}+{b2}*connections)
mat b=e(b)
matlist b
local b0=b[1,1]
local b1=b[1,2]
local b2=b[1,3]

egen totc=sum(compounds)
gen wa=compounds/totc
drop totc

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
keep cost_b1 proportion
gen atc=1
gen sd=1
save `temp'/atc_sd.dta, replace
append using `temp'/atc_s.dta
sort proportion
append using `data'/costs.dta
gen costscatter=1 if designed!=.

keep proportion cost_b1 atc s sd costscatter designed
save `temp'/cost_data.dta, replace

twoway ///
	(scatter cost_b1 proportion if designed==0 & costscatter==1, msize(medsmall) mlcolor(gs10) mfcolor(gs12) msymbol(O)) ///
	(scatter cost_b1 proportion if designed==1 & costscatter==1, msize(medsmall) mlcolor(gs10) mfcolor(white) msymbol(O)) ///
	(line cost_b1 proportion if atc==1 & s==1, sort lcolor(gs11) lwidth(med) lpattern(solid)) ///
	(line cost_b1 proportion if atc==1 & sd==1, sort lcolor(gs8) lwidth(thick) lpattern(solid)), ///
	ytitle("ATC per connection (USD)", size(medium) margin(medsmall)) ///
	xtitle("Community coverage (%)", size(medium) margin(medsmall)) ///
	ylabel(0(1000)6000, labsize(medsmall)) ///
	xlabel(0(20)100, labsize(medsmall)) ///
	xsize(9) ysize(11) ///
	legend(size(medsmall) label(1 "Sample communities") label(2 "Designed communities") ///
	label(3 "ATC curve (Sample)") label(4 "ATC curve (Sample & Designed)") ///
	order(3 4 1 2) symx(9) region(lstyle(none)) position(1) ring(0) rows(4)) ///
	scheme(s1mono)
	graph export `appDir'/figb9ab.png, replace	

rm `temp'/atc_s.dta
rm `temp'/atc_sd.dta

********************************************************************************
* Figure B9B: Natural monopoly: Alternative functional forms
********************************************************************************

***************
*** Panel A ***
***************

* Note: This is Figure 4, Panel A

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
	graph export `appDir'/figb9ba.png, replace

***************
*** Panel B ***
***************

use `data'/costs.dta, clear
gen costscatter=1
replace connections=connect_d if connections==.

nl (cost_b1={b0}+{b1}/connections)
mat b=e(b)
matlist b
local b0=b[1,1]
local b1=b[1,2]

egen totc=sum(compounds)
gen wa=compounds/totc
drop totc

forval n = 1/100 {
	gen x_`n'=`n'/100*ucompounds
	gen co_`n'=`b0'+`b1'/x_`n'
	gen mc_`n'=`b0'
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
append using `data'/costs.dta
gen costscatter=1 if designed!=.
keep proportion cost_b1 mc_b1 atc mc costscatter designed
save `temp'/cost_data.dta, replace

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
	xlabel(0(10)100, labsize(medsmall)) ///
	xsize(9) ysize(11) ///
	plotr(margin(0.3 0.3 0.4 3)) ///
	legend(size(medsmall) label(1 "Community takeup") label(2 "Community ATC (actual)") ///
	label(3 "Community ATC (designedd)") label(4 "ATC curve") label(5 "MC curve") ///
	label(6 "Demand curve") order(6 4 5) symx(9) region(lstyle(none)) position(1) ring(0) rows(5)) ///
	scheme(s1mono)
	graph export `appDir'/figb9bb.png, replace

***************
*** Panel C ***
***************

use `data'/costs.dta, clear
gen costscatter=1
replace connections=connect_d if connections==.
lowess cost_b1 connections, gen(newvar) nog
nl exp3: newvar connection
mat b=e(b)
matlist b
local b0=b[1,1]
local b1=b[1,2]
local b2=b[1,3]

egen totc=sum(compounds)
gen wa=compounds/totc
drop totc

forval n = 1/100 {
	gen x_`n'=`n'/100*ucompounds
	gen co_`n'=`b0'+`b1'*`b2'^(x_`n')
	gen mc_`n'=`b0'+`b1'*`b2'^(x_`n')+(x_`n')*ln(`b2')*`b1'*`b2'^(x_`n')
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
append using `data'/costs.dta
gen costscatter=1 if designed!=.
keep proportion cost_b1 mc_b1 atc mc costscatter designed
save `temp'/cost_data.dta, replace

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

use `temp'/cost_data.dta, clear
append using `temp'/demand_data.dta
replace price=cost_b1 if atc==1 | costscatter==1
replace price=mc_b1 if mc==1

twoway ///
	(scatter price proportion if demandscatter==1 & proportion <=100, msize(small) mlcolor(gs6) mfcolor(gs8) msymbol(O)) ///
	(scatter price proportion if costscatter==1 & price<=3000 & designed==0 & proportion <=100, msize(small) mlcolor(gs10) mfcolor(gs13) msymbol(O)) ///
	(scatter price proportion if costscatter==1 & price<=3000 & designed==1 & proportion <=100, msize(small) mlcolor(gs10) mfcolor(white) msymbol(O)) ///
	(line price proportion if atc==1 & price<=3000  & proportion <=100, lcolor(gs8) lwidth(thick) lpattern(solid)) ///
	(line price proportion if mc==1 & price<=3000  & proportion <=100, lcolor(gs4) lwidth(thin) lpattern(dash) msymbol(none)) ///
	(connected price proportion if demandcurve==1 & proportion <=100, lcolor(gs5) lwidth(thick) lpattern(solid) msize(medsmall) mcolor(gs2) msymbol(square)), ///
	xtitle("Take-up, community coverage (%)", size(medium) margin(medsmall)) ///
	ytitle("Connection price, ATC per connection (USD)", size(medium) margin(medsmall)) ///
	ylabel(0(500)3000, labsize(medsmall)) ///
	xlabel(0(10)100, labsize(medsmall)) ///
	xsize(9) ysize(11) ///
	plotr(margin(0.3 0.3 0.4 3)) ///
	legend(size(medsmall) label(1 "Community takeup") label(2 "Community ATC (actual)") ///
	label(3 "Community ATC (designed)") label(4 "ATC curve") label(5 "MC curve") ///
	label(6 "Demand curve") order(6 4 5) symx(9) region(lstyle(none)) position(1) ring(0) rows(5)) ///
	scheme(s1mono)
	graph export `appDir'/figb9bc.png, replace

********************************************************************************
* Figure B9C: Confidence intervals
********************************************************************************

use `data'/demand.dta, clear
drop if connected==1
gen temp=price/`rate'
gen p=round(temp,1)
keep p siteno takeup

gen temp=takeup*100
drop takeup
rename temp takeup
bysort siteno: egen proportion=mean(takeup)
duplicates drop siteno, force
keep p proportion
sort p

* CI for community-level means
gen lb=.
gen ub=.
foreach p in 0 171 284 398 {
	ci means proportion if p==`p'
	local l=`r(lb)'
	local u=`r(ub)'
	replace lb=`l' if p==`p'
	replace ub=`u' if p==`p'
}
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

* confidence intervals for nl curve
use `data'/costs.dta, clear
gen costscatter=1
replace connections=connect_d if connections==.

nl (cost_b1={b0}/connections+{b1}+{b2}*connections)
mat b=e(b)
matlist b
local b0=b[1,1]
local b1=b[1,2]
local b2=b[1,3]

predictnl fit = _b[/b0]/connections+_b[/b1]+_b[/b2]*connections, ci(lb ub) level(95)

nl (lb={b0}/connections+{b1}+{b2}*connections)
mat b=e(b)
matlist b
local l0=b[1,1]
local l1=b[1,2]
local l2=b[1,3]

nl (ub={b0}/connections+{b1}+{b2}*connections)
mat b=e(b)
matlist b
local u0=b[1,1]
local u1=b[1,2]
local u2=b[1,3]

* prep for constructing weighted average predicted nl curve (by compounds)
egen totc=sum(compounds)
gen wa=compounds/totc
drop totc

forval n = 1/100 {
	gen x_`n'=`n'/100*ucompounds
	gen co_`n'=`b0'/x_`n'+`b1'+`b2'*x_`n'
	gen w_co_`n'=co_`n'*wa
	egen cuco_`n'=sum(w_co_`n')
	drop co_`n' w_co_`n'

	gen co_`n'_u=`u0'/x_`n'+`u1'+`u2'*x_`n'
	gen w_co_`n'_u=co_`n'_u*wa
	egen cucou_`n'=sum(w_co_`n'_u)
	drop co_`n'_u w_co_`n'_u
	
	gen co_`n'_l=`l0'/x_`n'+`l1'+`l2'*x_`n'
	gen w_co_`n'_l=co_`n'_l*wa
	egen cucol_`n'=sum(w_co_`n'_l)
	drop co_`n'_l w_co_`n'_l
}

duplicates drop cuco_1, force
keep cuco_* cucol_* cucou_*
gen id=1

reshape long cuco_ cucou_ cucol_, i(id)
drop id
rename _j proportion

rename cuco_ cost_b1

preserve
keep cost_b1 proportion
gen atc=1
save `temp'/atc.dta, replace
restore
preserve
keep cucou_ proportion
rename cucou_ cost_b1
gen atc_u=1
save `temp'/atc_u.dta, replace
restore
keep cucol_ proportion
rename cucol_ cost_b1
gen atc_l=1
append using `temp'/atc_u.dta
append using `temp'/atc.dta
append using `data'/costs.dta
gen costscatter=1 if designed!=.
keep proportion cost_b1 atc atc_l atc_u designed costscatter
append using `temp'/demand_data.dta
replace price=cost_b1 if price==.

replace price=2499.139 if atc_u==1 & proportion==2

twoway ///
	(area price proportion if atc_u==1 & price<=3000, fcolor(gs12) lpattern(solid) lcolor(gs12) lwidth(none)) ///
 	(area price proportion if atc_l==1 & price<=3000, fcolor(white) lpattern(solid) lcolor(white) lwidth(none)) ///
	(area ub price if demandscatter==1 & proportion<=100, fcolor(gs12) lpattern(solid) lcolor(gs12) lwidth(none) hor) ///
 	(area lb price if demandscatter==1 & proportion<=100, fcolor(white) lpattern(solid) lcolor(white) lwidth(none) hor) ///
	(scatter price proportion if demandscatter==1 & proportion <=100, msize(small) mlcolor(gs6) mfcolor(gs8) msymbol(O)) ///
	(scatter price proportion if costscatter==1 & price<=3000 & designed==0 & proportion <=100, msize(small) mlcolor(gs10) mfcolor(gs13) msymbol(O)) ///
	(scatter price proportion if costscatter==1 & price<=3000 & designed==1 & proportion <=100, msize(small) mlcolor(gs10) mfcolor(white) msymbol(O)) ///
	(line price proportion if atc==1 & price<=3000 & proportion <=100, lcolor(gs8) lwidth(thick) lpattern(solid)) ///
	(line price proportion if demandcurve==1 & proportion <=100, lcolor(gs5) lwidth(thick) lpattern(solid) msize(medsmall) mcolor(gs2) msymbol(square)), ///
	xtitle("Take-up, community coverage (%)", size(medium) margin(medsmall)) ///
	ytitle("Connection price, ATC per connection (USD)", size(medium) margin(medsmall)) ///
	ylabel(0(500)3000, labsize(medsmall)) ///
	xlabel(0(20)100, labsize(medsmall)) ///
	xsize(9) ysize(11) ///
	plotr(margin(0.3 0.3 0.4 3)) ///
	legend(size(medsmall) label(1 "95% CI") label(8 "ATC curve") ///
	label(9 "Demand curve") label(6 "Sample communities") label(7 "Designed communities") ///
	order(9 8 1 6 7) symx(9) region(lstyle(none)) position(1) ring(0) rows(5)) ///
	scheme(s1mono)
	graph export `appDir'/figb9c.png, replace

rm `temp'/cost_data.dta
rm `temp'/demand_data.dta
rm `temp'/atc.dta
rm `temp'/atc_u.dta
rm `temp'/temp.dta

********************************************************************************
* Figure B10: How ATC per connection varies along different dimensions
********************************************************************************

***************
*** Panel A ***
***************

local var "population"

use `data'/costs.dta, clear
gen costscatter=1
replace connections=connect_d if connections==.

sum `var', detail
local median = `r(p50)'

nl (cost_b1={b0}/connections+{b1}+{b2}*connections) if `var'<`median'
mat b=e(b)
matlist b
local b0=b[1,1]
local b1=b[1,2]
local b2=b[1,3]

nl (cost_b1={b0}/connections+{b1}+{b2}*connections) if `var'>`median'
mat b=e(b)
matlist b
local c0=b[1,1]
local c1=b[1,2]
local c2=b[1,3]

egen totc=sum(compounds)
gen wa=compounds/totc
drop totc

forval n = 1/100 {
	gen x_`n'=`n'/100*ucompounds
	gen co_1_`n'=`b0'/x_`n'+`b1'+`b2'*x_`n'
	gen co_2_`n'=`c0'/x_`n'+`c1'+`c2'*x_`n'
	gen w_co_1_`n'=co_1_`n'*wa
	gen w_co_2_`n'=co_2_`n'*wa
	egen cuco_1_`n'=sum(w_co_1_`n')
	egen cuco_2_`n'=sum(w_co_2_`n')
	drop x_`n' co_*_`n' w_co_*_`n'
}

duplicates drop cuco_1_1, force
keep cuco_*
gen id=1

reshape long cuco_1_ cuco_2_, i(id)
drop id
rename _j proportion
rename cuco_1_ cost_b1_1
rename cuco_2_ cost_b1_2
keep cost_b1_1 cost_b1_2 proportion
gen atc=1
preserve
keep cost_b1_1 proportion atc
rename cost_b1_1 cost_b1
gen curve=0
save `temp'/temp.dta, replace
restore
keep cost_b1_2 proportion atc
rename cost_b1_2 cost_b1
gen curve=1
append using `temp'/temp.dta
append using `data'/costs.dta
gen costscatter=1 if designed!=.
keep proportion cost_b1 atc costscatter designed curve `var'

twoway ///
	(scatter cost_b1 proportion if costscatter==1 & `var'<`median', msize(medsmall) mlcolor(gs4) mfcolor(gs7) msymbol(O)) ///
	(scatter cost_b1 proportion if costscatter==1 & `var'>`median', msize(medsmall) mlcolor(gs8) mfcolor(gs11) msymbol(O)) ///
	(line cost_b1 proportion if atc==1 & curve==0, sort lcolor(gs7) lwidth(thick) lpattern(solid)) ///
	(line cost_b1 proportion if atc==1 & curve==1, sort lcolor(gs11) lwidth(thick) lpattern(solid)), ///
	ytitle("ATC per connection (USD)", size(medium) margin(medsmall)) ///
	xtitle("Community coverage (%)", size(medium) margin(medsmall)) ///
	ylabel(0(1000)6000, labsize(medsmall)) ///
	xlabel(0(10)100, labsize(medsmall)) ///
	xsize(9) ysize(11) ///
	legend(size(medsmall) label(1 "Low population") label(2 "High population") ///
	label(3 "Low population ATC curve") label(4 "High population ATC curve") ///
	order(3 4 1 2) symx(9) region(lstyle(none)) position(1) ring(0) rows(4)) ///
	scheme(s1mono)
	graph export `appDir'/figb10a.png, replace

***************
*** Panel B ***
***************

local var "terrain_2"

use `data'/costs.dta, clear
gen costscatter=1
replace connections=connect_d if connections==.

sum `var', detail
local median = `r(p50)'

nl (cost_b1={b0}/connections+{b1}+{b2}*connections) if `var'<`median'
mat b=e(b)
matlist b
local b0=b[1,1]
local b1=b[1,2]
local b2=b[1,3]

nl (cost_b1={b0}/connections+{b1}+{b2}*connections) if `var'>`median'
mat b=e(b)
matlist b
local c0=b[1,1]
local c1=b[1,2]
local c2=b[1,3]

egen totc=sum(compounds)
gen wa=compounds/totc
drop totc

forval n = 1/100 {
	gen x_`n'=`n'/100*ucompounds
	gen co_1_`n'=`b0'/x_`n'+`b1'+`b2'*x_`n'
	gen co_2_`n'=`c0'/x_`n'+`c1'+`c2'*x_`n'
	gen w_co_1_`n'=co_1_`n'*wa
	gen w_co_2_`n'=co_2_`n'*wa
	egen cuco_1_`n'=sum(w_co_1_`n')
	egen cuco_2_`n'=sum(w_co_2_`n')
	drop x_`n' co_*_`n' w_co_*_`n'
}

duplicates drop cuco_1_1, force
keep cuco_*
gen id=1

reshape long cuco_1_ cuco_2_, i(id)
drop id
rename _j proportion
rename cuco_1_ cost_b1_1
rename cuco_2_ cost_b1_2
keep cost_b1_1 cost_b1_2 proportion
gen atc=1
preserve
keep cost_b1_1 proportion atc
rename cost_b1_1 cost_b1
gen curve=0
save `temp'/temp.dta, replace
restore
keep cost_b1_2 proportion atc
rename cost_b1_2 cost_b1
gen curve=1
append using `temp'/temp.dta
append using `data'/costs.dta
gen costscatter=1 if designed!=.
keep proportion cost_b1 atc costscatter designed curve `var'

twoway ///
	(scatter cost_b1 proportion if costscatter==1 & `var'<`median', msize(medsmall) mlcolor(gs4) mfcolor(gs7) msymbol(O)) ///
	(scatter cost_b1 proportion if costscatter==1 & `var'>`median', msize(medsmall) mlcolor(gs8) mfcolor(gs11) msymbol(O)) ///
	(line cost_b1 proportion if atc==1 & curve==0, sort lcolor(gs7) lwidth(thick) lpattern(solid)) ///
	(line cost_b1 proportion if atc==1 & curve==1, sort lcolor(gs11) lwidth(thick) lpattern(solid)), ///
	ytitle("ATC per connection (USD)", size(medium) margin(medsmall)) ///
	xtitle("Community coverage (%)", size(medium) margin(medsmall)) ///
	ylabel(0(1000)6000, labsize(medsmall)) ///
	xlabel(0(10)100, labsize(medsmall)) ///
	xsize(9) ysize(11) ///
	legend(size(medsmall) label(1 "Low gradient") label(2 "High gradient") ///
	label(3 "Low gradient ATC curve") label(4 "High gradient ATC curve") ///
	order(3 4 1 2) symx(9) region(lstyle(none)) position(1) ring(0) rows(4)) ///
	scheme(s1mono)
	graph export `appDir'/figb10b.png, replace

***************
*** Panel C ***
***************
	
local var "electrification"

use `data'/costs.dta, clear
gen costscatter=1
replace connections=connect_d if connections==.

sum `var', detail
local median = `r(p50)'

nl (cost_b1={b0}/connections+{b1}+{b2}*connections) if `var'<`median'
mat b=e(b)
matlist b
local b0=b[1,1]
local b1=b[1,2]
local b2=b[1,3]

nl (cost_b1={b0}/connections+{b1}+{b2}*connections) if `var'>`median'
mat b=e(b)
matlist b
local c0=b[1,1]
local c1=b[1,2]
local c2=b[1,3]

egen totc=sum(compounds)
gen wa=compounds/totc
drop totc

forval n = 1/100 {
	gen x_`n'=`n'/100*ucompounds
	gen co_1_`n'=`b0'/x_`n'+`b1'+`b2'*x_`n'
	gen co_2_`n'=`c0'/x_`n'+`c1'+`c2'*x_`n'
	gen w_co_1_`n'=co_1_`n'*wa
	gen w_co_2_`n'=co_2_`n'*wa
	egen cuco_1_`n'=sum(w_co_1_`n')
	egen cuco_2_`n'=sum(w_co_2_`n')
	drop x_`n' co_*_`n' w_co_*_`n'
}

duplicates drop cuco_1_1, force
keep cuco_*
gen id=1

reshape long cuco_1_ cuco_2_, i(id)
drop id
rename _j proportion
rename cuco_1_ cost_b1_1
rename cuco_2_ cost_b1_2
keep cost_b1_1 cost_b1_2 proportion
gen atc=1
preserve
keep cost_b1_1 proportion atc
rename cost_b1_1 cost_b1
gen curve=0
save `temp'/temp.dta, replace
restore
keep cost_b1_2 proportion atc
rename cost_b1_2 cost_b1
gen curve=1
append using `temp'/temp.dta
append using `data'/costs.dta
gen costscatter=1 if designed!=.
keep proportion cost_b1 atc costscatter designed curve `var'

twoway ///
	(scatter cost_b1 proportion if costscatter==1 & `var'<`median', msize(medsmall) mlcolor(gs4) mfcolor(gs7) msymbol(O)) ///
	(scatter cost_b1 proportion if costscatter==1 & `var'>`median', msize(medsmall) mlcolor(gs8) mfcolor(gs11) msymbol(O)) ///
	(line cost_b1 proportion if atc==1 & curve==0, sort lcolor(gs7) lwidth(thick) lpattern(solid)) ///
	(line cost_b1 proportion if atc==1 & curve==1, sort lcolor(gs11) lwidth(thick) lpattern(solid)), ///
	ytitle("ATC per connection (USD)", size(medium) margin(medsmall)) ///
	xtitle("Community coverage (%)", size(medium) margin(medsmall)) ///
	ylabel(0(1000)6000, labsize(medsmall)) ///
	xlabel(0(10)100, labsize(medsmall)) ///
	xsize(9) ysize(11) ///
	legend(size(medsmall) label(1 "Low electrification") label(2 "High electrification") ///
	label(3 "Low electrification ATC curve") label(4 "High electrification ATC curve") ///
	order(3 4 1 2) symx(9) region(lstyle(none)) position(1) ring(0) rows(4)) ///
	scheme(s1mono)
	graph export `appDir'/figb10c.png, replace

***************
*** Panel D ***
***************
	
local var "busia"

use `data'/costs.dta, clear
gen costscatter=1
replace connections=connect_d if connections==.

nl (cost_b1={b0}/connections+{b1}+{b2}*connections) if `var'==0
mat b=e(b)
matlist b
local b0=b[1,1]
local b1=b[1,2]
local b2=b[1,3]

nl (cost_b1={b0}/connections+{b1}+{b2}*connections) if `var'==1
mat b=e(b)
matlist b
local c0=b[1,1]
local c1=b[1,2]
local c2=b[1,3]

egen totc=sum(compounds)
gen wa=compounds/totc
drop totc

forval n = 1/100 {
	gen x_`n'=`n'/100*ucompounds
	gen co_1_`n'=`b0'/x_`n'+`b1'+`b2'*x_`n'
	gen co_2_`n'=`c0'/x_`n'+`c1'+`c2'*x_`n'
	gen w_co_1_`n'=co_1_`n'*wa
	gen w_co_2_`n'=co_2_`n'*wa
	egen cuco_1_`n'=sum(w_co_1_`n')
	egen cuco_2_`n'=sum(w_co_2_`n')
	drop x_`n' co_*_`n' w_co_*_`n'
}

duplicates drop cuco_1_1, force
keep cuco_*
gen id=1

reshape long cuco_1_ cuco_2_, i(id)
drop id
rename _j proportion
rename cuco_1_ cost_b1_1
rename cuco_2_ cost_b1_2
keep cost_b1_1 cost_b1_2 proportion
gen atc=1
preserve
keep cost_b1_1 proportion atc
rename cost_b1_1 cost_b1
gen curve=0
save `temp'/temp.dta, replace
restore
keep cost_b1_2 proportion atc
rename cost_b1_2 cost_b1
gen curve=1
append using `temp'/temp.dta
append using `data'/costs.dta
gen costscatter=1 if designed!=.
keep proportion cost_b1 atc costscatter designed curve `var'

twoway ///
	(scatter cost_b1 proportion if costscatter==1 & `var'==0, msize(medsmall) mlcolor(gs4) mfcolor(gs7) msymbol(O)) ///
	(scatter cost_b1 proportion if costscatter==1 & `var'==1, msize(medsmall) mlcolor(gs8) mfcolor(gs11) msymbol(O)) ///
	(line cost_b1 proportion if atc==1 & curve==0, sort lcolor(gs7) lwidth(thick) lpattern(solid)) ///
	(line cost_b1 proportion if atc==1 & curve==1, sort lcolor(gs11) lwidth(thick) lpattern(solid)), ///
	ytitle("ATC per connection (USD)", size(medium) margin(medsmall)) ///
	xtitle("Community coverage (%)", size(medium) margin(medsmall)) ///
	ylabel(0(1000)6000, labsize(medsmall)) ///
	xlabel(0(10)100, labsize(medsmall)) ///
	xsize(9) ysize(11) ///
	legend(size(medsmall) label(1 "Siaya county") label(2 "Busia county") ///
	label(3 "Siaya county ATC curve") label(4 "Busia county ATC curve") ///
	order(3 4 1 2) symx(9) region(lstyle(none)) position(1) ring(0) rows(4)) ///
	scheme(s1mono)
	graph export `appDir'/figb10d.png, replace

***************
*** Panel E ***
***************
	
local var "market"

use `data'/costs.dta, clear
gen costscatter=1
replace connections=connect_d if connections==.

nl (cost_b1={b0}/connections+{b1}+{b2}*connections) if `var'==0
mat b=e(b)
matlist b
local b0=b[1,1]
local b1=b[1,2]
local b2=b[1,3]

nl (cost_b1={b0}/connections+{b1}+{b2}*connections) if `var'==1
mat b=e(b)
matlist b
local c0=b[1,1]
local c1=b[1,2]
local c2=b[1,3]

egen totc=sum(compounds)
gen wa=compounds/totc
drop totc

forval n = 1/100 {
	gen x_`n'=`n'/100*ucompounds
	gen co_1_`n'=`b0'/x_`n'+`b1'+`b2'*x_`n'
	gen co_2_`n'=`c0'/x_`n'+`c1'+`c2'*x_`n'
	gen w_co_1_`n'=co_1_`n'*wa
	gen w_co_2_`n'=co_2_`n'*wa
	egen cuco_1_`n'=sum(w_co_1_`n')
	egen cuco_2_`n'=sum(w_co_2_`n')
	drop x_`n' co_*_`n' w_co_*_`n'
}

duplicates drop cuco_1_1, force
keep cuco_*
gen id=1

reshape long cuco_1_ cuco_2_, i(id)
drop id
rename _j proportion
rename cuco_1_ cost_b1_1
rename cuco_2_ cost_b1_2
keep cost_b1_1 cost_b1_2 proportion
gen atc=1
preserve
keep cost_b1_1 proportion atc
rename cost_b1_1 cost_b1
gen curve=0
save `temp'/temp.dta, replace
restore
keep cost_b1_2 proportion atc
rename cost_b1_2 cost_b1
gen curve=1
append using `temp'/temp.dta
append using `data'/costs.dta
gen costscatter=1 if designed!=.
keep proportion cost_b1 atc costscatter designed curve `var'

twoway ///
	(scatter cost_b1 proportion if costscatter==1 & `var'==0, msize(medsmall) mlcolor(gs4) mfcolor(gs7) msymbol(O)) ///
	(scatter cost_b1 proportion if costscatter==1 & `var'==1, msize(medsmall) mlcolor(gs8) mfcolor(gs11) msymbol(O)) ///
	(line cost_b1 proportion if atc==1 & curve==0, sort lcolor(gs7) lwidth(thick) lpattern(solid)) ///
	(line cost_b1 proportion if atc==1 & curve==1, sort lcolor(gs11) lwidth(thick) lpattern(solid)), ///
	ytitle("ATC per connection (USD)", size(medium) margin(medsmall)) ///
	xtitle("Community coverage (%)", size(medium) margin(medsmall)) ///
	ylabel(0(1000)6000, labsize(medsmall)) ///
	xlabel(0(10)100, labsize(medsmall)) ///
	xsize(9) ysize(11) ///
	legend(size(medsmall) label(1 "Non-market community") label(2 "Market community") ///
	label(3 "Non-market community ATC curve") label(4 "Market community ATC curve") ///
	order(3 4 1 2) symx(9) region(lstyle(none)) position(1) ring(0) rows(4)) ///
	scheme(s1mono)
	graph export `appDir'/figb10e.png, replace

***************
*** Panel F ***
***************
	
local var "funded"

use `data'/costs.dta, clear
gen costscatter=1
replace connections=connect_d if connections==.

nl (cost_b1={b0}/connections+{b1}+{b2}*connections) if `var'==1
mat b=e(b)
matlist b
local b0=b[1,1]
local b1=b[1,2]
local b2=b[1,3]

nl (cost_b1={b0}/connections+{b1}+{b2}*connections) if `var'==0
mat b=e(b)
matlist b
local c0=b[1,1]
local c1=b[1,2]
local c2=b[1,3]

egen totc=sum(compounds)
gen wa=compounds/totc
drop totc

forval n = 1/100 {
	gen x_`n'=`n'/100*ucompounds
	gen co_1_`n'=`b0'/x_`n'+`b1'+`b2'*x_`n'
	gen co_2_`n'=`c0'/x_`n'+`c1'+`c2'*x_`n'
	gen w_co_1_`n'=co_1_`n'*wa
	gen w_co_2_`n'=co_2_`n'*wa
	egen cuco_1_`n'=sum(w_co_1_`n')
	egen cuco_2_`n'=sum(w_co_2_`n')
	drop x_`n' co_*_`n' w_co_*_`n'
}

duplicates drop cuco_1_1, force
keep cuco_*
gen id=1

reshape long cuco_1_ cuco_2_, i(id)
drop id
rename _j proportion
rename cuco_1_ cost_b1_1
rename cuco_2_ cost_b1_2
keep cost_b1_1 cost_b1_2 proportion
gen atc=1
preserve
keep cost_b1_1 proportion atc
rename cost_b1_1 cost_b1
gen curve=0
save `temp'/temp.dta, replace
restore
keep cost_b1_2 proportion atc
rename cost_b1_2 cost_b1
gen curve=1
append using `temp'/temp.dta
append using `data'/costs.dta
gen costscatter=1 if designed!=.
keep proportion cost_b1 atc costscatter designed curve `var'

twoway ///
	(scatter cost_b1 proportion if costscatter==1 & `var'==1, msize(medsmall) mlcolor(gs4) mfcolor(gs7) msymbol(O)) ///
	(scatter cost_b1 proportion if costscatter==1 & `var'==0, msize(medsmall) mlcolor(gs8) mfcolor(gs11) msymbol(O)) ///
	(line cost_b1 proportion if atc==1 & curve==0, sort lcolor(gs7) lwidth(thick) lpattern(solid)) ///
	(line cost_b1 proportion if atc==1 & curve==1, sort lcolor(gs11) lwidth(thick) lpattern(solid)), ///
	ytitle("ATC per connection (USD)", size(medium) margin(medsmall)) ///
	xtitle("Community coverage (%)", size(medium) margin(medsmall)) ///
	ylabel(0(1000)6000, labsize(medsmall)) ///
	xlabel(0(10)100, labsize(medsmall)) ///
	xsize(9) ysize(11) ///
	legend(size(medsmall) label(1 "Community funded early on") label(2 "Community funded later") ///
	label(3 "Community funded early on ATC curve") label(4 "Community funded later ATC curve") ///
	order(3 4 1 2) symx(9) region(lstyle(none)) position(1) ring(0) rows(4)) ///
	scheme(s1mono)
	graph export `appDir'/figb10f.png, replace

rm `temp'/temp.dta
	
********************************************************************************
* Figure B11: Social surplus implications under various demand curve assumptions
********************************************************************************

***************
*** Panel A ***
***************

use `data'/costs.dta, clear
gen costscatter=1
replace connections=connect_d if connections==.

nl (cost_b1={b0}/connections+{b1}+{b2}*connections)
mat b=e(b)
matlist b
local b0=b[1,1]
local b1=b[1,2]
local b2=b[1,3]

egen totc=sum(compounds)
gen wa=compounds/totc
drop totc

forval n = 1/100 {
	gen x_`n'=`n'/100*ucompounds
	gen co_`n'=`b0'/x_`n'+`b1'+`b2'*x_`n'
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
gen atc=1
save `temp'/atc.dta, replace
append using `data'/costs.dta
gen costscatter=1 if designed!=.
keep proportion cost_b1 atc costscatter designed
save `temp'/cost_data.dta, replace

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

sum cost_b1 if proportion==100 & atc==1, detail
* ATC at 100% is $ 739.29
local line=`r(mean)'
gen masscost=`line'
di `line'*84.7

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
	graph export `rootDir'/appendix/figb11a.png, replace

***************
*** Panel B ***
***************

use `temp'/cost_data.dta, clear
append using `temp'/demand_data.dta
replace price=cost_b1 if price==.

sum cost_b1 if proportion==100 & atc==1, detail
* ATC at 100% is $ 739.2913
local line=`r(mean)'
gen masscost=`line'
display `line'*84.7
* [0, 100] area under the lowest cost is $62,618 at 84.7 average community compound density.

* Extend demand curve through [0, 1.3] with intercept at $3000

display (95.52631-23.7467)*170.5708*0.5 + ///
	(23.7467-7.105263)*(170.5708) + ///
	(23.7467-7.105263)*(284.2847-170.5708)*0.5 + ///
	(7.105263-1.304348)*(397.9986-284.2847)*0.5 + ///
	(7.105263-1.304348)*(284.2847) + ///
	(1.304348)*(3000-397.9986)*0.5 + ///
	(1.304348)*(397.9986)
		
* Measured [1.3, 100] area under the demand curve is $11885.411
* Implied [0, 100] area is $14101.497.
* for now extend to zero

* What does the intercept need to be to overturn the results?

display ((95.52631-23.7467)*170.5708*0.5 + ///
	(23.7467-7.105263)*(170.5708) + ///
	(23.7467-7.105263)*(284.2847-170.5708)*0.5 + ///
	(7.105263-1.304348)*(397.9986-284.2847)*0.5 + ///
	(7.105263-1.304348)*(284.2847) + ///
	(1.304348)*(3000-397.9986)*0.5 - 62617.973)/-(1.304348)

* In order to overturn the results (e.g. costs > demand surplus), the y-intercept needs to be $37,594.

gen shadedcost=demandcurve
set obs 334
replace proportion=100 in 331
replace masscost=`line' in 331
replace shadedcost=1 in 331
replace proportion=0 in 332
replace masscost=`line' in 332
replace shadedcost=1 in 332

gen shadeddemand=demandcurve
replace proportion=0 in 334
replace price=3000 in 334
replace shadeddemand=1 in 334

sort proportion shadedcost shadeddemand

twoway ///
	(area masscost proportion if shadedcost==1, lpattern(solid) lstyle(unextended) lwidth(vthin) lcolor(gs13) fcolor(gs13)) ///
	(area price proportion if shadeddemand==1, lpattern(solid)  lstyle(unextended) lwidth(vthin) lcolor(gs11) fcolor(gs11)) ///
	(line price proportion if atc==1 & price<=3000  & proportion <=100, lcolor(gs8) lwidth(thick) lpattern(solid)) ///
	(line price proportion if demandcurve==1 & proportion <=100, lcolor(gs5) lwidth(thick) lpattern(solid)), ///
	xtitle("Take-up, community coverage (%)", size(medium) margin(medsmall)) ///
	ytitle("Connection price, ATC per connection (USD)", size(medium) margin(medsmall)) ///
	ylabel(0(500)3000, labsize(medsmall)) ///
	xlabel(0(10)100, labsize(medsmall)) ///
	xsize(9) ysize(11) ///
	plotr(margin(0.3 0.3 0.4 3)) ///
	text(650 70 "TC = $62,618", just(left) place(e) size(medsmall)) ///
	text(80 2 "CS = $14,101", just(left) place(e) size(medsmall)) ///
	legend(size(medsmall) label(1 "Area 1") label(2 "Area 2") label(3 "ATC curve") ///
	label(4 "Demand curve") order(4 3) symx(9) region(lstyle(none)) position(1) ring(0) rows(5)) ///
	scheme(s1mono)
	graph export `appDir'/figb11b.png, replace

display (62617.973-14101.497)
* $48,516 in total
display (62617.973-14101.497)/84.7
* $573 per household

***************
*** Panel C ***
***************

* Extend demand curve through [0, 1.3] with intercept at $3000; Step-shaped demand curve

display (95.52631-23.7467)*170.5708 + ///
	(23.7467-7.105263)*284.2847 + ///
	(7.105263-1.304348)*397.9986 + ///
	(1.304348)*3000

* Implied [0, 100] area is $23196.211

drop if demandcurve==1

set obs 440

replace proportion=95.52631 in 432
replace demandcurve=1 in 432
replace price=0 in 432

replace proportion=95.52631 in 433
replace demandcurve=1 in 433
replace price=170.5708 in 433

replace proportion=23.7467 in 434
replace demandcurve=1 in 434
replace price=170.5708 in 434

replace proportion=23.7467 in 435
replace demandcurve=1 in 435
replace price=284.2847 in 435

replace proportion=7.105263 in 436
replace demandcurve=1 in 436
replace price=284.2847 in 436

replace proportion=7.105263 in 437
replace demandcurve=1 in 437
replace price=397.9986 in 437

replace proportion=1.304348 in 438
replace demandcurve=1 in 438
replace price=397.9986 in 438

replace proportion=1.304348 in 439
replace demandcurve=1 in 439
replace price=3000 in 439

replace proportion=0 in 440
replace demandcurve=1 in 440
replace price=3000 in 440

drop shadeddemand
gen shadeddemand=demandcurve

twoway ///
	(area masscost proportion if shadedcost==1, lpattern(solid) lstyle(unextended) lwidth(vthin) lcolor(gs13) fcolor(gs13)) ///
	(area price proportion if shadeddemand==1, lpattern(solid)  lstyle(unextended) lwidth(vthin) lcolor(gs11) fcolor(gs11)) ///
	(line price proportion if atc==1 & price<=3000 & proportion <=100, lcolor(gs8) lwidth(thick) lpattern(solid)) ///
	(line price proportion if demandcurve==1 & proportion <=100, lcolor(gs5) lwidth(thick) lpattern(solid)), ///
	xtitle("Take-up, community coverage (%)", size(medium) margin(medsmall)) ///
	ytitle("Connection price, ATC per connection (USD)", size(medium) margin(medsmall)) ///
	ylabel(0(500)3000, labsize(medsmall)) ///
	xlabel(0(10)100, labsize(medsmall)) ///
	xsize(9) ysize(11) ///
	plotr(margin(0.3 0.3 0.4 3)) ///
	text(650 70 "TC = $62,618", just(left) place(e) size(medsmall)) ///
	text(80 2 "CS = $23,196", just(left) place(e) size(medsmall)) ///
	legend(size(medsmall) label(1 "Area 1") label(2 "Area 2") label(3 "ATC curve") ///
	label(4 "Demand curve") order(4 3) symx(9) region(lstyle(none)) position(1) ring(0) rows(5)) ///
	scheme(s1mono)
	graph export `appDir'/figb11c.png, replace

display (62617.973-23196.211)
* $39,422 in total
display (62617.973-23196.211)/84.7
* $465 per household

rm `temp'/atc.dta
rm `temp'/demand_data.dta
rm `temp'/cost_data.dta

********************************************************************************
* Figure B12: Social surplus implications of a government program
********************************************************************************

use `data'/costs.dta, clear
gen costscatter=1
replace connections=connect_d if connections==.

nl (cost_b1={b0}/connections+{b1}+{b2}*connections)
mat b=e(b)
matlist b
local b0=b[1,1]
local b1=b[1,2]
local b2=b[1,3]

egen totc=sum(compounds)
gen wa=compounds/totc
drop totc

forval n = 1/100 {
	gen x_`n'=`n'/100*ucompounds
	gen co_`n'=`b0'/x_`n'+`b1'+`b2'*x_`n'
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
gen atc=1
save `temp'/atc.dta, replace
append using `data'/costs.dta
gen costscatter=1 if designed!=.
keep proportion cost_b1 atc costscatter designed
save `temp'/cost_data.dta, replace

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
replace price=cost_b1 if price==.

sum cost_b1 if proportion==100 & atc==1, detail
local line=`r(mean)'
gen masscost=`line'
display `line'*84.7

gen shadedcost=demandcurve
set obs 334
replace proportion=100 in 331
replace masscost=`line' in 331
replace shadedcost=1 in 331
replace proportion=0 in 332
replace masscost=`line' in 332
replace shadedcost=1 in 332

gen shadeddemand=demandcurve
replace proportion=0 in 334
replace price=423.56741 in 334
replace shadeddemand=1 in 334

sort proportion shadedcost shadeddemand
* Demand curve at $170.5708 = 23.7467%
* ATC at 23.7467% is approximately $1057.3588
di 1063.599 - (1055.242-1063.599)/(24-23)*23
* 1255.81
di (1055.242-1063.599)/(24-23)*23.7467+1255.81

gen lmcpcost=1057.3588
set obs 436
replace proportion=23.7467 in 436
replace shadedcost=1 in 436

display (23.7467-7.105263)*(170.5708) + ///
	(23.7467-7.105263)*(284.2847-170.5708)*0.5 + ///
	(7.105263-1.304348)*(397.9986-284.2847)*0.5 + ///
	(7.105263-1.304348)*(284.2847) + ///
	(1.304348)*(423.56741-397.9986)*0.5 + ///
	(1.304348)*(397.9986)
* Implied [0, 23.7] area is 6299.4623

display 23.7467*1057.3588
* Implied total cost is $25,109

twoway ///
	(area lmcpcost proportion if shadedcost==1 & proportion<=24, lpattern(solid) lstyle(unextended) lwidth(vthin) lcolor(gs13) fcolor(gs13)) ///
	(area price proportion if shadeddemand==1 & proportion<=24, lpattern(solid)  lstyle(unextended) lwidth(vthin) lcolor(gs11) fcolor(gs11)) ///
	(line price proportion if atc==1 & price<=3000 & proportion <=100, lcolor(gs8) lwidth(thick) lpattern(solid)) ///
	(line price proportion if demandcurve==1 & proportion <=100, lcolor(gs5) lwidth(thick) lpattern(solid)), ///
	xtitle("Take-up, community coverage (%)", size(medium) margin(medsmall)) ///
	ytitle("Connection price, ATC per connection (USD)", size(medium) margin(medsmall)) ///
	ylabel(0(500)3000, labsize(medsmall)) ///
	xlabel(0(20)100, labsize(medsmall)) ///
	xsize(9) ysize(11) ///
	plotr(margin(0.3 0.3 0.4 3)) ///
	text(980 1 "TC = $25,109", just(left) place(e) size(small)) ///
	text(70 1 "CS = $6,299", just(left) place(e) size(small)) ///
	legend(size(medsmall) label(1 "Area 1") label(2 "Area 2") label(3 "ATC curve") ///
	label(4 "Demand curve") order(4 3) symx(9) region(lstyle(none)) position(1) ring(0) rows(5)) ///
	scheme(s1mono)
	graph export `appDir'/figb12.png, replace

display (25108.782-6299.4623)
* $18,809 gap between cost and demand
display (25108.782-6299.4623)/(84.7*0.237467)
* $935.16 per household

rm `temp'/atc.dta
rm `temp'/demand_data.dta
rm `temp'/cost_data.dta

********************************************************************************
* Figure B13: Demand, comparing households with bank accounts and high quality walls
********************************************************************************

***************
*** Panel A ***
***************

* CV question 1 (16b): "Would you be willing to pay XX Ksh for an electricity connection?"
use `data'/demand.dta, clear
keep if wall==0 & bank==0
bysort WTP_amt: egen takeup1=mean(WTP_r1)
drop if WTP_amt==94 | WTP_amt==95 | WTP_amt==5000 | WTP_amt==30000
* not enough observations (we asked about a couple different data points in the first week of surveying (randomly) before settling on set of price points)
* 119 observations in this group
duplicates drop WTP_amt, force
keep WTP_amt takeup1
gen cv1=takeup1*100
drop takeup1
rename WTP_amt kprice
sort kprice
save `temp'/cv1.dta, replace

* CV question 2 (16c): "16c. Imagine that you were offered an electricity connection at this price today, and you were given 6 weeks to complete the payment. Would you accept the offer?"
use `data'/demand.dta, clear
keep if wall==0 & bank==0
bysort WTP_amt: egen takeup2=mean(WTP_r2)
drop if WTP_amt==94 | WTP_amt==95 | WTP_amt==5000 | WTP_amt==30000
* not enough observations (we asked about a couple different data points in the first week of surveying (randomly) before settling on set of price points)
* 119 observations in this group
duplicates drop WTP_amt, force
keep WTP_amt takeup2
gen cv2=takeup2*100
drop takeup2
rename WTP_amt kprice
sort kprice
save `temp'/cv2.dta, replace

* Experimental results
use `data'/demand.dta, clear
keep if wall==0 & bank==0
bysort price: egen takeup3=mean(takeup)
* 236 observations in this group
duplicates drop price, force
drop if price==.
keep price takeup3
gen rp=takeup3*100
drop takeup3
rename price kprice
sort kprice
save `temp'/rp.dta, replace

* Merge files
use `temp'/cv1.dta, clear
merge kprice using `temp'/cv2.dta
drop _merge 
sort kprice
merge kprice using `temp'/rp.dta
gen price=kprice/`rate'
gen experimentpoint=1
drop _merge
sort kprice
save `temp'/experimental_takeup.dta, replace
rm `temp'/cv1.dta
rm `temp'/cv2.dta
rm `temp'/rp.dta

twoway ///
		(connected price cv1, lcolor(gs2) lwidth(medthin) lpattern(longdash) msize(small) mcolor(gs1) msymbol(square)) ///
		(connected price cv2, lcolor(gs10) lwidth(medthin) lpattern(longdash) msize(small) mcolor(gs8) msymbol(square)) ///
		(connected price rp if experimentpoint==1, lcolor(gs5) lwidth(thick) lpattern(solid) msize(small) mcolor(gs2) msymbol(square)), ///
		xtitle("Take-up (%)", size(medium) margin(medsmall)) ///
		ytitle("Connection price (USD)", size(medium) margin(medsmall)) ///
		xlabel(0(20)100, labsize(medsmall)) ///
		ylabel(0(100)900, labsize(medsmall)) ///
		xsize(9) ysize(11) ///
		scheme(s1mono) ///
		legend(size(medsmall) label( 1 "WTP") ///
		label(2 "WTP (time-limited offer)") ///
		label( 3 "Demand curve") ///
		symx(9) order(3 1 2) region(lstyle(none)) margin(zero) position(1) ring(0) rows(3))
		graph export `appDir'/figb13a.png, replace

***************
*** Panel B ***
***************

* CV question 1 (16b): "Would you be willing to pay XX Ksh for an electricity connection?"
use `data'/demand.dta, clear
keep if wall==1 & bank==1
bysort WTP_amt: egen takeup1=mean(WTP_r1)
drop if WTP_amt==94 | WTP_amt==95 | WTP_amt==5000 | WTP_amt==30000
* not enough observations (we asked about a couple different data points in the first week of surveying (randomly) before settling on set of price points)
* 119 observations in this group
duplicates drop WTP_amt, force
keep WTP_amt takeup1
gen cv1=takeup1*100
drop takeup1
rename WTP_amt kprice
sort kprice
save `temp'/cv1.dta, replace

* CV question 2 (16c): "16c. Imagine that you were offered an electricity connection at this price today, and you were given 6 weeks to complete the payment. Would you accept the offer?"
use `data'/demand.dta, clear
keep if wall==1 & bank==1
bysort WTP_amt: egen takeup2=mean(WTP_r2)
drop if WTP_amt==94 | WTP_amt==95 | WTP_amt==5000 | WTP_amt==30000
* not enough observations (we asked about a couple different data points in the first week of surveying (randomly) before settling on set of price points)
* 119 observations in this group
duplicates drop WTP_amt, force
keep WTP_amt takeup2
gen cv2=takeup2*100
drop takeup2
rename WTP_amt kprice
sort kprice
save `temp'/cv2.dta, replace

* Experimental results
use `data'/demand.dta, clear
keep if wall==1 & bank==1
bysort price: egen takeup3=mean(takeup)
* 236 observations in this group
duplicates drop price, force
drop if price==.
keep price takeup3
gen rp=takeup3*100
drop takeup3
rename price kprice
sort kprice
save `temp'/rp.dta, replace

* Merge files
use `temp'/cv1.dta, clear
merge kprice using `temp'/cv2.dta
drop _merge 
sort kprice
merge kprice using `temp'/rp.dta
gen price=kprice/`rate'
gen experimentpoint=1
drop _merge
sort kprice
save `temp'/experimental_takeup.dta, replace
rm `temp'/cv1.dta
rm `temp'/cv2.dta
rm `temp'/rp.dta

twoway ///
		(connected price cv1, lcolor(gs2) lwidth(medthin) lpattern(longdash) msize(small) mcolor(gs1) msymbol(square)) ///
		(connected price cv2, lcolor(gs10) lwidth(medthin) lpattern(longdash) msize(small) mcolor(gs8) msymbol(square)) ///
		(connected price rp if experimentpoint==1, lcolor(gs5) lwidth(thick) lpattern(solid) msize(small) mcolor(gs2) msymbol(square)), ///
		xtitle("Take-up (%)", size(medium) margin(medsmall)) ///
		ytitle("Connection price (USD)", size(medium) margin(medsmall)) ///
		xlabel(0(20)100, labsize(medsmall)) ///
		ylabel(0(100)900, labsize(medsmall)) ///
		xsize(9) ysize(11) ///
		scheme(s1mono) ///
		legend(size(medsmall) label( 1 "WTP") ///
		label(2 "WTP (time-limited offer)") ///
		label( 3 "Demand curve") ///
		symx(9) order(3 1 2) region(lstyle(none)) margin(zero) position(1) ring(0) rows(3))
		graph export `appDir'/figb13b.png, replace

rm `temp'/experimental_takeup.dta

********************************************************************************
* Figure B15: Discrepancies in costs and poles, by contractor
********************************************************************************

use `data'/costs.dta, clear
gen costscatter=1
drop if contractor==""
bysort contractor: egen design_poles=sum(poles_d)
bysort contractor: egen actual_poles=sum(poles_a)
bysort contractor: egen budget_cost=sum(sitecost_b)
bysort contractor: egen invoice_cost=sum(sitecost_a)
bysort contractor: egen connects=sum(connect_a)
gen pole_discrepancy=(actual_poles-design_poles)/design_poles*100
gen cost_discrepancy=(invoice_cost-budget_cost)/budget_cost*100

egen tot_design_poles=sum(poles_d)
egen tot_actual_poles=sum(poles_a)
egen tot_budget_cost=sum(sitecost_b)
egen tot_invoice_cost=sum(sitecost_a)
gen tot_pole_discrepancy=(tot_actual_poles-tot_design_poles)/tot_design_poles*100
gen tot_cost_discrepancy=(tot_invoice_cost-tot_budget_cost)/tot_budget_cost*100

duplicates drop contractor, force
keep contractor design_poles actual_poles invoice_cost budget_cost connects pole_discrepancy cost_discrepancy tot_*

gen rea_owned=0
replace rea_owned=1 if	contractor=="Contractor 4" | contractor=="Contractor 9" | contractor=="Contractor 11" |  ///
						contractor=="Contractor 2" | contractor=="Contractor 10" | contractor=="Contractor 14" | ///
						contractor=="Contractor 12"

twoway ///
	(scatter pole_discrepancy cost_discrepancy [weight=connects] if rea_owned==0, msymbol(O) msize(medlarge) mlcolor(gs2) mfcolor(gs12)) ///
	(scatter pole_discrepancy cost_discrepancy [weight=connects] if rea_owned==1, msymbol(O) msize(medlarge) mlcolor(gs2) mfcolor(gs12)), ///
	ytitle("Difference between actual and budgeted poles (%)", size(medsmall) margin(medsmall)) ///
	xtitle("Difference between invoiced and budgeted costs (%)", size(medsmall) margin(medsmall)) ///
	yline(0, lpattern(solid) lwidth(thin) lcolor(gs1)) xline(0, lpattern(solid) lwidth(thin) lcolor(gs1)) ///
	yline(-21.3, lpattern(dash) lwidth(thin) lcolor(gs6)) xline(1.7, lpattern(dash) lwidth(thin) lcolor(gs6)) ///
	xlabel(-50(10)20, labsize(small)) ///
	ylabel(-50(10)20, labsize(small)) ///
	xsize(18) ysize(15) ///
	text(-15 -49.6 "Average" "discrepancy" "in poles: -21.3%", color(gs6) just(left) place(e) size(small)) ///
	text(6 3.5 "Average" "discrepancy" "in costs: +1.7%", color(gs6) just(left) place(e) size(small)) ///
	legend(off) ///
	scheme(s1mono)
	graph export `appDir'/figb15.png, replace

		
