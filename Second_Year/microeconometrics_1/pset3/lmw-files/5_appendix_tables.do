clear all
set more off
local rootDir = "/Users/Research/LMW2019"

local data "`rootDir'/data/stata"
local temp "`rootDir'/data/temp"
local appDir "`rootDir'/appendix"

local rate=87.94

********************************************************************************
* Appendix Tables
* Last updated: 23 Apr 2019
********************************************************************************

********************************************************************************
* Table A1A: Impact of scale on ATC per connection, sample
********************************************************************************

use `data'/costs.dta, clear
drop if connect_a==.

gen atc_usd=cost_b1
drop if atc_usd==.

gen m1=connect_d
gen m2=m1^2
gen m3=m1^3
gen q1=proportion
gen q2=q1^2
gen q3=q1^3

egen mean_terrain_2=mean(terrain_2)
gen terrain=terrain_2-mean_terrain_2

egen mean_population=mean(population)
gen population2=(population-mean_population)
drop population
rename population2 population

egen mean_compounds=mean(compounds)
gen compounds1=(compounds-mean_compounds)
drop compounds
rename compounds1 compounds

gen interact_t1=terrain*m1
gen interact_t2=terrain*m2
gen variable_t=terrain

gen interact_p1=population*m1
gen interact_p2=population*m2/100
gen variable_p=population

gen interact_c1=compounds*m1
gen interact_c2=compounds*m2/100
gen variable_c=compounds

local controlsList ///
	busia market funded electrification population round_trip terrain

* model 1
reg atc_usd m1 m2, vce(robust)
	qui estadd local controls "No", replace
	qui estadd ysumm
	qui estimates store m1

* model 2
reg atc_usd m1 m2 ///
	busia market funded electrification population round_trip terrain, vce(robust)
	qui estadd local controls "Yes", replace
	qui estadd ysumm
	qui estimates store m2

* model 3
reg atc_usd q1 q2, vce(robust)
	qui estadd local controls "No", replace
	qui estadd ysumm
	qui estimates store m3

* model 4
reg atc_usd q1 q2 ///
	busia market funded electrification population round_trip terrain, vce(robust)
	qui estadd local controls "Yes", replace
	qui estadd ysumm
	qui estimates store m4

* Model 5
ivregress 2sls atc_usd (m1 m2 = tx2 tx3), vce(robust)
	qui estadd local controls "No", replace
	qui estadd ysumm
	qui estimates store m5

* Model 6
ivregress 2sls atc_usd (m1 m2 = tx2 tx3) ///
	busia market funded electrification population round_trip terrain, vce(robust)
	qui estadd local controls "Yes", replace
	qui estadd ysumm
	qui estimates store m6

* Model 7
ivregress 2sls atc_usd (q1 q2 = tx2 tx3), vce(robust)
	qui estadd local controls "No", replace
	qui estadd ysumm
	qui estimates store m7

* Model 8
ivregress 2sls atc_usd (q1 q2 = tx2 tx3) ///
	busia market funded electrification population round_trip terrain, vce(robust)
	qui estadd local controls "Yes", replace
	qui estadd ysumm
	qui estimates store m8
/*
esttab m1 m2 m3 m4 m5 m6 m7 m8 using `appDir'/taba1a.tex, ///
	replace b(%10.1f) se scalars( ///
	"ymean Mean of dep. variable (USD)" ///
	"N Observations" "r2 R$^2$") ///
	sfmt(%10.0f %10.0fc %10.2f) ///
	star(* 0.10 ** 0.05 *** 0.01) noobs label nomtitles ///
	align(cccccc) width(\hsize) nonotes ///
	varlabels(m1 "Number of connections (M)" m2 "M$^2$" ///
	q1 "Community coverage (Q)" q2 "Q$^2$" ///
	terrain "Land gradient" interact_t1 "Land gradient $\times$ M" interact_t2 "Land gradient $\times$ M$^2$" ///
	compounds "Households" interact_c1 "Households $\times$ M" interact_c2 "Households $\times$ M$^2$ / 100" ///
	population "Community population" interact_p1 "Community population $\times$ M" interact_p2 "Community population $\times$ M$^2$ / 100" ///
	busia "Busia=1" market "Market transformer=1" funded "Transformer funded early on=1" ///
	electrification "Community electrification rate" ///
	round_trip "Round-trip distance to REA (km)") ///	
	order(m1 m2 q1 q2 busia market funded electrification population round_trip terrain) ///
	keep(m1 m2 q1 q2 busia market funded electrification population round_trip terrain) ///
	title(\textemdash Impact of scale on average total cost (ATC) per connection, sample communities}\vspace{-0.5em)
* Note: Latex file is manually edited for formatting purposes
*/

********************************************************************************
* Table A1B: Impact of scale on ATC per connection, sample & designed
********************************************************************************

use `data'/costs.dta, clear

gen atc_usd=cost_b1
drop if atc_usd==.

gen m1=connect_d
gen m2=m1^2
gen m3=m1^3
gen q1=proportion
gen q2=q1^2
gen q3=q1^3

egen mean_terrain_2=mean(terrain_2)
gen terrain=terrain_2-mean_terrain_2

egen mean_population=mean(population)
gen population2=(population-mean_population)
drop population
rename population2 population

egen mean_compounds=mean(compounds)
gen compounds1=(compounds-mean_compounds)
drop compounds
rename compounds1 compounds

gen interact_t1=terrain*m1
gen interact_t2=terrain*m2
gen variable_t=terrain

gen interact_p1=population*m1
gen interact_p2=population*m2/100
gen variable_p=population

gen interact_c1=compounds*m1
gen interact_c2=compounds*m2/100
gen variable_c=compounds

local controlsList ///
	busia market funded electrification population round_trip terrain

* model 1
reg atc_usd m1 m2, vce(robust)
	qui estadd local controls "No", replace
	qui estadd ysumm
	qui estimates store m1

* model 2
reg atc_usd m1 m2 ///
	busia market funded electrification population round_trip terrain, vce(robust)
	qui estadd local controls "Yes", replace
	qui estadd ysumm
	qui estimates store m2

* model 3
reg atc_usd q1 q2, vce(robust)
	qui estadd local controls "No", replace
	qui estadd ysumm
	qui estimates store m3

* model 4
reg atc_usd q1 q2 ///
	busia market funded electrification population round_trip terrain, vce(robust)
	qui estadd local controls "Yes", replace
	qui estadd ysumm
	qui estimates store m4
/*
esttab m1 m2 m3 m4 using `appDir'/taba1b.tex, ///
	replace b(%10.1f) se scalars( ///
	"ymean Mean of dep. variable (USD)" ///
	"N Observations" "r2 R$^2$") ///
	sfmt(%10.0f %10.0fc %10.2f) ///
	star(* 0.10 ** 0.05 *** 0.01) noobs label nomtitles ///
	align(cccccc) width(\hsize) nonotes ///
	varlabels(m1 "Number of connections (M)" m2 "M$^2$" ///
	q1 "Community coverage (Q)" q2 "Q$^2$" ///
	terrain "Land gradient" interact_t1 "Land gradient $\times$ M" interact_t2 "Land gradient $\times$ M$^2$" ///
	compounds "Households" interact_c1 "Households $\times$ M" interact_c2 "Households $\times$ M$^2$ / 100" ///
	population "Community population" interact_p1 "Community population $\times$ M" interact_p2 "Community population $\times$ M$^2$ / 100" ///
	busia "Busia=1" market "Market transformer=1" funded "Transformer funded early on=1" ///
	electrification "Community electrification rate" ///
	round_trip "Round-trip distance to REA (km)") ///	
	order(m1 m2 q1 q2 busia market funded electrification population round_trip terrain) ///
	keep(m1 m2 q1 q2 busia market funded electrification population round_trip terrain) ///
	title(\textemdash Impact of scale on average total cost (ATC) per connection, sample and designed communities}\vspace{-0.5em)
* Note: Latex file is manually edited for formatting purposes
*/
********************************************************************************
* Table A1C: Impact of scale on ATC per connection
********************************************************************************

use `data'/costs.dta, clear
drop if connect_d==.
gen sample1=0
replace sample1=1 if siteno<=2520

gen atc_usd=cost_b1
drop if atc_usd==.

gen m1=connect_d
gen m2=m1^2
gen m3=m1^3
gen q1=proportion
gen q2=q1^2
gen q3=q1^3

egen mean_terrain_2=mean(terrain_2)
gen terrain=terrain_2-mean_terrain_2

egen mean_population=mean(population)
gen population2=(population-mean_population)
drop population
rename population2 population

egen mean_compounds=mean(compounds)
gen compounds1=(compounds-mean_compounds)
drop compounds
rename compounds1 compounds

gen interact_t1=terrain*m1
gen interact_t2=terrain*m2
gen variable_t=terrain

gen interact_p1=population*m1
gen interact_p2=population*m2/100
gen variable_p=population

gen interact_c1=compounds*m1
gen interact_c2=compounds*m2/100
gen variable_c=compounds

local controlsList ///
	busia market funded electrification population round_trip terrain

* model 1
reg atc_usd m1 m2 ///
	busia market funded electrification population round_trip terrain, vce(robust)
	qui estadd local controls "Yes", replace
	qui estadd ysumm
	qui estimates store m1

* model 2 - population interaction
reg atc_usd m1 m2 ///
	variable_p interact_p1 interact_p2 ///
	busia market funded electrification round_trip terrain, vce(robust)
	qui estadd local controls "Yes", replace
	qui estadd ysumm
	qui estimates store m2

* model 3 - terrain interaction
reg atc_usd m1 m2 ///
	variable_t interact_t1 interact_t2 ///
	busia market funded electrification population round_trip, vce(robust)
	qui estadd local controls "Yes", replace
	qui estadd ysumm
	qui estimates store m3

* model 4 - compounds interaction
reg atc_usd m1 m2 ///
	variable_c interact_c1 interact_c2 ///
	busia market funded electrification population round_trip terrain, vce(robust)
	qui estadd local controls "Yes", replace
	qui estadd ysumm
	qui estimates store m4
/*
esttab m1 m2 m3 m4 using `appDir'/taba1c.tex, ///
	replace b(%10.1f) se scalars( ///
	"controls Community controls" ///
	"ymean Mean of dep. variable (USD)" ///
	"N Observations" "r2 R$^2$") ///
	sfmt(%6.0g %10.0f %10.0fc %10.2f) ///
	star(* 0.10 ** 0.05 *** 0.01) noobs label nomtitles ///
	align(cccccc) width(\hsize) nonotes ///
	varlabels(m1 "Number of connections (M)" m2 "M$^2$" ///
	variable_t "Land gradient" interact_t1 "Land gradient $\times$ M" interact_t2 "Land gradient $\times$ M$^2$" ///
	variable_c "Households" interact_c1 "Households $\times$ M" interact_c2 "Households $\times$ M$^2$ / 100" ///
	variable_p "Community population" interact_p1 "Community population $\times$ M" interact_p2 "Community population $\times$ M$^2$ / 100" ///
	) ///	
	order(m1 m2 variable_p interact_p1 interact_p2 variable_t interact_t1 interact_t2 variable_c interact_c1 interact_c2) ///
	keep(m1 m2 variable_p interact_p1 interact_p2 variable_t interact_t1 interact_t2 variable_c interact_c1 interact_c2) ///
	title(\textemdash Impact of scale on average total cost (ATC) per connection, sample and designed communities}\vspace{-0.5em)
* Note: Latex file is manually edited for formatting purposes
*/

********************************************************************************
* Table B2: Baseline summary statistics and randomization balance check
********************************************************************************

* Note: Results are produced and stored into a matrix.

use `data'/demand.dta, clear
drop if connected==1

local varList ///
	female age senior someschool spouse notfarmer employed awareness bank earn ///
	hhsize hhyouth wall ownland distance spend ///
	bednets bikes sofapieces chickens cattle radio television ///
	electrification population

* Column 1: Means for the control groups
preserve
keep if anytx==0
foreach var of local varList {
	qui tabstat `var', stat(mean sd) save
	qui matrix `var'=r(StatTotal)
}
matrix column1 = ///
	female \ age \ senior \ someschool \ spouse \ notfarmer \ employed \ awareness \ bank \ earn \ ///
	hhsize \ hhyouth \ wall \ ownland \ distance \ spend \ ///
	bednets \ bikes \ sofapieces \ chickens \ cattle \ radio \ television \ ///
	electrification \ population

*mat2txt, matrix(column1) saving(`appDir'/tabb2a.csv) replace
restore

mat list column1

* Columns 2-5: Regressions of collection of dependent variables on treatment indicators 
foreach var of local varList {
	qui reg `var' tx1 tx2 tx3 busia funded market, vce(cluster siteno)
	qui testnl _b[tx3] = _b[tx2] = _b[tx1] = 0	
	qui ereturn list
	qui mat co=e(b)
	qui mat va=e(V)
	qui matrix r1_`var' = co[1,1], co[1,2], co[1,3], r(F), e(N)
	qui matrix r2_`var' = va[1,1]^0.5, va[2,2]^0.5, va[3,3]^0.5, r(p), .	
	qui mat list r1_`var'
	qui mat list r2_`var'	
}

matrix regs1 = 	r1_female \ r2_female \ ///
				r1_age \ r2_age \ ///
				r1_senior \ r2_senior \ ///
				r1_someschool \ r2_someschool \ ///
				r1_spouse \ r2_spouse \ ///
				r1_notfarmer \ r2_notfarmer \ ///
				r1_employed \ r2_employed \ ///
				r1_awareness \ r2_awareness \ ///
				r1_bank \ r2_bank \ ///
				r1_earn \ r2_earn

matrix regs2 = 	r1_hhsize \ r2_hhsize \ ///
				r1_hhyouth \ r2_hhyouth \ ///
				r1_wall \ r2_wall \ ///
				r1_ownland \ r2_ownland \ ///
				r1_distance \ r2_distance \ ///
				r1_spend \ r2_spend

matrix regs3 = 	r1_bednets \ r2_bednets \ ///
				r1_bikes \ r2_bikes \ ///
				r1_sofapieces \ r2_sofapieces \ ///
				r1_chickens \ r2_chickens \ ///
				r1_cattle \ r2_cattle \ ///
				r1_radio \ r2_radio \ ///
				r1_television \ r2_television

matrix regs4 = 	r1_electrification \ r2_electrification \ ///
				r1_population \ r2_population
			
matrix regs = regs1 \ regs2 \ regs3 \ regs4	
		
*mat2txt, matrix(regs) format(%9.3f %9.3f %9.3f %9.3f %9.0f) saving(`appDir'/tabb2b.csv) replace
mat list regs

sureg ///
	(female tx1 tx2 tx3 busia funded market) ///
	(age tx1 tx2 tx3 busia funded market) ///
	(senior tx1 tx2 tx3 busia funded market) ///
	(someschool tx1 tx2 tx3 busia funded market) ///
	(spouse tx1 tx2 tx3 busia funded market) ///
	(notfarmer tx1 tx2 tx3 busia funded market) ///
	(employed tx1 tx2 tx3 busia funded market) ///
	(awareness tx1 tx2 tx3 busia funded market) ///
	(bank tx1 tx2 tx3 busia funded market) ///
	(earn tx1 tx2 tx3 busia funded market) ///
	(hhsize tx1 tx2 tx3 busia funded market) ///
	(hhyouth tx1 tx2 tx3 busia funded market) ///
	(wall tx1 tx2 tx3 busia funded market) ///
	(ownland tx1 tx2 tx3 busia funded market) ///
	(distance tx1 tx2 tx3 busia funded market) ///
	(spend tx1 tx2 tx3 busia funded market) ///
	(bednets tx1 tx2 tx3 busia funded market) ///
	(bikes tx1 tx2 tx3 busia funded market) ///
	(sofapieces tx1 tx2 tx3 busia funded market) ///
	(chickens tx1 tx2 tx3 busia funded market) ///
	(cattle tx1 tx2 tx3 busia funded market) ///
	(radio tx1 tx2 tx3 busia funded market) ///
	(television tx1 tx2 tx3 busia funded market) ///
	(electrification tx1 tx2 tx3 busia funded market) ///
	(population tx1 tx2 tx3 busia funded market)	

test _b[tx1] = _b[tx2] = _b[tx3] = 0

sureg ///
	(hhsize tx1 tx2 tx3 busia funded market) ///
	(wall tx1 tx2 tx3 busia funded market) ///
	(chickens tx1 tx2 tx3 busia funded market) ///
	(age tx1 tx2 tx3 busia funded market) ///
	(someschool tx1 tx2 tx3 busia funded market) ///
	(notfarmer tx1 tx2 tx3 busia funded market) ///
	(bank tx1 tx2 tx3 busia funded market) ///
	(employed tx1 tx2 tx3 busia funded market) ///
	(senior tx1 tx2 tx3 busia funded market) ///	
	(electrification tx1 tx2 tx3 busia funded market) ///
	(population tx1 tx2 tx3 busia funded market)	

test _b[tx1] = _b[tx2] = _b[tx3] = 0
* note that p-value on F-stat is only 0.0656 for the select 9 variables


********************************************************************************
* Table B3: Characteristics of households taking-up electricity by treatment arm
********************************************************************************

use `data'/demand.dta, clear

local varList female age senior someschool spouse notfarmer employed awareness bank wall radio television 
foreach v of local varList {
	gen `v'1=`v'*100
}

global dem "female1 age senior1 someschool1 spouse1 notfarmer1 employed1 awareness1 bank1 earn"
global hh "hhsize hhyouth wall1 ownland distance spendl"
global ass "bednets sofapieces chickens radio1 television1"

estpost su $dem if tx3==1 & takeup==1
est store C
estpost su $dem if tx2==1 & takeup==1
est store B
estpost su $dem if tx1==1 & takeup==1
est store A
estpost su $dem if anytx==0 & takeup==1 & connected==0
est store D
/*
esttab C B A D using `appDir'/tabb3a.tex, replace ///
	varlabels (female1 "\hspace{1.5mm} Female (\%)" age "\hspace{1.5mm} Age (years)" senior1 "\hspace{1.5mm} Senior citizen (\%)" ///
	someschool1 "\hspace{1.5mm} Attended secondary school (\%)" spouse1 "\hspace{1.5mm} Married (\%)" ///
	notfarmer1 "\hspace{1.5mm} Not a farmer (\%)" employed1 "\hspace{1.5mm} Employed (\%)" awareness1 "\hspace{1.5mm} Basic political awareness (\%)" ///
	bank1 "\hspace{1.5mm} Has bank account (\%)" earn "\hspace{1.5mm} Monthly earnings (USD)") ///
	nomtitles ///
	cells(mean(fmt(2))) label booktabs nonum collabels(none) nogaps f noobs plain
*/
estpost su $hh if tx3==1 & takeup==1
est store C
estpost su $hh if tx2==1 & takeup==1
est store B
estpost su $hh if tx1==1 & takeup==1
est store A
estpost su $hh if anytx==0 & takeup==1 & connected==0
est store D
/*
esttab C B A D using `appDir'/tabb3b.tex, replace ///
	varlabels (hhsize "\hspace{1.5mm} Number of members" ///
	hhyouth "\hspace{1.5mm} Youth members" wall1 "\hspace{1.5mm} High-quality walls (\%)" ///
	ownland "\hspace{1.5mm} Land (acres)" distance "\hspace{1.5mm} Distance to transformer (m)" ///
	spendl "\hspace{1.5mm} Monthly (non-charcoal) energy (USD)") ///
	nomtitles ///
	cells(mean(fmt(2))) label booktabs nonum collabels(none) nogaps f noobs plain
*/
estpost su $ass if tx3==1 & takeup==1
est store C
estpost su $ass if tx2==1 & takeup==1
est store B
estpost su $ass if tx1==1 & takeup==1
est store A
estpost su $ass if anytx==0 & takeup==1 & connected==0
est store D
/*
esttab C B A D using `appDir'/tabb3c.tex, replace ///
	varlabels (bednets "\hspace{1.5mm} Bednets" ///
	sofapieces "\hspace{1.5mm} Sofa pieces" chickens "\hspace{1.5mm} Chickens" ///
	radio1 "\hspace{1.5mm} Owns radio (\%)" television1 "\hspace{1.5mm} Owns television (\%)") ///
	nomtitles ///
	cells(mean(fmt(2))) label booktabs nonum collabels(none) nogaps f noobs plain
*/
* Note: Next step is to compare columns 2, 3, and 4 to column 1. For any differences, significance stars are manually recorded into LaTeX.

* Checking for statistically significant differences.
use `data'/demand.dta, clear

drop if connected==1
local varList female age senior someschool spouse notfarmer employed awareness bank earn wall radio television
foreach v of local varList {
	gen `v'1=`v'*100
}

* comparing tx3 to tx2
gen mark=.
replace mark=1 if tx3==1 & takeup==1
replace mark=2 if tx2==1 & takeup==1
estpost ttest female1 age senior1 someschool1 spouse1 notfarmer1 employed1 awareness1 bank1 earn ///
	hhsize hhyouth wall1 ownland distance spendl ///
	bednets sofapieces chickens radio1 television1 if mark!=., by(mark)
estimates store ttest_high_med
drop mark

* comparing tx3 to tx1
gen mark=.
replace mark=1 if tx3==1 & takeup==1
replace mark=2 if tx1==1 & takeup==1
estpost ttest female1 age senior1 someschool1 spouse1 notfarmer1 employed1 awareness1 bank1 earn ///
	hhsize hhyouth wall1 ownland distance spendl ///
	bednets sofapieces chickens radio1 television1 if mark!=., by(mark)
estimates store ttest_high_low
drop mark

* comparing tx3 to control
gen mark=.
replace mark=1 if tx3==1 & takeup==1
replace mark=2 if anytx==0 & takeup==1
estpost ttest female1 age senior1 someschool1 spouse1 notfarmer1 employed1 awareness1 bank1 earn ///
	hhsize hhyouth wall1 ownland distance spendl ///
	bednets sofapieces chickens radio1 television1 if mark!=., by(mark)
estimates store ttest_high_control
drop mark

********************************************************************************
* Table B4: Impact of connection subsidy on take-up
********************************************************************************

*************************************************
* A. Interactions with community-level variables
*************************************************

use `data'/demand.dta, clear
gen takeuptemp=takeup*100
drop takeup
rename takeuptemp takeup

local int1List ///
	busia market funded electrification population wall age someschool employed senior hhsize chickens notfarmer bank black

foreach v in `int1List' {
	gen t1_`v'=tx1*`v'
	gen t2_`v'=tx2*`v'
	gen t3_`v'=tx3*`v'
}

gen interaction1=.
gen interaction2=.
gen interaction3=.
gen variable=.

* Regression 1
reg takeup tx1 tx2 tx3 ///
	hhsize wall age someschool senior chickens notfarmer bank employed ///
	busia funded market electrification population ///
	if connected==0, vce(cluster siteno)
	qui estadd ysumm
	qui estadd local controls "Yes", replace
	qui estimates store m1

* Regression 2
replace interaction1=t1_busia
replace interaction2=t2_busia
replace interaction3=t3_busia
replace variable=busia

reg takeup tx1 tx2 tx3 ///
	interaction1 interaction2 interaction3 ///
	hhsize wall age someschool senior chickens notfarmer bank employed ///
	variable funded market electrification population ///
	if connected==0, vce(cluster siteno)
	qui estadd ysumm
	qui estadd local controls "Yes", replace
	qui estimates store m2

* Regression 3
replace interaction1=t1_funded
replace interaction2=t2_funded
replace interaction3=t3_funded
replace variable=funded

reg takeup tx1 tx2 tx3 ///
	interaction1 interaction2 interaction3 ///
	hhsize wall age someschool senior chickens notfarmer bank employed ///
	busia variable market electrification population ///
	if connected==0, vce(cluster siteno)
	qui estadd ysumm
	qui estadd local controls "Yes", replace
	qui estimates store m3
	
* Regression 4
replace interaction1=t1_market
replace interaction2=t2_market
replace interaction3=t3_market
replace variable=market

reg takeup tx1 tx2 tx3 ///
	interaction1 interaction2 interaction3 ///
	hhsize wall age someschool senior chickens notfarmer bank employed  ///
	busia funded variable electrification population ///
	if connected==0, vce(cluster siteno)
	qui estadd ysumm
	qui estadd local controls "Yes", replace
	qui estimates store m4	
	
* Regression 5
replace interaction1=t1_population
replace interaction2=t2_population
replace interaction3=t3_population
replace variable=population

reg takeup tx1 tx2 tx3 ///
	interaction1 interaction2 interaction3 ///
	hhsize wall age someschool senior chickens notfarmer bank employed ///
	busia funded market variable electrification ///
	if connected==0, vce(cluster siteno)
	qui estadd ysumm
	qui estadd local controls "Yes", replace
	qui estimates store m5
/*
esttab m1 m2 m3 m4 m5 using `appDir'/tabb4a.tex, ///
	replace b(%10.1f) se scalars( ///
	"N Observations" "r2 R-squared") ///
	sfmt(%6.0g %10.2f %10.0fc %10.2f) ///
	star(* 0.10 ** 0.05 *** 0.01) noobs label nomtitles ///
	align(cccccc) width(\hsize) nonotes ///
	keep(tx1 tx2 tx3 variable interaction1 interaction2 interaction3) ///
	order(_cons tx1 tx2 tx3 variable interaction1 interaction2 interaction3) ///
	varlabels(_cons "Control (intercept)" tx1 "T1: Low subsidy\textemdash 29\% discount" ///
	tx2 "T2: Medium subsidy\textemdash 57\% discount" tx3 "T3: High subsidy\textemdash 100\% discount" ///
	variable "Interacted variable" ///
	interaction1 "T1 $\times$ interacted variable" ///
	interaction2 "T2 $\times$ interacted variable" ///
	interaction3 "T3 $\times$ interacted variable") ///
	title(\textemdash Impact of connection subsidy on take-up: Interactions with community-level variables}\vspace{-0.3em)
*/
* Note: Latex file is manually edited for formatting purposes

*************************************************
* B. Interactions with household-level variables
*************************************************

use `data'/demand.dta, clear
gen takeuptemp=takeup*100
drop takeup
rename takeuptemp takeup

local int1List ///
	busia market funded electrification population wall age someschool employed senior hhsize chickens notfarmer bank

foreach v in `int1List' {
	gen t1_`v'=tx1*`v'
	gen t2_`v'=tx2*`v'
	gen t3_`v'=tx3*`v'
}

gen interaction1=.
gen interaction2=.
gen interaction3=.
gen variable=.

* Regression 2
replace interaction1=t1_hhsize
replace interaction2=t2_hhsize
replace interaction3=t3_hhsize
replace variable=hhsize

reg takeup tx1 tx2 tx3 ///
	interaction1 interaction2 interaction3 ///
	variable wall age someschool senior chickens notfarmer bank employed ///
	busia funded market electrification population ///
	if connected==0, vce(cluster siteno)
	qui estadd ysumm
	qui estadd local controls "Yes", replace
	qui estimates store m2

* Regression 3
replace interaction1=t1_age
replace interaction2=t2_age
replace interaction3=t3_age
replace variable=age

reg takeup tx1 tx2 tx3 ///
	interaction1 interaction2 interaction3 ///
	hhsize wall variable someschool senior chickens notfarmer bank employed ///
	busia funded market electrification population ///
	if connected==0, vce(cluster siteno)
	qui estadd ysumm
	qui estadd local controls "Yes", replace
	qui estimates store m3
	
* Regression 4
replace interaction1=t1_senior
replace interaction2=t2_senior
replace interaction3=t3_senior
replace variable=senior

reg takeup tx1 tx2 tx3 ///
	interaction1 interaction2 interaction3 ///
	hhsize wall age someschool variable chickens notfarmer bank employed  ///
	busia funded market electrification population ///
	if connected==0, vce(cluster siteno)
	qui estadd ysumm
	qui estadd local controls "Yes", replace
	qui estimates store m4

* Regression 5
replace interaction1=t1_chickens
replace interaction2=t2_chickens
replace interaction3=t3_chickens
replace variable=chickens

reg takeup tx1 tx2 tx3 ///
	interaction1 interaction2 interaction3 ///
	hhsize wall age someschool senior variable notfarmer bank employed ///
	busia funded market electrification population ///
	if connected==0, vce(cluster siteno)
	qui estadd ysumm
	qui estadd local controls "Yes", replace
	qui estimates store m5

* Regression 6
replace interaction1=t1_bank
replace interaction2=t2_bank
replace interaction3=t3_bank
replace variable=employed

reg takeup tx1 tx2 tx3 ///
	interaction1 interaction2 interaction3 ///
	hhsize wall age someschool senior chickens notfarmer variable employed  ///
	busia funded market electrification population ///
	if connected==0, vce(cluster siteno)
	qui estadd ysumm
	qui estadd local controls "Yes", replace
	qui estimates store m6

* Regression 7
replace interaction1=t1_notfarmer
replace interaction2=t2_notfarmer
replace interaction3=t3_notfarmer
replace variable=notfarmer

reg takeup tx1 tx2 tx3 ///
	interaction1 interaction2 interaction3 ///
	hhsize wall age someschool senior chickens variable bank employed ///
	busia funded market electrification population ///
	if connected==0, vce(cluster siteno)
	qui estadd ysumm
	qui estadd local controls "Yes", replace
	qui estimates store m7
/*
esttab m2 m3 m4 using `appDir'/tabb4b.tex, ///
	replace b(%10.1f) se scalars( ///
	"N Observations" "r2 R-squared") ///
	sfmt(%6.0g %10.2f %10.0fc %10.2f) ///
	star(* 0.10 ** 0.05 *** 0.01) noobs label nomtitles ///
	align(cccccc) width(\hsize) nonotes ///
	keep(tx1 tx2 tx3 variable interaction1 interaction2 interaction3) ///
	order(_cons tx1 tx2 tx3 variable interaction1 interaction2 interaction3) ///
	varlabels(_cons "Control (intercept)" tx1 "T1: Low subsidy\textemdash 29\% discount" ///
	tx2 "T2: Medium subsidy\textemdash 57\% discount" tx3 "T3: High subsidy\textemdash 100\% discount" ///
	variable "Interacted variable" ///
	interaction1 "T1 $\times$ interacted variable" ///
	interaction2 "T2 $\times$ interacted variable" ///
	interaction3 "T3 $\times$ interacted variable") ///
	title(\textemdash Impact of connection subsidy on take-up: Interactions with household-level variables}\vspace{-0.3em)
*/
* Note: Latex file is manually edited for formatting purposes

*************************************************
* C. Interactions with household-level variables
*************************************************

/*
esttab m5 m6 m7 using `appDir'/tabb4c.tex, ///
	replace b(%10.1f) se scalars( ///
	"N Observations" "r2 R-squared") ///
	sfmt(%6.0g %10.2f %10.0fc %10.2f) ///
	star(* 0.10 ** 0.05 *** 0.01) noobs label nomtitles ///
	align(cccccc) width(\hsize) nonotes ///
	keep(tx1 tx2 tx3 variable interaction1 interaction2 interaction3) ///
	order(_cons tx1 tx2 tx3 variable interaction1 interaction2 interaction3) ///
	varlabels(_cons "Control (intercept)" tx1 "T1: Low subsidy\textemdash 29\% discount" ///
	tx2 "T2: Medium subsidy\textemdash 57\% discount" tx3 "T3: High subsidy\textemdash 100\% discount" ///
	variable "Interacted variable" ///
	interaction1 "T1 $\times$ interacted variable" ///
	interaction2 "T2 $\times$ interacted variable" ///
	interaction3 "T3 $\times$ interacted variable") ///
	title(\textemdash Impact of connection subsidy on take-up: Interactions with household-level variables}\vspace{-0.3em)
*/
* Note: Latex file is manually edited for formatting purposes

*************************************************
* D. All controls
*************************************************

use `data'/demand.dta, clear
gen takeuptemp=takeup*100
drop takeup
rename takeuptemp takeup

reg takeup tx1 tx2 tx3 ///
	if connected==0, vce(cluster siteno)
	qui estimates store m1

reg takeup tx1 tx2 tx3 ///
	hhsize wall age senior someschool chickens notfarmer bank employed ///
	busia funded market electrification population ///
	if connected==0, vce(cluster siteno)
	qui estimates store m2

replace distance=distance/100
reg takeup tx1 tx2 tx3 ///
	female age senior someschool spouse notfarmer awareness bank employed earn ///
	hhsize hhyouth wall ownland distance spendl ///
	bednets bikes sofapieces chickens cattle radio television ///
	busia funded market electrification population ///
	if connected==0, vce(cluster siteno)
	qui estimates store m8
/*
esttab m1 m2 m8 using `appDir'/tabb4d.tex, ///
	replace b(%10.1f) se scalars( ///
	"N Observations" "r2 R-squared") ///
	sfmt(%6.0g %10.2f %10.0fc %10.2f) ///
	star(* 0.10 ** 0.05 *** 0.01) noobs label nomtitles ///
	align(cccccc) width(\hsize) nonotes ///
	keep( ///
	_cons tx1 tx2 tx3 ///
	female age someschool spouse notfarmer awareness bank earn ///
	senior employed ///
	hhsize hhyouth wall ownland distance spendl ///
	bednets bikes sofapieces chickens cattle radio television ///
	electrification population ///
	busia funded market ///
	) ///
	order( ///
	_cons tx1 tx2 tx3 ///
	female age senior someschool spouse notfarmer employed awareness bank earn ///
	hhsize hhyouth wall ownland distance spendl ///
	bednets bikes sofapieces chickens cattle radio television ///
	electrification population ///
	busia funded market ///
	) ///
	varlabels(_cons "Control (intercept)" tx1 "T1: Low subsidy\textemdash 29\% discount" ///
	tx2 "T2: Medium subsidy\textemdash 57\% discount" tx3 "T3: High subsidy\textemdash 100\% discount" ///
	female "Female=1" age "Age (years)" someschool "Attended secondary school=1" ///
	spouse "Married=1" notfarmer "Not a farmer=1" awareness "Basic political awareness=1" ///
	bank "Has bank account=1" earn "Monthly earnings (USD)" senior "Senior citizen=1" ///
	employed "Employed=1" hhsize "Number of members" ///
	hhyouth "Youth members (age $\leq$ 18)" ///
	wall "High-quality walls=1" ownland "Land (acres)" ///
	distance "Distance to transformer (m)" spendl "Monthly (non-charcoal) energy (USD)" ///
	bednets "Number of bednets" bikes "Number of bicycles" /// 
	sofapieces "Number of sofa pieces" chickens "Number of chickens" ///
	cattle "Number of cattle" radio "Owns radio=1" ///
	television "Owns television=1" electrification "Electrification rate (%)" ///
	population "Community population" busia "Busia=1" ///
	funded "Funded and installed early on=1" market "Market status=1") ///
	title(\textemdash Impact of connection subsidy on take-up: All controls}\vspace{-0.3em)
*/
*************************************************
* E. Price
*************************************************

use `data'/demand.dta, clear

gen takeuptemp=takeup*100
drop takeup
rename takeuptemp takeup

gen p1=price/`rate'
gen p2=p1^2 /1000

reg takeup p1 p2 ///
	if connected==0, vce(cluster siteno)
	qui estimates store m8

reg takeup p1 p2 ///
	hhsize wall age someschool senior chickens notfarmer bank employed  ///
	busia funded market electrification population ///
	if connected==0, vce(cluster siteno)
	qui estimates store m9
/*
esttab m8 m9 using `appDir'/tabb4e.tex, ///
	replace b(%10.1f) se scalars( ///
	"N Observations" "r2 R-squared") ///
	sfmt(%6.0g %10.2f) ///
	star(* 0.10 ** 0.05 *** 0.01) noobs label nomtitles ///
	align(cccccc) width(\hsize) nonotes ///
	keep( ///
	p1 p2 ///
	hhsize wall age someschool senior chickens notfarmer bank employed  ///
	busia funded market electrification population ///
	) ///
	order( ///
	p1 p2 ///
	age senior someschool notfarmer bank employed hhsize wall chickens ///
	busia funded market electrification population ///
	) ///
	varlabels( ///
	p1 "Price" ///
	p2 "Price$^2$ $\times$ 1000" ///
	hhsize "Number of members" ///
	wall "High-quality walls=1" ///
	age "Age (years)" ///
	someschool "Attended secondary school=1" ///
	senior "Senior citizen=1" ///
	chickens "Number of chickens" ///
	notfarmer "Not a farmer=1" ///
	bank "Has bank account=1" ///
	employed  "Employed=1" ///
	busia "Busia=1" ///
	funded "Funded early on=1" ///
	market "Market status=1" ///
	electrification "Electrification rate (\%)" ///
	population "Community population" ///
	) ///
	title(\textemdash Impact of price on take-up}\vspace{-0.3em)
*/
* Note: Latex file is manually edited for formatting purposes


********************************************************************************
* Table B5: Actual vs fitted total cost and ATC values
********************************************************************************

* Panel B
use `data'/costs.dta, clear
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

local n=2
	gen x_`n'=2.1/100*ucompounds
	gen co_`n'=`b0'/x_`n'+`b1'+`b2'*x_`n'
	gen w_co_`n'=co_`n'*wa
	egen cuco_`n'=sum(w_co_`n')
	drop x_`n' co_`n' w_co_`n'

local n=5
	gen x_`n'=4.8/100*ucompounds
	gen co_`n'=`b0'/x_`n'+`b1'+`b2'*x_`n'
	gen w_co_`n'=co_`n'*wa
	egen cuco_`n'=sum(w_co_`n')
	drop x_`n' co_`n' w_co_`n'

local n=17
	gen x_`n'=17.1/100*ucompounds
	gen co_`n'=`b0'/x_`n'+`b1'+`b2'*x_`n'
	gen w_co_`n'=co_`n'*wa
	egen cuco_`n'=sum(w_co_`n')
	drop x_`n' co_`n' w_co_`n'
	
foreach n in 25 50 75 100 {
	gen x_`n'=`n'/100*ucompounds
	gen co_`n'=`b0'/x_`n'+`b1'+`b2'*x_`n'
	gen w_co_`n'=co_`n'*wa
	egen cuco_`n'=sum(w_co_`n')
	drop x_`n' co_`n' w_co_`n'
}

duplicates drop cuco_2, force
keep cuco_*
gen id=1

reshape long cuco_, i(id)
drop id
rename _j proportion
rename cuco_ atc
keep atc proportion
replace proportion = 2.1 if proportion==2
replace proportion = 4.8 if proportion==5
replace proportion = 17.1 if proportion==17
gen tc=atc*proportion/100*84.7

* Panel C
use `data'/costs.dta, clear

gen atc_usd=cost_b1
drop if atc_usd==.

gen m1=connect_d
gen m2=m1^2

reg atc_usd m1 m2, vce(robust)
mat b=e(b)
matlist b
local a=b[1,3]
local b1=b[1,1]
local b2=b[1,2]

egen totc=sum(compounds)
gen wa=compounds/totc
drop totc

local n=2
	gen x_`n'=2.1/100*ucompounds
	gen co_`n'=`a'+`b1'*x_`n'+`b2'*x_`n'^2	
	gen w_co_`n'=co_`n'*wa
	egen cuco_`n'=sum(w_co_`n')
	drop x_`n' co_`n' w_co_`n'

local n=5
	gen x_`n'=4.8/100*ucompounds
	gen co_`n'=`a'+`b1'*x_`n'+`b2'*x_`n'^2	
	gen w_co_`n'=co_`n'*wa
	egen cuco_`n'=sum(w_co_`n')
	drop x_`n' co_`n' w_co_`n'

local n=17
	gen x_`n'=17.1/100*ucompounds
	gen co_`n'=`a'+`b1'*x_`n'+`b2'*x_`n'^2	
	gen w_co_`n'=co_`n'*wa
	egen cuco_`n'=sum(w_co_`n')
	drop x_`n' co_`n' w_co_`n'
	
foreach n in 25 50 75 100 {
	gen x_`n'=`n'/100*ucompounds
	gen co_`n'=`a'+`b1'*x_`n'+`b2'*x_`n'^2	
	gen w_co_`n'=co_`n'*wa
	egen cuco_`n'=sum(w_co_`n')
	drop x_`n' co_`n' w_co_`n'
}

duplicates drop cuco_2, force
keep cuco_*
gen id=1

reshape long cuco_, i(id)
drop id
rename _j proportion
rename cuco_ atc
keep atc proportion
replace proportion = 2.1 if proportion==2
replace proportion = 4.8 if proportion==5
replace proportion = 17.1 if proportion==17
gen tc=atc*proportion/100*84.7

********************************************************************************
* Tables B6A, B6B, B6C: See "3_impacts.do"
********************************************************************************

********************************************************************************
* Table B8A: Impact of randomized offers on hypothetical and actual take-up
********************************************************************************

use `data'/demand.dta, clear
gen wtp1=WTP_r1*100
gen wtp2=WTP_r2*100
gen takeuptemp=takeup*100
drop takeup
rename takeuptemp takeup
gen cvsample=0

replace cvsample=1 if WTP_amt==75000 | WTP_amt==35000 | WTP_amt==25000 | WTP_amt==20000 | WTP_amt==15000 | WTP_amt==10000 | WTP_amt==0

foreach v in 0 10000 15000 20000 25000 35000 75000 {
	gen wtp_`v'=0
	replace wtp_`v'=1 if WTP_amt==`v'
	}

local int1List ///
	busia market funded electrification population wall someschool employed senior hhsize chickens notfarmer bank

foreach v in `int1List' {
	gen wtp_0_`v'=wtp_0*`v'
	gen wtp_10000_`v'=wtp_10000*`v'
	gen wtp_15000_`v'=wtp_15000*`v'
	gen wtp_20000_`v'=wtp_20000*`v'
	gen wtp_25000_`v'=wtp_25000*`v'
	gen wtp_75000_`v'=wtp_75000*`v'
}

* Regression 1
reg wtp1 wtp_0 wtp_10000 wtp_15000 wtp_20000 wtp_25000 wtp_75000 ///
	hhsize wall age someschool senior chickens notfarmer bank employed ///
	busia funded market electrification population ///
	if cvsample==1, vce(cluster siteno)
	qui estadd ysumm
	qui estadd local controls "Yes", replace
	qui estimates store m1

* Regression 2
reg wtp2 wtp_0 wtp_10000 wtp_15000 wtp_20000 wtp_25000 wtp_75000 ///
	hhsize wall age someschool senior chickens notfarmer bank employed ///
	busia funded market electrification population ///
	if cvsample==1, vce(cluster siteno)
	qui estadd ysumm
	qui estadd local controls "Yes", replace
	qui estimates store m2

* Regression 3
replace wtp_25000=tx1
replace wtp_15000=tx2
replace wtp_0=tx3

reg takeup wtp_0 wtp_15000 wtp_25000 ///
	hhsize wall age someschool senior chickens notfarmer bank employed ///
	busia funded market electrification population ///
	if connected==0, vce(cluster siteno)
	qui estadd ysumm
	qui estadd local controls "Yes", replace
	qui estimates store m3
/*
esttab m1 m2 m3 using `appDir'/tabb8a.tex, ///
	replace b(%10.1f) se scalars( ///
	"ymean Mean of dependent variable" ///
	"N Observations" "r2 R$^2$") ///
	sfmt(%10.1f %10.0fc %10.2f) ///
	star(* 0.10 ** 0.05 *** 0.01) noobs label nomtitles ///
	align(cccccc) width(\hsize) nonotes ///
	keep(wtp_0 wtp_10000 wtp_15000 wtp_20000 wtp_25000 wtp_75000 ///
	age senior someschool notfarmer employed bank hhsize wall chickens) ///
	order(wtp_75000 wtp_25000 wtp_20000 wtp_15000 wtp_10000 wtp_0 ///
	age senior someschool notfarmer employed bank hhsize wall chickens) ///
	varlabels( ///
	wtp_75000 "\\$853 offer" ///
	wtp_25000 "\\$284 offer / T1: Low subsidy\textemdash 29\% discount" ///
	wtp_20000 "\\$227 offer" ///
	wtp_15000 "\\$171 offer / T2: Medium subsidy\textemdash 57\% discount" ///
	wtp_10000 "\\$114 offer" ///
	wtp_0 "Free offer / T3: High subsidy\textemdash 100\% discount" ///
	hhsize "Number of household members" ///
	wall "High-quality walls=1" ///
	chickens "Number of chickens=1" ///
	age "Age (years)" ///
	someschool "Attended secondary school=1" ///
	senior "Senior citizen=1" ///
	notfarmer "Not a farmer=1" ///
	bank "Has bank account=1" ///
	employed "Employed=1") ///
	title(\textemdash Impact of WTP offer on stated take-up of electricity connections}\vspace{-0.3em)
*/
* Note: Latex file is manually edited for formatting purposes

*************************************************************************************************
* Table B8B: Impact of randomized offers on take-up, with interactions
*************************************************************************************************

use `data'/demand.dta, clear
gen wtp1=WTP_r1*100
gen wtp2=WTP_r2*100
gen takeuptemp=takeup*100
drop takeup
rename takeuptemp takeup
gen cvsample=0

replace cvsample=1 if WTP_amt==75000 | WTP_amt==35000 | WTP_amt==25000 | WTP_amt==20000 | WTP_amt==15000 | WTP_amt==10000 | WTP_amt==0

foreach v in 0 10000 15000 20000 25000 35000 75000 {
	gen wtp_`v'=0
	replace wtp_`v'=1 if WTP_amt==`v'
	}

local int1List ///
	busia market funded electrification population wall someschool employed senior hhsize chickens notfarmer bank

foreach v in `int1List' {
	gen wtp_0_`v'=wtp_0*`v'
	gen wtp_10000_`v'=wtp_10000*`v'
	gen wtp_15000_`v'=wtp_15000*`v'
	gen wtp_20000_`v'=wtp_20000*`v'
	gen wtp_25000_`v'=wtp_25000*`v'
	gen wtp_75000_`v'=wtp_75000*`v'
}

gen interaction1=.
gen interaction2=.
gen interaction3=.
gen interaction4=.
gen interaction5=.
gen interaction6=.
gen variable=.

* Regression 4
reg wtp2 wtp_0 wtp_10000 wtp_15000 wtp_20000 wtp_25000 wtp_75000 ///
	hhsize wall age someschool senior chickens notfarmer bank employed ///
	busia funded market electrification population ///
	if cvsample==1, vce(cluster siteno)
	qui estadd ysumm
	qui estadd local controls "Yes", replace
	qui estimates store m4

* Regression 5
replace interaction1=wtp_75000_wall
replace interaction2=wtp_25000_wall
replace interaction3=wtp_20000_wall
replace interaction4=wtp_15000_wall
replace interaction5=wtp_10000_wall
replace interaction6=wtp_0_wall
replace variable=wall

reg wtp2 wtp_0 wtp_10000 wtp_15000 wtp_20000 wtp_25000 wtp_75000 ///
	variable interaction1 interaction2 interaction3 interaction4 interaction5 interaction6 ///
	hhsize age someschool senior chickens notfarmer bank employed ///
	busia funded market electrification population ///
	if cvsample==1, vce(cluster siteno)
	qui estadd ysumm
	qui estadd local controls "Yes", replace
	qui estimates store m5

* Regression 6
replace interaction1=wtp_75000_bank
replace interaction2=wtp_25000_bank
replace interaction3=wtp_20000_bank
replace interaction4=wtp_15000_bank
replace interaction5=wtp_10000_bank
replace interaction6=wtp_0_bank
replace variable=bank

reg wtp2 wtp_0 wtp_10000 wtp_15000 wtp_20000 wtp_25000 wtp_75000 ///
	variable interaction1 interaction2 interaction3 interaction4 interaction5 interaction6 ///
	hhsize wall age someschool senior chickens notfarmer employed ///
	busia funded market electrification population ///
	if cvsample==1, vce(cluster siteno)
	qui estadd ysumm
	qui estadd local controls "Yes", replace
	qui estimates store m6

* Regression 7
replace interaction1=wtp_75000_someschool
replace interaction2=wtp_25000_someschool
replace interaction3=wtp_20000_someschool
replace interaction4=wtp_15000_someschool
replace interaction5=wtp_10000_someschool
replace interaction6=wtp_0_someschool
replace variable=someschool

reg wtp2 wtp_0 wtp_10000 wtp_15000 wtp_20000 wtp_25000 wtp_75000 ///
	variable interaction1 interaction2 interaction3 interaction4 interaction5 interaction6 ///
	hhsize wall age senior chickens notfarmer bank employed ///
	busia funded market electrification population ///
	if cvsample==1, vce(cluster siteno)
	qui estadd ysumm
	qui estadd local controls "Yes", replace
	qui estimates store m7
/*
esttab m4 m5 m6 m7 using `appDir'/tabb8b.tex, ///
	replace b(%10.1f) se scalars( ///
	"ymean Mean of dependent variable" ///
	"N Observations" "r2 R$^2$") ///
	sfmt(%10.1f %10.0fc %10.2f) ///
	star(* 0.10 ** 0.05 *** 0.01) noobs label nomtitles ///
	align(cccccc) width(\hsize) nonotes ///
	keep(wtp_0 wtp_10000 wtp_15000 wtp_20000 wtp_25000 wtp_75000 variable interaction6 interaction5 interaction4 interaction3 interaction2 interaction1) ///
	order(wtp_75000 wtp_25000 wtp_20000 wtp_15000 wtp_10000 wtp_0 variable interaction1 interaction2 interaction3 interaction4 interaction5 interaction6) ///
	varlabels( ///
	wtp_75000 "\\$853 offer" ///
	wtp_25000 "\\$284 offer / T1: Low subsidy\textemdash 29\% discount" ///
	wtp_20000 "\\$227 offer" ///
	wtp_15000 "\\$171 offer / T2: Medium subsidy\textemdash 57\% discount" ///
	wtp_10000 "\\$114 offer" ///
	wtp_0 "Free offer / T3: High subsidy\textemdash 100\% discount" ///
	variable "Interacted variable" ///
	interaction1 "\\$853 offer $\times$ interacted variable" ///
	interaction2 "\\$284 offer $\times$ interacted variable" ///
	interaction3 "\\$227 offer $\times$ interacted variable" ///
	interaction4 "\\$171 offer $\times$ interacted variable" ///
	interaction5 "\\$114 offer $\times$ interacted variable" ///
	interaction6 "Free offer $\times$ interacted variable") ///
	title(\textemdash Impact of WTP offer on stated take-up of electricity connections}\vspace{-0.3em)
* Note: Latex file is manually edited for formatting purposes
*/
*************************************************************************************************
* Table B8C: Predictors of credit constraints
*************************************************************************************************

use `data'/demand.dta, clear

gen switch=.
replace switch=100 if WTP_r1==1 & WTP_r2==0
replace switch=0 if WTP_r1==1 & WTP_r2==1
* replace switch=0 if WTP_r1==0 & WTP_r2==0
keep if switch==100 | switch==0
* This leaves us with 2,277 / 2,289 observations. 12 observations with missing WTP data.

gen cvsample=0
replace cvsample=1 if WTP_amt==75000 | WTP_amt==35000 | WTP_amt==25000 | WTP_amt==20000 | WTP_amt==15000 | WTP_amt==10000 | WTP_amt==0

foreach v in 0 10000 15000 20000 25000 35000 75000 {
	gen wtp_`v'=0
	replace wtp_`v'=1 if WTP_amt==`v'
	}

* Regression 1
reg switch wtp_10000 wtp_15000 wtp_20000 wtp_25000 wtp_35000 wtp_75000 ///
	if cvsample==1, vce(cluster siteno)
	qui estadd ysumm
	qui estadd local controls "Yes", replace
	qui estimates store m1

* Regression 2
reg switch wtp_10000 wtp_15000 wtp_20000 wtp_25000 wtp_35000 wtp_75000 ///
	hhsize wall age someschool senior chickens notfarmer bank employed ///
	if cvsample==1, vce(cluster siteno)
	qui estadd ysumm
	qui estadd local controls "Yes", replace
	qui estimates store m2
/*	
esttab m1 m2 using `appDir'/tabb8c.tex, ///
	replace b(%10.1f) se scalars( ///
	"ymean Mean of dependent variable" ///
	"N Observations" "r2 R$^2$") ///
	sfmt(%10.1f %10.0fc %10.2f) ///
	star(* 0.10 ** 0.05 *** 0.01) noobs label nomtitles ///
	align(cccccc) width(\hsize) nonotes ///
	keep(wtp_10000 wtp_15000 wtp_20000 wtp_25000 wtp_35000 wtp_75000 ///
	age senior someschool notfarmer employed bank hhsize wall chickens) ///
	order(wtp_75000 wtp_35000 wtp_25000 wtp_20000 wtp_15000 wtp_10000 ///
	age senior someschool notfarmer employed bank hhsize wall chickens) ///
	varlabels( ///
	wtp_75000 "\\$853 offer" ///
	wtp_35000 "\\$398 offer / Existing fixed price" ///
	wtp_25000 "\\$284 offer / T1: Low subsidy\textemdash 29\% discount" ///
	wtp_20000 "\\$227 offer" ///
	wtp_15000 "\\$171 offer / T2: Medium subsidy\textemdash 57\% discount" ///
	wtp_10000 "\\$114 offer" ///
	wtp_0 "Free offer / T3: High subsidy\textemdash 100\% discount" ///
	hhsize "Number of household members" ///
	wall "High-quality walls=1" ///
	chickens "Number of chickens=1" ///
	age "Age (years)" ///
	someschool "Attended secondary school=1" ///
	senior "Senior citizen=1" ///
	notfarmer "Not a farmer=1" ///
	bank "Has bank account=1" ///
	employed "Employed=1") ///
	title(\textemdash Predictors of credit constraints}\vspace{-0.3em)
* Note: Latex file is manually edited for formatting purposes
*/


