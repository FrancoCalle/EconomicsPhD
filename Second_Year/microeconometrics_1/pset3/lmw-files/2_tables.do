clear all
set more off
local rootDir = "/Users/Research/LMW2019"

local data "`rootDir'/data/stata"
local temp "`rootDir'/data/temp"
local tabDir "`rootDir'/tables"

local rate=87.94

********************************************************************************
* Tables A
* Last updated: 23 Apr 2019
********************************************************************************

********************************************************************************
* Table 1: Differences between unconnected and grid connected households at baseline
********************************************************************************

use `data'/demand.dta, clear

gen connected1=1
replace connected1=2 if connected==1

local varList female senior someschool spouse notfarmer employed awareness bank wall radio television
foreach v of local varList {
	gen `v'1=`v'*100
}

estpost ttest female1 age senior1 someschool1 spouse1 notfarmer1 employed1 awareness1 bank1 earn ///
	hhsize hhyouth wall1 ownland distance spendl ///
	bednets sofapieces chickens radio1 television1, by(connected1)
estimates store tab1
/*
esttab tab1 using `tabDir'/tab1.tex, ///
	cells("mu_1(fmt(1) label(Unconnected)) mu_2(fmt(1) label(Connected)) p(fmt(2) label(p-value of diff.))") ///
	replace noobs label ///
	varlabels (female1 "\hspace{2.5mm} Female (\%)" age "\hspace{2.5mm} Age (years)" senior1 "\hspace{2.5mm} Senior citizen (\%)" ///
	someschool1 "\hspace{2.5mm} Attended secondary schooling (\%)" spouse1 "\hspace{2.5mm} Married (\%)" ///
	notfarmer1 "\hspace{2.5mm} Not a farmer (\%)" employed1 "\hspace{2.5mm} Employed (\%)" awareness1 "\hspace{2.5mm} Basic political awareness (\%)" ///
	bank1 "\hspace{2.5mm} Has bank account (\%)" earn "\hspace{2.5mm} Monthly earnings (USD)" hhsize "\hspace{2.5mm} Number of members" ///
	hhyouth "\hspace{2.5mm} Youth members (age $\leq$ 18)" wall1 "\hspace{2.5mm} High-quality walls (\%)" ///
	ownland "\hspace{2.5mm} Land (acres)" distance "\hspace{2.5mm} Distance to transformer (m)" ///
	spendl "\hspace{2.5mm} Monthly (non-charcoal) energy (USD)" bednets "\hspace{2.5mm} Bednets" ///
	sofapieces "\hspace{2.5mm} Sofa pieces" chickens "\hspace{2.5mm} Chickens" ///
	radio1 "\hspace{2.5mm} Owns radio (\%)" television1 "\hspace{2.5mm} Owns television (\%)") ///
	nonumbers nostar alignment(ccc) gaps width(\hsize) ///
	title(Differences between electricity grid unconnected vs. connected households at baseline)

* Note: Latex file is manually edited for formatting purposes
*/
********************************************************************************
* Table 2: Impact of grid connection subsidy on take-up of electricity connections
********************************************************************************

use `data'/demand.dta, clear
gen takeuptemp=takeup*100
drop takeup
rename takeuptemp takeup

local m=200
replace prop`m'm=0 if prop`m'm==.
gen prop`m'm2=prop`m'm*100

local int1List ///
	busia market funded electrification population wall someschool employed senior hhsize chickens notfarmer bank earn black prop`m'm2

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
	if connected==0, vce(cluster siteno)
	qui estadd ysumm
	qui estadd local controls "No", replace
	qui estimates store m1

* Regression 2
reg takeup tx1 tx2 tx3 ///
	hhsize wall age someschool senior chickens notfarmer bank employed ///
	busia funded market electrification population ///
	if connected==0, vce(cluster siteno)
	qui estadd ysumm
	qui estadd local controls "Yes", replace
	qui estimates store m2

* Regression 3
replace interaction1=t1_wall
replace interaction2=t2_wall
replace interaction3=t3_wall
replace variable=wall

reg takeup tx1 tx2 tx3 ///
	interaction1 interaction2 interaction3 ///
	hhsize variable age someschool senior chickens notfarmer bank employed ///
	busia funded market electrification population ///
	if connected==0, vce(cluster siteno)
	qui estadd ysumm
	qui estadd local controls "Yes", replace
	qui estimates store m3

* Regression 4
replace interaction1=t1_someschool
replace interaction2=t2_someschool
replace interaction3=t3_someschool
replace variable=someschool

reg takeup tx1 tx2 tx3 ///
	interaction1 interaction2 interaction3 ///
	hhsize wall age variable senior chickens notfarmer bank employed ///
	busia funded market electrification population ///
	if connected==0, vce(cluster siteno)
	qui estadd ysumm
	qui estadd local controls "Yes", replace
	qui estimates store m4
	
* Note that interaction with farmer moved to appendix

* Regression 5
replace interaction1=t1_electrification
replace interaction2=t2_electrification
replace interaction3=t3_electrification
replace variable=electrification

reg takeup tx1 tx2 tx3 ///
	interaction1 interaction2 interaction3 ///
	hhsize wall age someschool senior chickens notfarmer bank employed ///
	busia funded market variable population ///
	if connected==0, vce(cluster siteno)
	qui estadd ysumm
	qui estadd local controls "Yes", replace
	qui estimates store m5

* Regression 6
replace interaction1=t1_prop`m'm2
replace interaction2=t2_prop`m'm2
replace interaction3=t3_prop`m'm2
replace variable=prop`m'm2

reg takeup tx1 tx2 tx3 ///
	interaction1 interaction2 interaction3 ///
	hhsize wall age someschool senior chickens notfarmer bank employed ///
	busia funded market electrification population variable ///
	if connected==0, vce(cluster siteno)
	qui estadd ysumm
	qui estadd local controls "Yes", replace
	qui estimates store m6

* Regression 7
replace interaction1=t1_black
replace interaction2=t2_black
replace interaction3=t3_black
replace variable=black

reg takeup tx1 tx2 tx3 ///
	interaction1 interaction2 interaction3 ///
	hhsize wall age someschool senior chickens notfarmer bank employed ///
	busia funded market electrification population variable ///
	if connected==0, vce(cluster siteno)
	qui estadd ysumm
	qui estadd local controls "Yes", replace
	qui estimates store m7

* Regression 8
replace interaction1=t1_earn
replace interaction2=t2_earn
replace interaction3=t3_earn
replace variable=earn

reg takeup tx1 tx2 tx3 ///
	interaction1 interaction2 interaction3 ///
	hhsize wall age someschool senior chickens notfarmer bank employed ///
	busia funded market variable population ///
	if connected==0, vce(cluster siteno)
	qui estadd ysumm
	qui estadd local controls "Yes", replace
	qui estimates store m8
/*
esttab m1 m2 m3 m8 m4 m5 m6 m7 using `tabDir'/tab2.tex, ///
	replace b(%10.1f) se scalars( ///
	"N Observations" "r2 R$^2$") ///
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
	title(\textemdash Impact of grid connection subsidy on take-up of electricity connections}\vspace{-0.3em)

* Note: Latex file is manually edited for formatting purposes
*/

********************************************************************************
* Table 3: See "3_impacts.do"
********************************************************************************

********************************************************************************
* Tables 4, 5: See "Tables 4, 5 - Calculations - 07-01-18.xlsx" for calculations
********************************************************************************

********************************************************************************
* Linearity test (see: Footnote 19)
********************************************************************************

use `data'/demand.dta, clear

gen takeuptemp=takeup*100
drop takeup
rename takeuptemp takeup

reg takeup tx1 tx2 tx3 ///
	hhsize wall age someschool senior chickens notfarmer bank employed ///
	busia funded market electrification population ///
	if connected==0, vce(cluster siteno)
test (_b[tx3]-_b[tx2])/15 = (_b[tx2]-_b[tx1])/10 = (_b[tx1]-_b[_cons])/10

* Note: F(2, 149)=23.03, Prob > F=0.0000 --> Reject linearity, stick to non-parametric regressions

