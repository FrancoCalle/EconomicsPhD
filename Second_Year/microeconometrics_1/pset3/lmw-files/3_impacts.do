clear all
set more off
local rootDir = "/Users/Research/LMW2019"

local data "`rootDir'/data/stata"
local temp "`rootDir'/data/temp"

* annual average exchange rates
local xr_14 = 87.94
local xr_16 = 101.53
local xr_17 = 103.41

* annual average inflation
local i_15 = 0.0801
local i_16 = 0.0635
local i_17 = 0.0450

local x_16 = (1+`i_16')*(1+`i_15')*`xr_14'
di `x_16'
* 101.015 KES/USD
local x_17 = (1+`i_17')*(1+`i_16')*(1+`i_15')*`xr_14'
di `x_17'
* 105.561

********************************************************************************
* Table 3, B6A, B6B, B6C - Impacts
* Last updated: 23 April 2019
********************************************************************************

********************************************************************************
* 1. Append data
********************************************************************************

use `data'/impacts_1.dta, clear
append using `data'/impacts_s.dta
local base base_energyspending_ksh base_asset_value_ksh
local r1 gridelec_spending_ksh kerosene_spending_ksh energy_spending_ksh asset_value_ksh percapcons_ksh
foreach var in `base' {
	gen `var'_x=`var'/`xr_14'
	qui sum `var', detail
	local label : var lab `var'
	la var `var'_x "`label' (converted to USD)"
}
foreach var in `r1' {
	gen `var'_x=`var'/`x_16'
	qui sum `var', detail
	local label : var lab `var'
	la var `var'_x "`label' (converted to USD)"
}
save `temp'/temp.dta, replace

use `data'/impacts_2.dta, clear
local base base_energyspending_ksh base_asset_value_ksh
local r2 gridelec_spending_ksh tot_hh_month_earn_ksh kerosene_spending_ksh energy_spending_ksh asset_value_ksh percapcons_ksh
foreach var in `base' {
	gen `var'_x=`var'/`xr_14'
	qui sum `var', detail
	local label : var lab `var'
	la var `var'_x "`label' (converted to USD)"
}
foreach var in `r2' {
	gen `var'_x=`var'/`x_17'
	qui sum `var', detail
	local label : var lab `var'
	la var `var'_x "`label' (converted to USD)"
}
append using `temp'/temp.dta
rm `temp'/temp.dta
gen round2=(round=="r2")
la var round2 "Round 2 (=1)"
keep reppid siteno hhid child treat_low treat_med treat_full treated_all ///
	busia base_market base_population base_connected_rate base_transearly ///
	base_connected base_housing base_bank base_energyspending_ksh_x base_female base_age base_educ base_asset_value_ksh_x ///
	female age_resp child_age child_female sibs ///
	connected connected_r gridelec_spending_ksh_x tot_hh_month_earn_ksh_x fraction_employed_all hours_worked asset_value_ksh_x percapcons_ksh_x symptoms_index life_index st_studenttest st_score_kcpe knowledge_index r_crime_index ///
	electric_lighting number_appliances mobilephone radio television iron kerosene_spending_ksh_x energy_spending_ksh_x solar_shs ///
	EMI NEMI moved round round2 metered survey

order reppid siteno hhid child treat_low treat_med treat_full treated_all ///
	busia base_market base_population base_connected_rate base_transearly ///
	base_connected base_housing base_bank base_energyspending_ksh_x base_female base_age base_educ base_asset_value_ksh_x ///
	female age_resp child_age child_female sibs ///
	connected connected_r gridelec_spending_ksh_x tot_hh_month_earn_ksh_x fraction_employed_all hours_worked asset_value_ksh_x percapcons_ksh_x symptoms_index life_index st_studenttest st_score_kcpe knowledge_index r_crime_index ///
	electric_lighting number_appliances mobilephone radio television iron kerosene_spending_ksh_x energy_spending_ksh_x solar_shs ///
	EMI NEMI moved round round2 metered survey

foreach var in base_energyspending_ksh_x base_asset_value_ksh_x gridelec_spending_ksh_x tot_hh_month_earn_ksh_x asset_value_ksh_x percapcons_ksh_x kerosene_spending_ksh_x energy_spending_ksh_x  {
   	local newname = subinstr("`var'", "_ksh_x", "",1)
   	rename `var' `newname'
}

save `temp'/impacts.dta, replace


********************************************************************************
* 2. Energy outcomes
********************************************************************************

foreach set in r1 r2 s pool {
	set more off
	use `temp'/impacts.dta, clear
	local round = "`set'" 
	if "`round'"!="pool" {
		keep if round=="`round'" & child==0
	}
	if "`round'"=="pool" {
		keep if round=="r1" & child==0 | round=="r2" & child==0
		replace round="pool"
	}
save `temp'/temp.dta, replace

local treatments "treat_low treat_med treat_full"
local connected connected
local connected_r1 connected
local connected_r2 connected
local connected_pool connected
local connected_s connected_r
local energy gridelec_spending electric_lighting number_appliances mobilephone radio television iron kerosene_spending energy_spending solar_shs
local community busia base_market base_transearly base_connected_rate base_population 
local household_r1 female age base_educ base_bank base_housing base_asset_value base_energyspending
local household_r2 female age base_educ base_bank base_housing base_asset_value base_energyspending
local household_pool female age base_educ base_bank base_housing base_asset_value base_energyspending round2
local household_s female age

foreach x in connected electric_lighting solar_shs mobilephone radio television iron {
	replace `x'= 100 if `x'==1
}

replace connected_r=connected_r/100

* a) control, itt, tot
foreach x in `connected'{
	di ""
	di ""
	di "********************************************************************************************"
	di "DepVar: `x'"
	di "Round: `round'"
	di ""
	di "(1) Mean values for the control group"
	su `x' if treat_low==0 & treat_med==0 & treat_full==0 & round=="`round'"
	local m1_`x' = round(`r(mean)', .01)
	local sd1_`x' = round(`r(sd)', .01)
	di "m: `m1_`x''"
	di "sd: `sd1_`x''"
	di ""
	di "(2) ITT regressions w/ 0% and 100% subsidy groups"	
	reg `x' treat_full `community' `household_`round'' if treat_low==0 & treat_med==0 & round=="`round'", vce(cluster siteno)
	mat temp = r(table)
	local b1_`x' = round(temp[1,1], .01)
	local se1_`x' = round(temp[2,1], .001)
	local p1_`x' = round(temp[4,1], .001)
	local n1_`x' = e(N)
	di "b: `b1_`x''"
	di "se: `se1_`x''"
	di "pval: `p1_`x''"
	di "n: `n1_`x''"
	local b2_`x' = .
	local se2_`x' = .
	local p2_`x' = .
	local n2_`x' = .
	}

foreach x in `energy'{
	di ""
	di ""
	di "********************************************************************************************"
	di "DepVar: `x'"
	di "Round: `round'"
	di ""
	di "(1) Mean values for the control group"
	su `x' if treat_low==0 & treat_med==0 & treat_full==0 & round=="`round'"
	local m1_`x' = round(`r(mean)', .01)
	local sd1_`x' = round(`r(sd)', .001)
	di "m: `m1_`x''"
	di "sd: `sd1_`x''"
	di ""
	di "(2) ITT regressions w/ 0% and 100% subsidy groups"	
	reg `x' treat_full `community' `household_`round'' if treat_low==0 & treat_med==0 & round=="`round'", vce(cluster siteno)
	mat temp = r(table)
	local b1_`x' = round(temp[1,1], .01)
	local se1_`x' = round(temp[2,1], .001)
	local p1_`x' = round(temp[4,1], .001)
	local n1_`x' = e(N)
	di "b: `b1_`x''"
	di "se: `se1_`x''"
	di "pval: `p1_`x''"
	di "n: `n1_`x''"
	di ""
	di "(3) TOT regressions w/ 1 endog. var. and 3 instruments"	
	replace connected= 1 if connected==100
	ivregress 2sls `x' `community' `household_`round'' (`connected_`round'' = `treatments' `community' `household_`round'') if round=="`round'", vce(cluster siteno)
	mat temp = r(table)
	local b2_`x' = round(temp[1,1], .01)
	local se2_`x' = round(temp[2,1], .001)
	local p2_`x' = round(temp[4,1], .001)
	local n2_`x' = e(N)
	di "b: `b2_`x''"
	di "se: `se2_`x''"
	di "pval: `p2_`x''"
	di "n: `n2_`x''"
	replace connected= 100 if connected==1
}

* b) FDR q-values

use `temp'/temp.dta, clear

local treatments "treat_low treat_med treat_full"
local connected_r1 connected
local connected_r2 connected
local connected_pool connected
local connected_s connected_r
local energy electric_lighting number_appliances mobilephone radio television iron kerosene_spending energy_spending solar_shs
local community busia base_market base_transearly base_connected_rate base_population
local household female age base_educ base_bank base_housing base_asset_value base_energyspending round2

foreach x in electric_lighting solar_shs mobilephone radio television iron {
	replace `x'= 100 if `x'==1
}

foreach x in `energy' {
	ivregress 2sls `x' `community' `household_`round'' (`connected_`round'' = `treatments') if round=="`round'", vce(cluster siteno)
	mat temp 		= r(table)
	local p_`x' 	= temp[4,1]
	local b_`x'	= temp[1,1]
	di "pval: `p_`x''"
	di "b: `b_`x''"
}

drop _all
set obs 9
g var_name = ""
g b = .
g pval = .

local n = 1
	
foreach x in `energy' {
	replace var_name = "`x'" if _n == `n'
	replace b = `b_`x'' if _n == `n'
	replace pval = `p_`x'' if _n == `n'
	local n = `n'+1
}
	
* collect the total number of p-values tested
quietly sum pval
local totalpvals = r(N)
	
* sort the p-values in ascending order and generate a variable that codes each p-value's rank
quietly gen int original_sorting_order = _n
quietly sort pval
quietly gen int rank = _n if pval~=.

* set the initial counter to 1 
local qval = 1

* generate the variable that will contain the BH (1995) q-values
gen bh95_qval = 1 if pval~=.

* set up a loop that begins by checking which hypotheses are rejected at q = 1.000, then checks which hypotheses are rejected at q = 0.999, then checks which hypotheses are rejected at q = 0.998, etc.  The loop ends by checking which hypotheses are rejected at q = 0.001.
	
while `qval' > 0 {
	* generate value qr/M
	quietly gen fdr_temp = `qval'*rank/`totalpvals'
		
	* generate binary variable checking condition p(r) <= qr/M
	quietly gen reject_temp = (fdr_temp>=pval) if fdr_temp~=.
		
	* generate variable containing p-value ranks for all p-values that meet above condition
	quietly gen reject_rank = reject_temp*rank
		
	* record the rank of the largest p-value that meets above condition
	quietly egen total_rejected = max(reject_rank)
		
	* a p-value has been rejected at level q if its rank is less than or equal to the rank of the max p-value that meets the above condition
	replace bh95_qval = `qval' if rank <= total_rejected & rank~=.
		
	* reduce q by 0.001 and repeat loop
	quietly drop fdr_temp reject_temp reject_rank total_rejected
	local qval = `qval' - .001
}
	
quietly sort original_sorting_order
pause off
set more on

di "Benjamini Hochberg (1995) q-vals are in variable 'bh95_qval'"
di	"Sorting order is the same as the original vector of p-values"	

}

********************************************************************************
* 3. Primary outcomes
********************************************************************************

foreach set in r1 r2 s pool {
	set more off
	use `temp'/impacts.dta, clear
	local round = "`set'" 
	if "`round'"!="pool" {
		keep if round=="`round'"
	}
	if "`round'"=="pool" {
		keep if round=="r1" | round=="r2"
		replace round="pool"
	}
save `temp'/temp.dta, replace

local treatments "treat_low treat_med treat_full"
local connected connected
local connected_r1 connected
local connected_r2 connected
local connected_pool connected
local connected_s connected_r
local primary_r1a gridelec_spending fraction_employed_all hours_worked asset_value percapcons symptoms_index life_index
local primary_r1b st_studenttest
local primary_r1c knowledge_index EMI NEMI NEMI NEMI NEMI
local primary_sa gridelec_spending fraction_employed_all hours_worked asset_value percapcons symptoms_index life_index
local primary_sb st_studenttest
local primary_sc knowledge_index EMI NEMI NEMI NEMI NEMI
local primary_r2a gridelec_spending fraction_employed_all tot_hh_month_earn hours_worked asset_value percapcons symptoms_index life_index
local primary_r2b st_score_kcpe
local primary_r2c r_crime_index EMI NEMI NEMI NEMI
local primary_poola gridelec_spending fraction_employed_all tot_hh_month_earn hours_worked asset_value percapcons symptoms_index life_index
local primary_poolb st_studenttest st_score_kcpe
local primary_poolc knowledge_index r_crime_index EMI NEMI
local community busia base_market base_transearly base_connected_rate base_population 
local household_r1 female age base_educ base_bank base_housing base_asset_value base_energyspending
local household_r2 female age base_educ base_bank base_housing base_asset_value base_energyspending
local household_pool female age base_educ base_bank base_housing base_asset_value base_energyspending round2
local household_s female age
local students_r1 base_educ base_bank base_housing base_asset_value base_energyspending child_female child_age sibs
local students_r2 base_educ base_bank base_housing base_asset_value base_energyspending child_female child_age sibs
local students_pool base_educ base_bank base_housing base_asset_value base_energyspending child_female child_age sibs round2
local students_s child_female child_age sibs

foreach x in connected {
	replace `x'= 100 if `x'==1
}
replace connected_r=connected_r/100

* a) control, itt, tot
foreach x in `connected'{
	di ""
	di ""
	di "********************************************************************************************"
	di "DepVar: `x'"
	di ""
	di "(1) Mean values for the control group"
	su `x' if treat_low==0 & treat_med==0 & treat_full==0 & child==0 & round=="`round'"
	local m1_`x' = round(`r(mean)', .01)
	local sd1_`x' = round(`r(sd)', .01)
	di "m: `m1_`x''"
	di "sd: `sd1_`x''"
	di ""
	di "(2) ITT regressions w/ 0% and 100% subsidy groups"	
	reg `x' treat_full `community' `household_`round'' if treat_low==0 & treat_med==0 & child==0 & round=="`round'", vce(cluster siteno)
	mat temp = r(table)
	local b1_`x' = round(temp[1,1], .01)
	local se1_`x' = round(temp[2,1], .001)
	local p1_`x' = round(temp[4,1], .001)
	local n1_`x' = e(N)
	di "b: `b1_`x''"
	di "se: `se1_`x''"
	di "pval: `p1_`x''"
	di "n: `n1_`x''"
	local b2_`x' = .
	local se2_`x' = .
	local p2_`x' = .
	local n2_`x' = .
	}

foreach x in `primary_`round'a'{
	di ""
	di ""
	di "********************************************************************************************"
	di "DepVar: `x'"
	di ""
	di "(1) Mean values for the control group"
	su `x' if treat_low==0 & treat_med==0 & treat_full==0 & child==0 & round=="`round'"
	local m1_`x' = round(`r(mean)', .01)
	local sd1_`x' = round(`r(sd)', .001)
	di "m: `m1_`x''"
	di "sd: `sd1_`x''"
	di ""
	di "(2) ITT regressions w/ 0% and 100% subsidy groups"	
	reg `x' treat_full `community' `household_`round'' if treat_low==0 & treat_med==0 & child==0 & round=="`round'", vce(cluster siteno)
	mat temp = r(table)
	local b1_`x' = round(temp[1,1], .01)
	local se1_`x' = round(temp[2,1], .001)
	local p1_`x' = round(temp[4,1], .001)
	local n1_`x' = e(N)
	di "b: `b1_`x''"
	di "se: `se1_`x''"
	di "pval: `p1_`x''"
	di "n: `n1_`x''"
	di ""
	di "(3) TOT regressions w/ 1 endog. var. and 3 instruments"	
	replace connected= 1 if connected==100
	ivregress 2sls `x' `community' `household_`round'' (`connected_`round'' = `treatments' `community' `household_`round'') if child==0 & round=="`round'", vce(cluster siteno)
	mat temp = r(table)
	local b2_`x' = round(temp[1,1], .01)
	local se2_`x' = round(temp[2,1], .001)
	local p2_`x' = round(temp[4,1], .001)
	local n2_`x' = e(N)
	di "b: `b2_`x''"
	di "se: `se2_`x''"
	di "pval: `p2_`x''"
	di "n: `n2_`x''"
	replace connected= 100 if connected==1
}

foreach x in `primary_`round'b'{
	di ""
	di ""
	di "********************************************************************************************"
	di "DepVar: `x'"
	di ""
	di "(1) Mean values for the control group"
	su `x' if treat_low==0 & treat_med==0 & treat_full==0 & child==1 & round=="`round'"
	local m1_`x' = round(`r(mean)', .01)
	local sd1_`x' = round(`r(sd)', .001)
	di "m: `m1_`x''"
	di "sd: `sd1_`x''"
	di ""
	di "(2) ITT regressions w/ 0% and 100% subsidy groups"	
	reg `x' treat_full `community' `students_`round'' if treat_low==0 & treat_med==0 & child==1 & round=="`round'", vce(cluster siteno)
	mat temp = r(table)
	local b1_`x' = round(temp[1,1], .01)
	local se1_`x' = round(temp[2,1], .001)
	local p1_`x' = round(temp[4,1], .001)
	local n1_`x' = e(N)
	di "b: `b1_`x''"
	di "se: `se1_`x''"
	di "pval: `p1_`x''"
	di "n: `n1_`x''"
	di ""
	di "(3) TOT regressions w/ 1 endog. var. and 3 instruments"	
	replace connected= 1 if connected==100
	ivregress 2sls `x' `community' `students_`round'' (`connected_`round'' = `treatments' `community' `students_`round'') if child==1 & round=="`round'", vce(cluster siteno)
	mat temp = r(table)
	local b2_`x' = round(temp[1,1], .01)
	local se2_`x' = round(temp[2,1], .001)
	local p2_`x' = round(temp[4,1], .001)
	local n2_`x' = e(N)
	di "b: `b2_`x''"
	di "se: `se2_`x''"
	di "pval: `p2_`x''"
	di "n: `n2_`x''"
	replace connected= 100 if connected==1
}

foreach x in `primary_`round'c'{
	di ""
	di ""
	di "********************************************************************************************"
	di "DepVar: `x'"
	di ""
	di "(1) Mean values for the control group"
	su `x' if treat_low==0 & treat_med==0 & treat_full==0 & child==0 & round=="`round'"
	local m1_`x' = round(`r(mean)', .01)
	local sd1_`x' = round(`r(sd)', .001)
	di "m: `m1_`x''"
	di "sd: `sd1_`x''"
	di ""
	di "(2) ITT regressions w/ 0% and 100% subsidy groups"	
	reg `x' treat_full `community' `household_`round'' if treat_low==0 & treat_med==0 & child==0 & round=="`round'", vce(cluster siteno)
	mat temp = r(table)
	local b1_`x' = round(temp[1,1], .01)
	local se1_`x' = round(temp[2,1], .001)
	local p1_`x' = round(temp[4,1], .001)
	local n1_`x' = e(N)
	di "b: `b1_`x''"
	di "se: `se1_`x''"
	di "pval: `p1_`x''"
	di "n: `n1_`x''"
	di ""
	di "(3) TOT regressions w/ 1 endog. var. and 3 instruments"	
	replace connected= 1 if connected==100
	ivregress 2sls `x' `community' `household_`round'' (`connected_`round'' = `treatments' `community' `household_`round'') if child==0 & round=="`round'", vce(cluster siteno)
	mat temp = r(table)
	local b2_`x' = round(temp[1,1], .01)
	local se2_`x' = round(temp[2,1], .001)
	local p2_`x' = round(temp[4,1], .001)
	local n2_`x' = e(N)
	di "b: `b2_`x''"
	di "se: `se2_`x''"
	di "pval: `p2_`x''"
	di "n: `n2_`x''"
	replace connected= 100 if connected==1
}

* b) FDR q-values
use `temp'/temp.dta, clear

local treatments "treat_low treat_med treat_full"
local connected connected
local connected connected
local connected_r1 connected
local connected_r2 connected
local connected_s connected_r
local primary_r1a fraction_employed_all hours_worked asset_value percapcons symptoms_index life_index
local primary_r1b st_studenttest
local primary_r1c knowledge_index
local primary_sa fraction_employed_all hours_worked asset_value percapcons symptoms_index life_index
local primary_sb st_studenttest
local primary_sc knowledge_index
local primary_r2a fraction_employed_all tot_hh_month_earn hours_worked asset_value percapcons symptoms_index life_index
local primary_r2b st_score_kcpe
local primary_r2c r_crime_index
local primary_poola fraction_employed_all tot_hh_month_earn hours_worked asset_value percapcons symptoms_index life_index
local primary_poolb st_studenttest st_score_kcpe
local primary_poolc knowledge_index r_crime_index
local community busia base_market base_transearly base_connected_rate base_population 
local household_r1 female age base_educ base_bank base_housing base_asset_value base_energyspending
local household_r2 female age base_educ base_bank base_housing base_asset_value base_energyspending
local household_pool female age base_educ base_bank base_housing base_asset_value base_energyspending round2
local household_s female age
local students_r1 base_educ base_bank base_housing base_asset_value base_energyspending child_female child_age sibs
local students_r2 base_educ base_bank base_housing base_asset_value base_energyspending child_female child_age sibs
local students_pool base_educ base_bank base_housing base_asset_value base_energyspending child_female child_age sibs round2
local students_s child_female child_age sibs

replace connected_r=connected_r/100

foreach x in `primary_`round'a' {
	ivregress 2sls `x' `community' `household_`round'' (`connected_`round'' = `treatments' `community' `household_`round'') if child==0 & round=="`round'", vce(cluster siteno)
	mat temp 		= r(table)
	local p_`x' 	= temp[4,1]
	local b_`x'	= temp[1,1]
	di "pval: `p_`x''"
	di "b: `b_`x''"
}

foreach x in `primary_`round'b' {
	ivregress 2sls `x' `community' `students_`round'' (`connected_`round'' = `treatments' `community' `students_`round'') if child==1 & round=="`round'", vce(cluster siteno)
	mat temp 		= r(table)
	local p_`x' 	= temp[4,1]
	local b_`x'	= temp[1,1]
	di "pval: `p_`x''"
	di "b: `b_`x''"
}

foreach x in `primary_`round'c' {
	ivregress 2sls `x' `community' `household_`round'' (`connected_`round'' = `treatments' `community' `household_`round'') if child==0 & round=="`round'", vce(cluster siteno)
	mat temp 		= r(table)
	local p_`x' 	= temp[4,1]
	local b_`x'	= temp[1,1]
	di "pval: `p_`x''"
	di "b: `b_`x''"
}

drop _all
set obs 11
g var_name = ""
g b = .
g pval = .

local n = 1
	
foreach x in `primary_`round'a' `primary_`round'b' `primary_`round'c'  {
	replace var_name = "`x'" if _n == `n'
	replace b = `b_`x'' if _n == `n'
	replace pval = `p_`x'' if _n == `n'
	local n = `n'+1
}
	
* collect the total number of p-values tested
quietly sum pval
local totalpvals = r(N)
	
* sort the p-values in ascending order and generate a variable that codes each p-value's rank
quietly gen int original_sorting_order = _n
quietly sort pval
quietly gen int rank = _n if pval~=.

* set the initial counter to 1 
local qval = 1

* generate the variable that will contain the BH (1995) q-values
gen bh95_qval = 1 if pval~=.

* set up a loop that begins by checking which hypotheses are rejected at q = 1.000, then checks which hypotheses are rejected at q = 0.999, then checks which hypotheses are rejected at q = 0.998, etc.  The loop ends by checking which hypotheses are rejected at q = 0.001.
	
while `qval' > 0 {
	* generate value qr/M
	quietly gen fdr_temp = `qval'*rank/`totalpvals'
		
	* generate binary variable checking condition p(r) <= qr/M
	quietly gen reject_temp = (fdr_temp>=pval) if fdr_temp~=.
		
	* generate variable containing p-value ranks for all p-values that meet above condition
	quietly gen reject_rank = reject_temp*rank
		
	* record the rank of the largest p-value that meets above condition
	quietly egen total_rejected = max(reject_rank)
		
	* a p-value has been rejected at level q if its rank is less than or equal to the rank of the max p-value that meets the above condition
	replace bh95_qval = `qval' if rank <= total_rejected & rank~=.
		
	* reduce q by 0.001 and repeat loop
	quietly drop fdr_temp reject_temp reject_rank total_rejected
	local qval = `qval' - .001
}
	
quietly sort original_sorting_order
pause off
set more on

di "Benjamini Hochberg (1995) q-vals are in variable 'bh95_qval'"
di	"Sorting order is the same as the original vector of p-values"	

}
rm `temp'/impacts.dta
rm `temp'/temp.dta









