clear
set mem 500m
set matsize 500 


global root = "C:\Users\cdippel\Dropbox\SharedFolders\NAI\empirics"  


use "$root\forcedcoexistence_webfinal.dta",clear 

************* // TABLE 1
sum logpcinc if FC == 0 & HC == 0		& year == 2000
sum logpcinc if FC == 0 & HC == 1		& year == 2000
sum logpcinc if FC == 1 & HC == 0		& year == 2000
sum logpcinc if FC == 1 & HC == 1		& year == 2000

count if FC == 0 & HC == 0		& year == 2000
count if FC == 0 & HC == 1		& year == 2000
count if FC == 1 & HC == 0		& year == 2000
count if FC == 1 & HC == 1		& year == 2000



************* // TABLE 2  (run with ivreg2 in order to do two-way clustering ;  note the "sm" adjustment at the end) 
quietly ivreg2 logpcinc_co 		FC						if year == 2000 , cluster(eaid statenumber ) sm		
	outreg2 FC using NAI_130606balance, tstat dec(3) excel replace
quietly ivreg2 logunempl_co 		FC						if year == 2000 , cluster(eaid statenumber ) sm
	outreg2 FC using NAI_130606balance, tstat dec(3) excel append
quietly ivreg2 logdist 		 	FC						if year == 2000 , cluster(eaid statenumber ) sm
	outreg2 FC using NAI_130606balance, tstat dec(3) excel append
quietly ivreg2 logruggedness 		FC						if year == 2000 , cluster(eaid statenumber )sm
	outreg2 FC using NAI_130606balance, tstat dec(3) excel append
quietly ivreg2 logresarea_sqkm  	FC						if year == 2000 , cluster(eaid statenumber ) sm
	outreg2 FC using NAI_130606balance, tstat dec(3) excel append
quietly ivreg2 HC 				FC						if year == 2000 , cluster(eaid statenumber ) sm
	outreg2 FC using NAI_130606balance, tstat dec(3) excel append
quietly ivreg2 ea_v5 			FC						if year == 2000 , cluster(eaid statenumber ) sm
	outreg2 FC using NAI_130606balance, tstat dec(3) excel append
quietly ivreg2 ea_v30 			FC						if year == 2000 , cluster(eaid statenumber ) sm
	outreg2 FC using NAI_130606balance, tstat dec(3) excel append
quietly ivreg2 ea_v32 			FC						if year == 2000 , cluster(eaid statenumber ) sm
	outreg2 FC using NAI_130606balance, tstat dec(3) excel append
quietly ivreg2 ea_v66 			FC						if year == 2000 , cluster(eaid statenumber ) sm
	outreg2 FC using NAI_130606balance, tstat dec(3) excel append
quietly ivreg2 logpop 			FC						if year == 2000 , cluster(eaid statenumber ) sm // & geo_id2 != 2430
	outreg2 FC using NAI_130606balance, tstat dec(3) excel append
quietly ivreg2 popadultshare 	FC						if year == 2000 , cluster(eaid statenumber ) sm
	outreg2 FC using NAI_130606balance, tstat dec(3) excel append
quietly ivreg2 casino 			FC						if year == 2000 , cluster(eaid statenumber ) sm
	outreg2 FC using NAI_130606balance, tstat dec(3) excel append
	
quietly ivreg2 logpcinc_co 		FC	HC					if year == 2000 , cluster(eaid statenumber ) sm		
	outreg2 FC using NAI_130606balance, tstat dec(3) excel append
quietly ivreg2 logunempl_co 		FC	HC						if year == 2000 , cluster(eaid statenumber ) sm
	outreg2 FC using NAI_130606balance, tstat dec(3) excel append
quietly ivreg2 logdist 		 	FC	HC					if year == 2000 , cluster(eaid statenumber ) sm
	outreg2 FC using NAI_130606balance, tstat dec(3) excel append
quietly ivreg2 logruggedness 	FC	HC						if year == 2000 , cluster(eaid statenumber )sm
	outreg2 FC using NAI_130606balance, tstat dec(3) excel append
quietly ivreg2 logresarea_sqkm  	FC	HC						if year == 2000 , cluster(eaid statenumber ) sm
	outreg2 FC using NAI_130606balance, tstat dec(3) excel append

quietly ivreg2 ea_v5 			FC	HC						if year == 2000 , cluster(eaid statenumber ) sm
	outreg2 FC using NAI_130606balance, tstat dec(3) excel append
quietly ivreg2 ea_v30 			FC	HC						if year == 2000 , cluster(eaid statenumber ) sm
	outreg2 FC using NAI_130606balance, tstat dec(3) excel append
quietly ivreg2 ea_v32 			FC	HC						if year == 2000 , cluster(eaid statenumber ) sm
	outreg2 FC using NAI_130606balance, tstat dec(3) excel append
quietly ivreg2 ea_v66 			FC	HC						if year == 2000 , cluster(eaid statenumber ) sm
	outreg2 FC using NAI_130606balance, tstat dec(3) excel append
quietly ivreg2 logpop 			FC	HC						if year == 2000 , cluster(eaid statenumber ) sm // 
	outreg2 FC using NAI_130606balance, tstat dec(3) excel append
quietly ivreg2 popadultshare 	FC	HC						if year == 2000 , cluster(eaid statenumber ) sm
	outreg2 FC using NAI_130606balance, tstat dec(3) excel append
quietly ivreg2 casino 			FC	HC						if year == 2000 , cluster(eaid statenumber ) sm
	outreg2 FC using NAI_130606balance, tstat dec(3) excel append

	
	
************* // TABLE 3  (not that with two-way clustering and fixed effects, STATA usually complains about the covariance matrix of moment cnoditions) 
xi 			i.statenumber	i.eaid		// generates fixed effects 

	set more off
quietly ivreg2 logpcinc FC HC 																																					if year == 2000 , cluster(eaid statenumber) sm
	outreg2 FC HC using NAI_T3, tstat dec(3) excel replace
quietly ivreg2 logpcinc FC HC 	logpcinc_co logunempl_co logdist logruggedness  logresarea_sqkm 																				if year == 2000	, cluster(eaid statenumber) sm
	outreg2 FC HC using NAI_T3, tstat dec(3) excel append
quietly ivreg2 logpcinc FC HC 	logpcinc_co logunempl_co logdist logruggedness  logresarea_sqkm ea_v5 ea_v30 ea_v32 ea_v66														if year == 2000	, cluster(eaid statenumber) sm
	outreg2 FC HC using NAI_T3, tstat dec(3) excel append
quietly ivreg2 logpcinc FC HC 	logpcinc_co logunempl_co logdist logruggedness  logresarea_sqkm ea_v5 ea_v30 ea_v32 ea_v66	logpop logpopsq popadultshare casino 				if year == 2000	, cluster(eaid statenumber) sm
	outreg2 FC HC using NAI_T3, tstat dec(3) excel append
quietly ivreg2 logpcinc FC HC 	logpcinc_co logunempl_co logdist logruggedness  logresarea_sqkm ea_v5 ea_v30 ea_v32 ea_v66	logpop logpopsq popadultshare casino _Ist*			if year == 2000	, cluster(eaid statenumber) sm
	outreg2 FC HC using NAI_T3, tstat dec(3) excel append
quietly ivreg2 logpcinc FC HC 																																			_Ie*	if year == 2000 , cluster(eaid statenumber) sm
	outreg2 FC  using NAI_T3, tstat dec(3) excel append
quietly ivreg2 logpcinc FC HC 	logpcinc_co logunempl_co logdist logruggedness  logresarea_sqkm 																		_Ie*	if year == 2000	, cluster(eaid statenumber) sm
	outreg2 FC  using NAI_T3, tstat dec(3) excel append
quietly ivreg2 logpcinc FC HC 	logpcinc_co logunempl_co logdist logruggedness  logresarea_sqkm ea_v5 ea_v30 ea_v32 ea_v66												_Ie*	if year == 2000	, cluster(eaid statenumber) sm
	outreg2 FC  using NAI_T3, tstat dec(3) excel append
quietly ivreg2 logpcinc FC HC 	logpcinc_co logunempl_co logdist logruggedness  logresarea_sqkm ea_v5 ea_v30 ea_v32 ea_v66	logpop logpopsq popadultshare casino 		_Ie*	if year == 2000	, cluster(eaid statenumber) sm
	outreg2 FC  using NAI_T3, tstat dec(3) excel append
quietly ivreg2 logpcinc FC HC 	logpcinc_co logunempl_co logdist logruggedness  logresarea_sqkm ea_v5 ea_v30 ea_v32 ea_v66	logpop logpopsq popadultshare casino _Ist*	_Ie*	if year == 2000	, cluster(eaid statenumber) sm
	outreg2 FC  using NAI_T3, tstat dec(3) excel append
	

	
************* // TABLES 4	+ 5

use "$root\forcedcoexistence_webfinal.dta",clear 

xi i.statenumber
gen instrument_precious = instrument_gold + intrument_silver

gen wprec_enviro = wgold_enviro + wsilver_enviro 


// FS
set more off
quietly ivreg2 FC instrument_gold  intrument_silver HC 																																													if year == 2000 , cl(eaid statenumber) sm
	outreg2 instrument_gold  intrument_silver using NAI_T4, tstat dec(3) excel replace
quietly ivreg2 FC instrument_gold  intrument_silver HC logpcinc_co logunempl_co logdist logruggedness logresarea_sqkm 																													if year == 2000 , cl(eaid statenumber)  sm
	outreg2 instrument_gold  intrument_silver using NAI_T4, tstat dec(3) excel append
quietly ivreg2 FC instrument_gold  intrument_silver HC logpcinc_co logunempl_co logdist logruggedness logresarea_sqkm ea_v5  ea_v30 ea_v32 ea_v66  																						if year == 2000 , cl(eaid statenumber)  sm
	outreg2 instrument_gold  intrument_silver using NAI_T4, tstat dec(3) excel append
quietly ivreg2 FC instrument_gold  intrument_silver HC logpcinc_co logunempl_co logdist logruggedness logresarea_sqkm ea_v5  ea_v30 ea_v32 ea_v66  	logpop  popadultshare casino														if year == 2000 , cl(eaid statenumber)  sm
	outreg2 instrument_gold  intrument_silver using NAI_T4, tstat dec(3) excel append
quietly ivreg2 FC instrument_gold  intrument_silver HC logpcinc_co logunempl_co logdist logruggedness logresarea_sqkm ea_v5  ea_v30 ea_v32 ea_v66  	logpop  popadultshare casino  _Ist*													if year == 2000 , cl(eaid statenumber)  sm // partial(i.statenumber) 
	outreg2 instrument_gold  intrument_silver using NAI_T4, tstat dec(3) excel append
quietly ivreg2 FC instrument_gold  intrument_silver HC logpcinc_co logunempl_co logdist logruggedness logresarea_sqkm ea_v5  ea_v30 ea_v32 ea_v66  	logpop  popadultshare casino  _Ist*	removal  homelandruggedness wprec_enviro				if year == 2000 , cl(eaid statenumber)  sm // partial(i.statenumber) 
	outreg2 instrument_gold  intrument_silver using NAI_T4, tstat dec(3) excel append


// RF
quietly ivreg2 logpcinc instrument_gold  intrument_silver HC 																																												if year == 2000 , cl(eaid statenumber)  sm
	outreg2 instrument_gold  intrument_silver using NAI_T4, tstat dec(3) excel append
quietly ivreg2 logpcinc instrument_gold  intrument_silver HC logpcinc_co logunempl_co logdist logruggedness logresarea_sqkm 																												if year == 2000 , cl(eaid statenumber)  sm
	outreg2 instrument_gold  intrument_silver using NAI_T4, tstat dec(3) excel append
quietly ivreg2 logpcinc instrument_gold  intrument_silver HC logpcinc_co logunempl_co logdist logruggedness logresarea_sqkm ea_v5  ea_v30 ea_v32 ea_v66  																					if year == 2000 , cl(eaid statenumber)  sm
	outreg2 instrument_gold  intrument_silver using NAI_T4, tstat dec(3) excel append
quietly ivreg2 logpcinc instrument_gold  intrument_silver HC logpcinc_co logunempl_co logdist logruggedness logresarea_sqkm ea_v5  ea_v30 ea_v32 ea_v66 logpop  popadultshare casino   														if year == 2000 , cl(eaid statenumber)  sm	 
	outreg2 instrument_gold  intrument_silver using NAI_T4, tstat dec(3) excel append
quietly ivreg2 logpcinc instrument_gold  intrument_silver HC logpcinc_co logunempl_co logdist logruggedness logresarea_sqkm ea_v5  ea_v30 ea_v32 ea_v66 logpop  popadultshare casino _Ist*													if year == 2000 , cl(eaid statenumber)  sm	// partial(i.statenumber) 
	outreg2 instrument_gold  intrument_silver using NAI_T4, tstat dec(3) excel append
quietly ivreg2 logpcinc instrument_gold  intrument_silver HC logpcinc_co logunempl_co logdist logruggedness logresarea_sqkm ea_v5  ea_v30 ea_v32 ea_v66 logpop  popadultshare casino _Ist*	removal homelandruggedness wprec_enviro				if year == 2000 , cl(eaid statenumber)  sm // partial(i.statenumber) 
	outreg2 instrument_gold  intrument_silver using NAI_T4, tstat dec(3) excel append


// 2S
set more off
quietly ivreg2 logpcinc (FC = instrument_gold  intrument_silver ) HC 																																												if year == 2000 , cl(eaid statenumber)  sm endog(FC ) partial(HC )
	outreg2 using NAI_130606FS, tstat dec(3) excel replace e(idstat idp widstat widstatp j jp estat estatp)
quietly ivreg2 logpcinc (FC = instrument_gold  intrument_silver ) HC logpcinc_co logunempl_co logdist logruggedness logresarea_sqkm  																												if year == 2000 , cl(eaid statenumber)  sm endog(FC )	partial(HC logpcinc_co logunempl_co logruggedness logresarea_sqkm)
	outreg2 using NAI_130606FS, tstat dec(3) excel append e(idstat idp widstat widstatp j jp estat estatp)
quietly ivreg2 logpcinc (FC = instrument_gold  intrument_silver ) HC logpcinc_co logunempl_co logdist logruggedness logresarea_sqkm ea_v5  ea_v30 ea_v32 ea_v66  	 																				if year == 2000 , cl(eaid statenumber)  sm endog(FC ) partial(HC logpcinc_co logunempl_co logruggedness logresarea_sqkm ea_v5  ea_v30 ea_v32 ea_v66)
	outreg2 using NAI_130606FS, tstat dec(3) excel append e(idstat idp widstat widstatp j jp estat estatp)
quietly ivreg2 logpcinc (FC = instrument_gold  intrument_silver ) HC logpcinc_co logunempl_co logdist logruggedness logresarea_sqkm ea_v5  ea_v30 ea_v32 ea_v66  logpop  popadultshare casino     													if year == 2000 ,  cl(eaid statenumber)  sm endog(FC ) partial(HC logpcinc_co logunempl_co logruggedness logresarea_sqkm ea_v5  ea_v30 ea_v32 ea_v66 logpop  popadultshare casino)
	outreg2 using NAI_130606FS, tstat dec(3) excel append e(idstat idp widstat widstatp j jp estat estatp)
quietly ivreg2 logpcinc (FC = instrument_gold  intrument_silver ) HC logpcinc_co logunempl_co logdist logruggedness logresarea_sqkm ea_v5  ea_v30 ea_v32 ea_v66  logpop  popadultshare casino _Ist*													if year == 2000 , cl(eaid statenumber)  sm endog(FC ) partial(HC logpcinc_co logunempl_co logruggedness logresarea_sqkm ea_v5  ea_v30 ea_v32 ea_v66 logpop  popadultshare casino _Ist*)
	outreg2 using NAI_130606FS, tstat dec(3) excel append e(idstat idp widstat widstatp j jp estat estatp)
quietly ivreg2 logpcinc (FC = instrument_gold  intrument_silver ) HC logpcinc_co logunempl_co logdist logruggedness logresarea_sqkm ea_v5  ea_v30 ea_v32 ea_v66  logpop  popadultshare casino _Ist*	removal homelandruggedness wprec_enviro					 if year == 2000 , cl(eaid statenumber)  sm endog(FC )   partial(HC logpcinc_co logunempl_co  logruggedness logresarea_sqkm ea_v5  ea_v30 ea_v32 ea_v66 logpop  popadultshare casino  logdist _Ist* removal homelandruggedness wprec_enviro			  )
	outreg2 using NAI_130606FS, tstat dec(3) excel append e(idstat idp widstat widstatp j jp estat estatp)

			quietly ivreg2 logpcinc (FC = instrument_precious ) HC 																																											   		if year == 2000 , cl(eaid statenumber)  sm endog(FC ) partial(HC )
				outreg2 using NAI_130606FS, tstat dec(3) excel append e(idstat idp widstat widstatp j jp estat estatp)
			quietly ivreg2 logpcinc (FC = instrument_precious ) HC logpcinc_co logunempl_co logdist logruggedness logresarea_sqkm  																										     		if year == 2000 , cl(eaid statenumber)  sm endog(FC )	partial(HC logpcinc_co logunempl_co logruggedness logresarea_sqkm)
				outreg2 using NAI_130606FS, tstat dec(3) excel append e(idstat idp widstat widstatp j jp estat estatp)
			quietly ivreg2 logpcinc (FC = instrument_precious ) HC logpcinc_co logunempl_co logdist logruggedness logresarea_sqkm ea_v5  ea_v30 ea_v32 ea_v66  	 																					if year == 2000 , cl(eaid statenumber)  sm endog(FC ) partial(HC logpcinc_co logunempl_co logruggedness logresarea_sqkm ea_v5  ea_v30 ea_v32 ea_v66)
				outreg2 using NAI_130606FS, tstat dec(3) excel append e(idstat idp widstat widstatp j jp estat estatp)
			quietly ivreg2 logpcinc (FC = instrument_precious ) HC logpcinc_co logunempl_co logdist logruggedness logresarea_sqkm ea_v5  ea_v30 ea_v32 ea_v66  	logpop  popadultshare casino  								   						if year == 2000 ,  cl(eaid statenumber)  sm endog(FC ) partial(HC logpcinc_co logunempl_co logruggedness logresarea_sqkm ea_v5  ea_v30 ea_v32 ea_v66 logpop  popadultshare casino)
				outreg2 using NAI_130606FS, tstat dec(3) excel append e(idstat idp widstat widstatp j jp estat estatp)
			quietly ivreg2 logpcinc (FC = instrument_precious) HC logpcinc_co logunempl_co logdist logruggedness logresarea_sqkm ea_v5  ea_v30 ea_v32 ea_v66  	logpop  popadultshare casino _Ist*													if year == 2000 , cl(eaid statenumber)  sm endog(FC ) partial(HC logpcinc_co logunempl_co logruggedness logresarea_sqkm ea_v5  ea_v30 ea_v32 ea_v66 logpop  popadultshare casino _Ist*)
				outreg2 using NAI_130606FS, tstat dec(3) excel append e(idstat idp widstat widstatp j jp estat estatp)
			quietly ivreg2 logpcinc (FC = instrument_precious) HC logpcinc_co logunempl_co logdist logruggedness logresarea_sqkm ea_v5  ea_v30 ea_v32 ea_v66  	logpop  popadultshare casino  _Ist*	removal homelandruggedness wprec_enviro			if year == 2000 , cl(eaid statenumber)  sm endog(FC )  partial(HC logpcinc_co logunempl_co  logruggedness logresarea_sqkm ea_v5  ea_v30 ea_v32 ea_v66 logpop  popadultshare casino  logdist _Ist* removal homelandruggedness wprec_enviro			  )
				outreg2 using NAI_130606FS, tstat dec(3) excel append e(idstat idp widstat widstatp j jp estat estatp)
				
************* // TABLE 6

use "$root\forcedcoexistence_webfinal.dta",clear 

gen d_1970 = (year == 1970)
gen d_1980 = (year == 1980)
gen d_1990 = (year == 1990)
gen d_2000 = (year == 2000)

// time-varying controls for columns 3-6
local tlist	d_1970 d_1980 d_1990 d_2000  		
local xlist  FC HC casino logpcinc_co logdist logunempl_co ea_v5 ea_v30 ea_v32 ea_v66 	
foreach T in `tlist' {
		foreach X in `xlist' {
	gen `T'_`X' = 	`T'*`X'
	}
}
// casino-environment interactions for column 6
local x2list  logpcinc_co logdist logunempl_co 
foreach X in `x2list' {
	gen d_2000_casino_`X' = 	d_2000_casino*`X'
}
*

	local environmentlist logpcinc_co d_1980_logpcinc_co d_1990_logpcinc_co d_2000_logpcinc_co logdist d_1980_logdist d_1990_logdist d_2000_logdist  logunempl_co d_1980_logunempl_co d_1990_logunempl_co d_2000_logunempl_co
	local tribelist HC d_1980_HC d_1990_HC d_2000_HC   ea_v5 d_1980_ea_v5 d_1990_ea_v5 d_2000_ea_v5 ea_v30 d_1980_ea_v30 d_1990_ea_v30 d_2000_ea_v30 ea_v32 d_1980_ea_v32 d_1990_ea_v32 d_2000_ea_v32 ea_v66 d_1980_ea_v66 d_1990_ea_v66 d_2000_ea_v66
	local interactionlist  d_2000_casino_logunempl_co d_2000_casino_logdist d_2000_casino_logpcinc_co
quietly areg logpcinc d_1980 d_1990 d_2000 																													, absorb(geo_id2)	cluster(geo_id2)	// 1 
	outreg2 d_1980 d_1990 d_2000 																															using NAI_T6, tstat dec(3) excel replace
quietly areg logpcinc d_1980 d_1990 d_2000 FC d_1980_FC d_1990_FC d_2000_FC 																				, absorb(geo_id2)	cluster(geo_id2)	// 2
	outreg2 d_1980 d_1990 d_2000 d_1980_FC d_1990_FC d_2000_FC 																								using NAI_T6, tstat dec(3) excel append
quietly areg logpcinc d_1980 d_1990 d_2000 FC d_1980_FC d_1990_FC d_2000_FC 	 							`environmentlist'	`tribelist'				   	, absorb(geo_id2)	cluster(geo_id2)	// 3
	outreg2 d_1980 d_1990 d_2000 d_1980_FC d_1990_FC d_2000_FC 																								using NAI_T6, tstat dec(3) excel append
		/*
		test d_1980_logdist d_1980_logunempl_co	 
		test d_1990_logdist d_1990_logunempl_co	
		test d_2000_logdist   d_2000_logunempl_co 	
		test d_1980_orgd d_1980_v5 d_1980_v30 d_1980_v32 d_1980_v66 
		test d_1990_orgd d_1990_v5 d_1990_v30 d_1990_v32 d_1990_v66
		test d_2000_orgd d_2000_v5 d_2000_v30 d_2000_v32 d_2000_v66
		*/
quietly areg logpcinc d_1980 d_1990 d_2000 FC d_1980_FC d_1990_FC d_2000_FC	 	 														 	 d_2000_casino	, absorb(geo_id2)	cluster(geo_id2)	// 4
	outreg2 d_1980 d_1990 d_2000 d_1980_FC d_1990_FC d_2000_FC d_2000_casino 																				using NAI_T6, tstat dec(3) excel append
quietly areg logpcinc d_1980 d_1990 d_2000 FC d_1980_FC d_1990_FC d_2000_FC	 	 	 	 					`environmentlist'	`tribelist'	 d_2000_casino	, absorb(geo_id2)	cluster(geo_id2)	// 5
	outreg2 d_1980 d_1990 d_2000 d_1980_FC d_1990_FC d_2000_FC d_2000_casino 																				using NAI_T6, tstat dec(3) excel append

quietly areg logpcinc d_1980 d_1990 d_2000 FC d_1980_FC d_1990_FC d_2000_FC	 		 `interactionlist'		`environmentlist'	`tribelist'		d_2000_casino	, absorb(geo_id2)	cluster(geo_id2)	// 5
	outreg2 d_1980 d_1990 d_2000 d_1980_FC d_1990_FC d_2000_FC d_2000_casino d_2000_casino_logdist d_2000_casino_logunempl_co d_2000_casino_logpcinc_co		using NAI_T6, tstat dec(3) excel append

drop d_*
	

************* // TABLE 7

use "$root\forcedcoexistence_webfinal.dta",clear 

	
gen d_2000 = (year == 2000)
gen d_1990 = (year == 1990)
local tlist	d_1990 d_2000  
local xlist  FC 	HC ea_v5 ea_v30 ea_v32 ea_v66 	casino //	count_1_all	logpcinc_co logdist  logunempl_co logresarea_sqkm logpop logpopsq 
foreach T in `tlist' {
		foreach X in `xlist' {
	gen `T'_`X' = 	`T'*`X'
	}
}
local extension governance notgovernance	
local type  5 7 // 5 
foreach z in `type'{			// 
foreach x in `extension' {		// 
count if share_`z'_`x' != . 
count if count_`z'_`x' != . 
}
}
*

set more off

poisson count_5_governance	 FC			d_2000  							HC ea_v5 ea_v30 ea_v32 ea_v66			logcount_1_all 		d_2000_casino	logpop logpopsq 			logpcinc_co  logunempl_co	if year > 1980	, vce(cl geo_id2)
poisson count_5_governance	 FC			d_2000  							, vce(cl geo_id2)


local i=0
local tribelist d_1990_HC d_2000_HC   	d_1990_ea_v5 d_2000_ea_v5 		d_1990_ea_v30 d_2000_ea_v30 		d_1990_ea_v32 d_2000_ea_v32 			d_1990_ea_v66 d_2000_ea_v66
local xlist 	logpop logpopsq 			logpcinc_co  logunempl_co	
set more off
	local extension governance notgovernance	
	local type  5 7 
foreach z in `type'{			
foreach x in `extension' {		 
quietly poisson count_`z'_`x'	 FC			d_2000  							HC ea_v5 ea_v30 ea_v32 ea_v66			logcount_1_all 		d_2000_casino	`xlist'			, vce(cl geo_id2)
		if `i'==0{
				outreg2 FC using NAI_T7PanelA, tstat dec(3) excel  replace
		}
		else {
				outreg2 FC using NAI_T7PanelA, tstat dec(3) excel  append
		}
quietly poisson count_`z'_`x'	 		    d_2000 d_1990_FC d_2000_FC			`tribelist'	 							logcount_1_all  	d_2000_casino	`xlist' 		, vce(cl geo_id2)
	outreg2 d_1990_FC d_2000_FC 			using NAI_T7PanelA, tstat dec(3) excel append addstat(Pseudo R2, e(r2_p))

	local i=1
}
}	
*

local i=0
local tribelist d_1990_HC d_2000_HC   	d_1990_ea_v5 d_2000_ea_v5 		d_1990_ea_v30 d_2000_ea_v30 		d_1990_ea_v32 d_2000_ea_v32 			d_1990_ea_v66 d_2000_ea_v66
local xlist 	logpop logpopsq 			logpcinc_co  logunempl_co	
set more off
	local extension governance notgovernance	
	local type  5 7 
foreach z in `type'{			
foreach x in `extension' {	
quietly reg share_`z'_`x'	 FC			d_2000  							HC ea_v5 ea_v30 ea_v32 ea_v66			logcount_1_all 		d_2000_casino	`xlist'			, vce(cl geo_id2)
		if `i'==0{
				outreg2 FC using NAI_T7PanelB, tstat dec(3) excel  replace
		}
		else {
				outreg2 FC using NAI_T7PanelB, tstat dec(3) excel  append
		}
quietly reg share_`z'_`x'	 			d_2000 d_1990_FC d_2000_FC	  		`tribelist'	 							logcount_1_all 		d_2000_casino	`xlist'			, vce(cl geo_id2)
	outreg2 d_1990_FC d_2000_FC 			using NAI_T7PanelB, tstat dec(3) excel append 

	local i=1
}
}	
drop d_*



************* // TABLE 8
	
use "$root\forcedcoexistence_webfinal.dta",clear 

xi  i.statenumber

local controls	logpcinc_co logunempl_co logdist logruggedness logresarea_sqkm  ea_v5  ea_v30 ea_v32 ea_v66  logpop  popadultshare casino
quietly ivreg2  logpcinc_wage	FC HC 	`controls'	_Ist* , cl(eaid statenumber)  sm
	outreg2 FC using NAI_130606T8, tstat dec(3) excel replace
quietly ivreg2  logavg_salary	FC HC 	`controls'	_Ist* , cl(eaid statenumber)  sm
	outreg2 FC using NAI_130606T8, tstat dec(3) excel append
quietly ivreg2  share_empl		FC HC 	`controls'	_Ist* , cl(eaid statenumber)  sm
	outreg2 FC using NAI_130606T8, tstat dec(3) excel append
quietly ivreg2  logpcinc_ssi		FC HC 	`controls'	_Ist* , cl(eaid statenumber)  sm
	outreg2 FC using NAI_130606T8, tstat dec(3) excel append
quietly ivreg2   logbiatotalpc		FC HC 	`controls'	 , cl(eaid statenumber)  sm
	outreg2 FC using NAI_130606T8, tstat dec(3) excel append
quietly ivreg2  share_workinco	FC HC 	`controls' _Ist* ,cl(eaid statenumber)  sm
	outreg2 FC using NAI_130606T8, tstat dec(3) excel append
quietly ivreg2  traveltimetowork	FC HC 	`controls'	_Ist* , cl(eaid statenumber)  sm
	outreg2 FC using NAI_130606T8, tstat dec(3) excel append
quietly ivreg2  share_english	FC HC 	`controls'	_Ist* , cl(eaid statenumber)  sm
	outreg2 FC using NAI_130606T8, tstat dec(3) excel append
quietly ivreg2  share_ba 		FC HC 	`controls'	_Ist* , cl(eaid statenumber)  sm
	outreg2 FC using NAI_130606T8, tstat dec(3) excel append
	
	


			
		
		
		
		