* P.HULL
* UP Grandfathering Analysis

clear all
cap log close
set more off
set trace off
set matsize 10000
set seed 42
adopath + "Z:\ado\plus\"

global lotterydata "Z:/charters/MA Lotteries/all lotto data/"
global combdata    "Z:/charters/Combined data/"
global cleandata   "X:/indistrict/Data/"
global bpsdata     "Z:/BPS/Peter/BPSVA/Data/"
global aypdata     "Z:/charters/MA DESE data/Accountability/"
global simsdata    "Z:/charters/MA DESE data/SIMS/"

global programs    "Z:/charters/indistrict/Programs/"

local simsdate "10 30 2014"

local covstems "female black hisp asian otherrace sped lep foodfred fem_min"
local controls "foodfred_bl lep_bl c_state_mrawsc_bl c_state_erawsc_bl y_*"

matrix cohorts       = J(200,10,.)
matrix balance       = J(200,10,.)
matrix gfing         = J(200,5,.)
matrix gfingnoFD     = J(200,5,.)
matrix choice        = J(20,20,.)
matrix switching     = J(200,10,.)
matrix persistence   = J(200,10,.)

local setup          = 1
local switchingsetup = 1

local balance        = 1 /* Tab 6                */
local gfing          = 1 /* Tab 7                */
local switching      = 1 /* Tab 10 (Column 2)    */
local persistence    = 1 /* Tab A3 (Columns 4-6) */
local cohorts        = 1 /* Tab B2 (Panel B)     */
local altgfing       = 1 /* Tab B7               */
local choice         = 1 /* Tab B10  (Panel A)   */

***********************************************
* Construct analysis files                    *             
***********************************************
if `setup'==1 {
use "${cleandata}simsmcaswide `simsdate'", clear

forval g = 3/8 {
		foreach t in m e {
			qui replace c_state_`t'rawsc`g' = c_boston_`t'rawsc`g'
		}
}

gen baselinegrade = 5

foreach v in boston year `covstems' c_state_mrawsc c_state_erawsc{
	gen `v'_bl = `v'5
}

gen white_bl = 1-(hisp_bl + black_bl + asian_bl + otherrace_bl)

gen hasbaselinedemos = !missing(female_bl, black_bl, hisp_bl, asian_bl, otherrace_bl, sped_bl, lep_bl, foodfred_bl)

compress
save "${cleandata}obsfile_2014", replace		
}
if `switchingsetup'==1 {
	use "${cleandata}simsmcaswide `simsdate'", clear
	forval g = 3/8 {
			foreach t in m e {
				qui replace c_state_`t'rawsc`g' = c_boston_`t'rawsc`g'
			}
	}
	keep sasid year7 masscode7 year8 masscode8 c_state_mrawsc5 c_state_erawsc5 sped5 lep5 foodfred5
	forval outgr = 7/8 {
		preserve
		collapse (mean) c_state_mrawsc5 c_state_erawsc5 sped5 lep5 foodfred5, by(year`outgr' masscode`outgr')
		gen peer_mscore_for`outgr' = c_state_mrawsc5
		gen peer_escore_for`outgr' = c_state_erawsc5
		gen peer_sped_for`outgr' = sped5
		gen peer_lep_for`outgr' = lep5
		gen peer_foodfred_for`outgr' = foodfred5
		
		drop if year`outgr' == . | masscode`outgr' == .
		keep masscode`outgr' year`outgr' peer_*
		tempfile peer_for`outgr'
		save `peer_for`outgr'', replace
		restore
	}
	
	use "${cleandata}obsfile_2014", clear
	forval outgr = 7/8 {
		merge m:1 year`outgr' masscode`outgr' using `peer_for`outgr''
		drop _merge
	}	
	
	save "${cleandata}obsfile_2014_wpeers", replace
}
***********************************************
* Counting match & follow-up                  *
***********************************************
if `cohorts' == 1 {
	use "${simsdata}simsoct11", clear
	keep if (grade=="07" | grade=="08") & school=="04800405" 
	keep sasid
	tempfile followup
	save `followup'
	use "${simsdata}simsoct10", clear
	keep if (grade=="06" | grade=="07") & school=="00350435" 
	merge 1:1 sasid using `followup', keep(master match)
	gen temp = _merge==3
	summ temp
	matrix cohorts[1,6] = r(mean)
	count
	local gf_elig = r(N)
	keep sasid grade
	destring grade, replace
	destring sasid, replace
	gen ingavin = 1
	reshape wide ingavin, i(sasid) j(grade)
	gen year6 = 2011
	gen year7 = 2011
	tempfile gavin_fa10
	save `gavin_fa10', replace

	use "${cleandata}obsfile_2014", clear 
	egen cell = group(year5 masscode5 hisp5 black5 asian5 otherrace5 female5 sped5 foodfred5)
	
	merge 1:1 sasid using `gavin_fa10', keep(master match)
	drop _merge

	forval g = 6/7 {
		gen gf_elig`g' = 1 if ingavin`g'==1 & year`g'==2011
		replace gf_elig`g' = 0 if ingavin`g'==. & year`g'==2011
		bys cell: egen temp = max(gf_elig`g') if cell!=.
		gen samp`g' = (temp==1 & gf_elig`g'!=.) & boston5==1 & bcharter5==0 & spedalt5==0
		drop temp
	}
	replace samp7 = 0 if samp6 == 1
	gen matchedsample = (samp6==1 | samp7==1)
	egen gf_eligR = rowmax(gf_elig6 gf_elig7)
	egen gf_eligwB = rowmax(gf_elig6 gf_elig7) if cell!=.
	egen gf_eligiD = rowmax(gf_elig6 gf_elig7) if cell!=. & boston5==1 & bcharter5==0 & spedalt5==0
	gen gf_eligiS = gf_elig6 if samp6
	replace gf_eligiS = gf_elig7 if samp7

	local row = 1 
	local col = 1
	count if gf_eligR == 1
		matrix cohorts[`row',`col'] = r(N)
	local ++col
	count if gf_eligwB == 1
		matrix cohorts[`row',`col'] = r(N)
	local ++col	
	count if gf_eligiD == 1
		matrix cohorts[`row',`col'] = r(N)
	local ++col	
	count if gf_eligiS == 1
		matrix cohorts[`row',`col'] = r(N)
	local ++col
	count if matchedsample == 1 & gf_eligiS == 0
		matrix cohorts[`row',`col'] = r(N)
	
	clear
	svmat cohorts
	br
	
}

***********************************************
* Grandfathering Balance                      *
***********************************************
if `balance' == 1 {

use "${simsdata}simsoct10", clear
keep if (grade=="06" | grade=="07") & school=="00350435" 
keep sasid grade
destring grade, replace
destring sasid, replace
gen ingavin = 1
reshape wide ingavin, i(sasid) j(grade)
gen year6 = 2011
gen year7 = 2011
tempfile gavin_fa10
save `gavin_fa10', replace

use if boston5==1 & spedalt5==0 & year5>=2006 & hasbaselinedemos using "${cleandata}obsfile_2014", clear 
egen cell = group(year5 masscode5 hisp5 black5 asian5 otherrace5 female5 sped5 foodfred5)

merge 1:1 sasid using `gavin_fa10', keep(master match)
drop _merge

forval g = 6/7 {
	gen gf_elig`g' = 1 if ingavin`g'==1 & year`g'==2011
	replace gf_elig`g' = 0 if ingavin`g'==. & year`g'==2011
	bys cell: egen temp = max(gf_elig`g') if cell!=.
	gen samp`g' = (temp==1 & gf_elig`g'!=.) & bcharter5==0 & spedalt5==0
	drop temp
}
replace samp7 = 0 if samp6 == 1
gen matchedsample = (samp6==1 | samp7==1)
gen gf_elig = gf_elig6 if samp6
replace gf_elig = gf_elig7 if samp7
replace cell = . if matchedsample == 0

gen compsample = 0
qui levelsof year5 if matchedsample==1, local(years)
qui count if matchedsample==1
local denom = r(N)
local prps ""
foreach y in `years' {
	qui count if year5==`y' & matchedsample==1
	local prop`y' = r(N)/`denom'
	local prps "`prps', `prop`y''"
}
local prps = substr("`prps'",2,.)
disp "`prps'"
local maxprop = max(`prps')
disp `maxprop'
foreach y in `years' {
	replace compsample = runiform()<(`prop`y''/`maxprop') if year5==`y' & compsamp == 0 & spedalt5==0
}
gen compsamp_charter = (charter6==1 | charter7==1 | charter8==1) if charter6!=. | charter7!=. | charter8!=.

keep if matchedsample == 1 | compsample == 1

foreach g in 6 7 8 {
	gen attend_UP`g' = masscode`g'  == 4800405
	replace attend_UP`g' = . if masscode`g' == .
}
gen years1 = attend_UP7 if samp6==1
replace years1 = attend_UP8 if samp7==1

gen any_attendUP = (attend_UP7==1 | attend_UP8==1) if (attend_UP7!=. | attend_UP8!=.)

gen no_attrit0 = .
replace no_attrit0 = (c_state_mrawsc6!=. | c_state_erawsc6!=.) if samp6==1 | (compsample==1 & samp7==0)
replace no_attrit0 = (c_state_mrawsc7!=. | c_state_erawsc7!=.) if samp7==1
gen no_attrit1 = .
replace no_attrit1 = (c_state_mrawsc7!=. | c_state_erawsc7!=.) if samp6==1 | (compsample==1 & samp7==0)
replace no_attrit1 = (c_state_mrawsc8!=. | c_state_erawsc8!=.) if samp7==1
gen no_attrit2 = .
replace no_attrit2 = (c_state_mrawsc8!=. | c_state_erawsc8!=.) if samp6==1 | (compsample==1 & samp7==0)

local row = 1
local col = 1

gen mgr_bl = c_state_mrawsc5 - c_state_mrawsc4
gen egr_bl = c_state_erawsc5 - c_state_erawsc4
qui tab cell, gen(d_)

foreach dep in hisp black white asian female sped foodfred lep c_state_mrawsc mgr c_state_erawsc egr {

	summ `dep'_bl if compsample==1
		matrix balance[`row',`col']   = r(mean)
		matrix balance[`row'+2,`col'] = r(N)
	summ `dep'_bl if compsample==1 & compsamp_charter==1
		matrix balance[`row',`col'+1]   = r(mean)
		matrix balance[`row'+2,`col'+1] = r(N)
	summ `dep'_bl if matchedsample==1 & any_attendUP==1
		matrix balance[`row',`col'+2]   = r(mean)
		matrix balance[`row'+2,`col'+2] = r(N)
	summ `dep'_bl if matchedsample==1 & gf_elig==1
		matrix balance[`row',`col'+3]   = r(mean)
		matrix balance[`row'+2,`col'+3] = r(N)

	qui areg gf_elig c_state_mrawsc4 c_state_erawsc4 if matchedsample==1 & `dep'_bl!=. & c_state_mrawsc4!=. & c_state_erawsc4!=. & years1!=., absorb(cell)
	predict pZ, xbd
	replace pZ = 0.000001 if pZ<=0
	replace pZ = 0.999999 if pZ>=1
	gen kappa = 1 - gf_elig*(1-years1)/pZ - (1-gf_elig)*years1/(1-pZ)
	summ `dep'_bl if matchedsample==1 & `dep'_bl!=. & c_state_mrawsc4!=. & c_state_erawsc4!=. & years1!=. [iw=kappa]
		matrix balance[`row',`col'+4]   = r(mean)
	summ kappa if matchedsample==1 & `dep'_bl!=. & c_state_mrawsc4!=. & c_state_erawsc4!=. & years1!=.
		matrix balance[`row'+2,`col'+4] = `=`=r(N)'*`=r(mean)''
	drop pZ kappa
	summ `dep'_bl if matchedsample==1 & gf_elig==0
		matrix balance[`row',`col'+5]   = r(mean)
		matrix balance[`row'+2,`col'+5] = r(N)
	
	if "`dep'" == "lep" | "`dep'" == "c_state_mrawsc" | "`dep'" == "c_state_erawsc" {
		areg `dep'_bl gf_elig if matchedsample==1, r absorb(cell)
			matrix balance[`row',`col'+6]   = _b[gf_elig]
			matrix balance[`row'+1,`col'+6] = _se[gf_elig]
			matrix balance[`row'+2,`col'+6] = e(N)
		areg `dep'_bl gf_elig c_state_mrawsc4 c_state_erawsc4 if matchedsample==1, r absorb(cell)
			matrix balance[`row',`col'+7]   = _b[gf_elig]
			matrix balance[`row'+1,`col'+7] = _se[gf_elig]
			matrix balance[`row'+2,`col'+7] = e(N)
		areg `dep'_bl gf_elig c_state_mrawsc4 c_state_erawsc4 if matchedsample==1 & no_attrit0==1 & no_attrit1==1, r absorb(cell)
			matrix balance[`row',`col'+8]   = _b[gf_elig]
			matrix balance[`row'+1,`col'+8] = _se[gf_elig]
			matrix balance[`row'+2,`col'+8] = e(N)
	}
	if "`dep'" == "mgr" | "`dep'" == "egr" {
		areg `dep'_bl gf_elig if matchedsample==1, r absorb(cell)
			matrix balance[`row',`col'+6]   = _b[gf_elig]
			matrix balance[`row'+1,`col'+6] = _se[gf_elig]
			matrix balance[`row'+2,`col'+6] = e(N)
		areg `dep'_bl gf_elig c_state_mrawsc3 c_state_erawsc3 if matchedsample==1, r absorb(cell)
			matrix balance[`row',`col'+7]   = _b[gf_elig]
			matrix balance[`row'+1,`col'+7] = _se[gf_elig]
			matrix balance[`row'+2,`col'+7] = e(N)
		areg `dep'_bl gf_elig c_state_mrawsc3 c_state_erawsc3 if matchedsample==1 & no_attrit0==1 & no_attrit1==1, r absorb(cell)
			matrix balance[`row',`col'+8]   = _b[gf_elig]
			matrix balance[`row'+1,`col'+8] = _se[gf_elig]
			matrix balance[`row'+2,`col'+8] = e(N)
	}
	
	local row = `row' + 3

} 

foreach dep in no_attrit0 no_attrit1 no_attrit2 {

	summ `dep' if compsample==1
		matrix balance[`row',`col']   = r(mean)
		matrix balance[`row'+2,`col'] = r(N)
	summ `dep' if compsample==1 & compsamp_charter==1
		matrix balance[`row',`col'+1]   = r(mean)
		matrix balance[`row'+2,`col'+1] = r(N)
	summ `dep' if matchedsample==1 & any_attendUP==1
		matrix balance[`row',`col'+2]   = r(mean)
		matrix balance[`row'+2,`col'+2] = r(N)
	summ `dep' if matchedsample==1 & gf_elig==1
		matrix balance[`row',`col'+3]   = r(mean)
		matrix balance[`row'+2,`col'+3] = r(N)
	summ `dep' if matchedsample==1 & gf_elig==0
		matrix balance[`row',`col'+5]   = r(mean)
		matrix balance[`row'+2,`col'+5] = r(N)

	areg `dep' gf_elig if matchedsample==1, r absorb(cell)
		matrix balance[`row',`col'+6]   = _b[gf_elig]
		matrix balance[`row'+1,`col'+6] = _se[gf_elig]
		matrix balance[`row'+2,`col'+6] = e(N)

	local row = `row' + 3

} 

clear
svmat balance
br

}

***********************************************
* Grandfathering IV                           *
***********************************************
if `gfing' == 1 {
use "${simsdata}simsoct10", clear
keep if (grade=="06" | grade=="07") & school=="00350435" 
keep sasid grade
destring grade, replace
destring sasid, replace
gen ingavin = 1
reshape wide ingavin, i(sasid) j(grade)
gen year6 = 2011
gen year7 = 2011
tempfile gavin_fa10
save `gavin_fa10', replace
	
use if boston5==1 & spedalt5==0 & year5>=2006 & hasbaselinedemos using "${cleandata}obsfile_2014", clear 
egen cell = group(year5 masscode5 hisp5 black5 asian5 otherrace5 female5 sped5 foodfred5)

merge 1:1 sasid using `gavin_fa10', keep(master match)
drop _merge

forval g = 6/7 {
	gen gf_elig`g' = 1 if ingavin`g'==1 & year`g'==2011
	replace gf_elig`g' = 0 if ingavin`g'==. & year`g'==2011
	bys cell: egen temp = max(gf_elig`g') if cell!=.
	gen samp`g' = (temp==1 & gf_elig`g'!=.) & bcharter5==0 & spedalt5==0
	drop temp
}

replace samp7 = 0 if samp6 == 1
gen matchedsample = (samp6==1 | samp7==1)
gen gf_elig = gf_elig6 if samp6
replace gf_elig = gf_elig7 if samp7
replace cell = . if matchedsample == 0

keep if matchedsample == 1

qui tab cell, gen(d_)

foreach g in 6 7 8 {
	gen attend_UP`g' = masscode`g'  == 4800405
	replace attend_UP`g' = . if masscode`g' == .
}

gen years1 = attend_UP7 if samp6==1
replace years1 = attend_UP8 if samp7==1
gen years2 = attend_UP8+attend_UP7*(1+repeat7) if samp6==1

forval out = 1/2 {	
	local outgr10 = 6+`out'
	local outgr09 = 6+`out'+1
	foreach s in m e {
		gen `s'score`out' = c_state_`s'rawsc`outgr10' if samp6==1
		replace `s'score`out' = c_state_`s'rawsc`outgr09' if samp7==1
	}
}

expand 2, gen(copy)
qui tab copy, gen(copy_)

gen years = years1 if copy==0
replace years = years2 if copy==1

gen grade = 7 if samp6==1 & copy==0
replace grade = 8 if (samp6==1 & copy==1) | (samp7==1 & copy==0)
qui tab grade, gen(outg_)

gen yearout = year7 if grade==7
replace yearout = year8 if grade==8
qui tab yearout, gen(y_)

local row = 1
local col = 1

foreach s in m e {

		gen `s'score  = `s'score1 if copy==0
		replace `s'score = `s'score2 if copy==1
		
		summ `s'score if gf_elig==0 & ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1))
			matrix gfingnoFD[`row',`col']   = r(mean)
			matrix gfingnoFD[`row'+2,`col'] = r(N)
		
		/* OLS */
		reg `s'score years d_* `controls' copy_* outg_* if ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)), cluster(sasid)
			matrix gfingnoFD[`row',`col'+1]   = _b[years]
			matrix gfingnoFD[`row'+1,`col'+1] = _se[years]
			matrix gfingnoFD[`row'+2,`col'+1] = e(N)	

		/* RF */
		areg `s'score gf_elig `controls' copy_* outg_* if years!=. & ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)), absorb(cell) cluster(sasid)
			matrix gfingnoFD[`row',`col'+2]   = _b[gf_elig]
			matrix gfingnoFD[`row'+1,`col'+2] = _se[gf_elig]
			matrix gfingnoFD[`row'+2,`col'+2] = e(N)
		
		/* FS */
		areg years gf_elig `controls' copy_* outg_* if `s'score!=. & ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)), absorb(cell) cluster(sasid)
			matrix gfingnoFD[`row',`col'+3]   = _b[gf_elig]
			matrix gfingnoFD[`row'+1,`col'+3] = _se[gf_elig]
			matrix gfingnoFD[`row'+2,`col'+3] = e(N)
			
		/* 2SLS */
		ivreg2 `s'score (years = gf_elig) d_* `controls' copy_* outg_* if ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)), partial(d_* copy_*) cluster(sasid)
			matrix gfingnoFD[`row',`col'+4]   = _b[years]
			matrix gfingnoFD[`row'+1,`col'+4] = _se[years]
			matrix gfingnoFD[`row'+2,`col'+4] = e(N)

	local row = `row' + 3
}

forval copy = 0/1 {

	foreach s in m e {
		
		summ `s'score if gf_elig==0 & ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)) & copy==`copy'
			matrix gfingnoFD[`row',`col']   = r(mean)
			matrix gfingnoFD[`row'+2,`col'] = r(N)
			
		/* OLS */
		reg `s'score years  d_* `controls' copy_* outg_* if ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)) & copy==`copy', cluster(sasid)
			matrix gfingnoFD[`row',`col'+1]   = _b[years]
			matrix gfingnoFD[`row'+1,`col'+1] = _se[years]
			matrix gfingnoFD[`row'+2,`col'+1] = e(N)
		
		/* RF */
		areg `s'score gf_elig `controls' copy_* outg_* if years!=. & ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)) & copy==`copy', absorb(cell) cluster(sasid)
			matrix gfingnoFD[`row',`col'+2]   = _b[gf_elig]
			matrix gfingnoFD[`row'+1,`col'+2] = _se[gf_elig]
			matrix gfingnoFD[`row'+2,`col'+2] = e(N)
		
		/* FS */
		areg years gf_elig `controls' copy_* outg_* if `s'score!=. & ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)) & copy==`copy', absorb(cell) cluster(sasid)
			matrix gfingnoFD[`row',`col'+3]   = _b[gf_elig]
			matrix gfingnoFD[`row'+1,`col'+3] = _se[gf_elig]
			matrix gfingnoFD[`row'+2,`col'+3] = e(N)
			
		/* 2SLS */
		ivreg2 `s'score (years = gf_elig) d_* `controls' copy_* outg_* if ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)) & copy==`copy', partial(d_* copy_*) cluster(sasid)
			matrix gfingnoFD[`row',`col'+4]   = _b[years]
			matrix gfingnoFD[`row'+1,`col'+4] = _se[years]
			matrix gfingnoFD[`row'+2,`col'+4] = e(N)

	local row = `row' + 3
	}
}

forval grade = 7/8 {

	foreach s in m e {
		
		summ `s'score if gf_elig==0 & ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)) & grade==`grade'
			matrix gfingnoFD[`row',`col']   = r(mean)
			matrix gfingnoFD[`row'+2,`col'] = r(N)
			
		/* OLS */
		reg `s'score years d_* `controls' copy_* outg_* if ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)) & grade==`grade',  cluster(sasid)
			matrix gfingnoFD[`row',`col'+1]   = _b[years]
			matrix gfingnoFD[`row'+1,`col'+1] = _se[years]
			matrix gfingnoFD[`row'+2,`col'+1] = e(N)	
		
		/* RF */
		areg `s'score gf_elig `controls' copy_* outg_* if years!=. & ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)) & grade==`grade', absorb(cell) cluster(sasid)
			matrix gfingnoFD[`row',`col'+2]   = _b[gf_elig]
			matrix gfingnoFD[`row'+1,`col'+2] = _se[gf_elig]
			matrix gfingnoFD[`row'+2,`col'+2] = e(N)
		
		/* FS */
		areg years gf_elig `controls' copy_* outg_* if `s'score!=. & ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)) & grade==`grade', absorb(cell) cluster(sasid)
			matrix gfingnoFD[`row',`col'+3]   = _b[gf_elig]
			matrix gfingnoFD[`row'+1,`col'+3] = _se[gf_elig]
			matrix gfingnoFD[`row'+2,`col'+3] = e(N)
			
		/* 2SLS */
		ivreg2 `s'score (years = gf_elig) d_* `controls' copy_* outg_* if ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)) & grade==`grade', partial(d_* copy_*) cluster(sasid)
			matrix gfingnoFD[`row',`col'+4]   = _b[years]
			matrix gfingnoFD[`row'+1,`col'+4] = _se[years]
			matrix gfingnoFD[`row'+2,`col'+4] = e(N)

	local row = `row' + 3
	}
}

/* linear exclusion restriction fix */			
local row = 1
foreach s in m e {
		
		gen `s'level = .
		replace `s'level = `s'score
		replace `s'score = `s'score - c_state_`s'rawsc6 if samp6==1
		replace `s'score = `s'score - c_state_`s'rawsc7 if samp7==1
		
		*MS added nonmissing restriction to keep consistent across tables
		summ `s'level if gf_elig==0 & ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1))
			matrix gfing[`row',`col']   = r(mean)
			matrix gfing[`row'+2,`col'] = r(N)
			
		/* OLS */
		reg `s'score years  d_* `controls' copy_* outg_* ,  cluster(sasid)
			matrix gfing[`row',`col'+1]   = _b[years]
			matrix gfing[`row'+1,`col'+1] = _se[years]
			matrix gfing[`row'+2,`col'+1] = e(N)	
		
		
		/* RF */
		areg `s'score gf_elig `controls' copy_* outg_* if years!=., absorb(cell) cluster(sasid)
			matrix gfing[`row',`col'+2]   = _b[gf_elig]
			matrix gfing[`row'+1,`col'+2] = _se[gf_elig]
			matrix gfing[`row'+2,`col'+2] = e(N)
		
		/* FS */
		areg years gf_elig `controls' copy_* outg_* if `s'score!=., absorb(cell) cluster(sasid)
			matrix gfing[`row',`col'+3]   = _b[gf_elig]
			matrix gfing[`row'+1,`col'+3] = _se[gf_elig]
			matrix gfing[`row'+2,`col'+3] = e(N)
			
		/* 2SLS */
		ivreg2 `s'score (years = gf_elig) d_* `controls' copy_* outg_* , partial(d_* copy_*) cluster(sasid)
			matrix gfing[`row',`col'+4]   = _b[years]
			matrix gfing[`row'+1,`col'+4] = _se[years]
			matrix gfing[`row'+2,`col'+4] = e(N)

	local row = `row' + 3
}

forval copy = 0/1 {
	
	foreach s in m e {
		
		*Summarize level instead of gain
		summ `s'level if gf_elig==0 & ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)) & copy==`copy'
			matrix gfing[`row',`col']   = r(mean)
			matrix gfing[`row'+2,`col'] = r(N)
			
		/* OLS */
		reg `s'score years  d_* `controls' copy_* outg_* if ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)) & copy==`copy', cluster(sasid)
			matrix gfing[`row',`col'+1]   = _b[years]
			matrix gfing[`row'+1,`col'+1] = _se[years]
			matrix gfing[`row'+2,`col'+1] = e(N)	
		
		/* RF */
		areg `s'score gf_elig `controls' copy_* outg_* if years!=. & ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)) & copy==`copy', absorb(cell) cluster(sasid)
			matrix gfing[`row',`col'+2]   = _b[gf_elig]
			matrix gfing[`row'+1,`col'+2] = _se[gf_elig]
			matrix gfing[`row'+2,`col'+2] = e(N)
		
		/* FS */
		areg years gf_elig `controls' copy_* outg_* if `s'score!=. & ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)) & copy==`copy', absorb(cell) cluster(sasid)
			matrix gfing[`row',`col'+3]   = _b[gf_elig]
			matrix gfing[`row'+1,`col'+3] = _se[gf_elig]
			matrix gfing[`row'+2,`col'+3] = e(N)
			
		/* 2SLS */
		ivreg2 `s'score (years = gf_elig) d_* `controls' copy_* outg_* if ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)) & copy==`copy', partial(d_* copy_*) cluster(sasid)
			matrix gfing[`row',`col'+4]   = _b[years]
			matrix gfing[`row'+1,`col'+4] = _se[years]
			matrix gfing[`row'+2,`col'+4] = e(N)

	local row = `row' + 3
	}
}

forval grade = 7/8 {
	
	foreach s in m e {
		
		*Summarize level instead of scores
		summ `s'level if gf_elig==0 & ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)) & grade==`grade'
			matrix gfing[`row',`col']   = r(mean)
			matrix gfing[`row'+2,`col'] = r(N)
			
		/* OLS */
		reg `s'score years d_* `controls' copy_* outg_* if ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)) & grade==`grade',  cluster(sasid)
			matrix gfing[`row',`col'+1]   = _b[years]
			matrix gfing[`row'+1,`col'+1] = _se[years]
			matrix gfing[`row'+2,`col'+1] = e(N)
		
		/* RF */
		areg `s'score gf_elig `controls' copy_* outg_* if years!=. & ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)) & grade==`grade', absorb(cell) cluster(sasid)
			matrix gfing[`row',`col'+2]   = _b[gf_elig]
			matrix gfing[`row'+1,`col'+2] = _se[gf_elig]
			matrix gfing[`row'+2,`col'+2] = e(N)
		
		/* FS */
		areg years gf_elig `controls' copy_* outg_* if `s'score!=. & ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)) & grade==`grade', absorb(cell) cluster(sasid)
			matrix gfing[`row',`col'+3]   = _b[gf_elig]
			matrix gfing[`row'+1,`col'+3] = _se[gf_elig]
			matrix gfing[`row'+2,`col'+3] = e(N)
			
		/* 2SLS */
		ivreg2 `s'score (years = gf_elig) d_* `controls' copy_* outg_* if ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)) & grade==`grade', partial(d_* copy_*) cluster(sasid)
			matrix gfing[`row',`col'+4]   = _b[years]
			matrix gfing[`row'+1,`col'+4] = _se[years]
			matrix gfing[`row'+2,`col'+4] = e(N)

	local row = `row' + 3
	}
}

clear
svmat gfingnoFD
svmat gfing
br

}
***********************************************
* Complier school choice                      *
***********************************************
if `choice' == 1 {
use "${simsdata}simsoct10", clear
keep if (grade=="06" | grade=="07") & school=="00350435" 
keep sasid grade
destring grade, replace
destring sasid, replace
gen ingavin = 1
reshape wide ingavin, i(sasid) j(grade)
gen year6 = 2011
gen year7 = 2011
tempfile gavin_fa10
save `gavin_fa10', replace
	
use if boston5==1 & spedalt5==0 & year5>=2006 & hasbaselinedemos using "${cleandata}obsfile_2014", clear 
egen cell = group(year5 masscode5 hisp5 black5 asian5 otherrace5 female5 sped5 foodfred5)

merge 1:1 sasid using `gavin_fa10', keep(master match)
drop _merge

forval g = 6/7 {
	gen gf_elig`g' = 1 if ingavin`g'==1 & year`g'==2011
	replace gf_elig`g' = 0 if ingavin`g'==. & year`g'==2011
	bys cell: egen temp = max(gf_elig`g') if cell!=.
	gen samp`g' = (temp==1 & gf_elig`g'!=.) & bcharter5==0 & spedalt5==0
	drop temp
}

replace samp7 = 0 if samp6 == 1
gen matchedsample = (samp6==1 | samp7==1)
gen gf_elig = gf_elig6 if samp6
replace gf_elig = gf_elig7 if samp7
replace cell = . if matchedsample == 0

gen attend_UP7 = masscode7  == 4800405
replace attend_UP7 = . if masscode7 == .
gen attend_UP8 = masscode8  == 4800405
replace attend_UP8 = . if masscode8 == .

gen years1 = attend_UP7 if samp6==1
replace years1 = attend_UP8 if samp7==1

keep if matchedsample == 1
qui tab cell, gen(d_)

gen grade = 6 if samp6
replace grade = 7 if samp7
qui tab grade, gen(outg_)

gen yearout = year7 if samp6
replace yearout = year8 if samp7
qui tab yearout, gen(y_)

forval out = 1/2 {	
	local outgr10 = 6+`out'
	local outgr09 = 6+`out'+1
	foreach s in m e {
		gen `s'score`out' = c_state_`s'rawsc`outgr10' if samp6==1
		replace `s'score`out' = c_state_`s'rawsc`outgr09' if samp7==1
	}
}

gen no_attrit1 = .
replace no_attrit1 = ((mscore1!=. & c_state_mrawsc6!=.) | (c_state_erawsc7!=. & c_state_erawsc6!=.)) if samp6==1
replace no_attrit1 = ((mscore1!=. & c_state_mrawsc7!=.) | (c_state_erawsc8!=. & c_state_erawsc7!=.)) if samp7==1

gen up          = attend_UP7>0 & attend_UP7!=. & no_attrit1==1 if samp6==1
replace up      = attend_UP8>0 & attend_UP8!=. & no_attrit1==1 if samp7==1
gen bcharter     = bcharter7==1 & attend_UP7==0 & no_attrit1==1 if samp6==1
replace bcharter = bcharter8==1 & attend_UP8==0 & no_attrit1==1 if samp7==1
gen bpublic      = bcharter7==0 & boston7==1 & no_attrit1==1 if samp6==1
replace bpublic  = bcharter8==0 & boston8==1 & no_attrit1==1 if samp7==1
gen massoth      = bcharter7==0 & boston7==0 & no_attrit1==1 if samp6==1
replace massoth  = bcharter8==0 & boston8==0 & no_attrit1==1 if samp7==1
gen leftmass     = charter7==. | no_attrit1==0 if samp6==1
replace leftmass = charter8==. | no_attrit1==0 if samp7==1

gen D = attend_UP7 if samp6==1
replace D = attend_UP8 if samp7==1
gen oneminusD = 1-D
replace oneminusD = 1 if oneminusD == .
gen Z = gf_elig
keep if no_attrit1==1
foreach c of varlist D gf_elig d_* outg_* `controls' {
	drop if `c' == .
}
qui probit Z d_* `controls'
predict pZ
gen kappa_0 = (1-D)*((1-Z)-(1-pZ))/(pZ*(1-pZ))

local row = 1
foreach s in up bcharter bpublic massoth leftmass {
	summ `s' if gf_elig==0
		matrix choice[`row',1] = r(mean)
		matrix choice[`row'+1,1] = r(N)
	summ `s' if gf_elig==1
		matrix choice[`row',2] = r(mean)	
		matrix choice[`row'+1,2] = r(N)
	
	gen Y = `s'*oneminusD

	cap noisily {
		summ `s' [iweight=kappa_0]
		matrix choice[`row',3] = r(mean)
	}
	drop Y
	local ++row

}
clear
svmat choice
br
}

***********************************************
* Legacy score persistence                    *
***********************************************
if `persistence' == 1 {
use "${simsdata}simsoct10", clear
keep if (grade=="06" | grade=="07") & school=="00350435" 
keep sasid grade
destring grade, replace
destring sasid, replace
gen ingavin = 1
reshape wide ingavin, i(sasid) j(grade)
gen year6 = 2011
gen year7 = 2011
tempfile gavin_fa10
save `gavin_fa10', replace
	
use if boston5==1 & spedalt5==0 & year5>=2006 & hasbaselinedemos using "${cleandata}obsfile_2014", clear 
egen cell = group(year5 masscode5 hisp5 black5 asian5 otherrace5 female5 sped5 foodfred5)

merge 1:1 sasid using `gavin_fa10', keep(master match)
drop _merge

forval g = 6/7 {
	gen gf_elig`g' = 1 if ingavin`g'==1 & year`g'==2011
	replace gf_elig`g' = 0 if ingavin`g'==. & year`g'==2011
	bys cell: egen temp = max(gf_elig`g') if cell!=.
	gen samp`g' = (temp==1 & gf_elig`g'!=.) & bcharter5==0 & spedalt5==0
	drop temp
}

replace samp7 = 0 if samp6 == 1
gen matchedsample = (samp6==1 | samp7==1)
gen gf_elig = gf_elig6 if samp6
replace gf_elig = gf_elig7 if samp7
replace cell = . if matchedsample == 0

keep if matchedsample == 1

qui tab cell, gen(d_)

foreach g in 6 7 8 {
	gen attend_UP`g' = masscode`g'  == 4800405
	replace attend_UP`g' = . if masscode`g' == .
}

gen years1 = attend_UP7 if samp6==1
replace years1 = attend_UP8 if samp7==1
gen years2 = attend_UP8+attend_UP7*(1+repeat7) if samp6==1

forval out = 1/2 {	
	local outgr10 = 6+`out'
	local outgr09 = 6+`out'+1
	foreach s in m e {
		gen `s'score`out' = c_state_`s'rawsc`outgr10' if samp6==1
		replace `s'score`out' = c_state_`s'rawsc`outgr09' if samp7==1
	}
}

expand 2, gen(copy)
qui tab copy, gen(copy_)

gen years = years1 if copy==0
replace years = years2 if copy==1

gen grade = 7 if samp6==1 & copy==0
replace grade = 8 if (samp6==1 & copy==1) | (samp7==1 & copy==0)
qui tab grade, gen(outg_)

gen yearout = year7 if grade==7
replace yearout = year8 if grade==8
qui tab yearout, gen(y_)

local row = 1
local col = 1

egen temp = group(sped_bl masscode5) 
qui tab temp, gen(sd1_)
drop temp
egen temp = group(samp6 samp7)
qui tab temp, gen(sd2_)
drop temp
qui tab year_bl, gen(sd3_)
foreach var1 of varlist sd*_* {
	foreach var2 of varlist gf_elig {
		gen `var2'_`var1' = `var2'*`var1'
	}
}

foreach s in m e {

		gen `s'score  = `s'score1 if copy==0
		replace `s'score = `s'score2 if copy==1
		gen `s'legacy = c_state_`s'rawsc6 if samp6==1
		replace `s'legacy = c_state_`s'rawsc7 if samp7==1
		gen `s'gain = `s'score - `s'legacy

		ivreg2 `s'legacy gf_elig d_* `controls' copy_* outg_* if `s'score!=. & copy==0, partial(d_* copy_* outg_*) cluster(sasid)
			matrix persistence[`row',`col']    = _b[gf_elig]
			matrix persistence[`row'+1,`col']  = _se[gf_elig]
			matrix persistence[`row'+11,`col'] = e(N)

		local ++col
		ivreg2 `s'score (years `s'legacy = gf_elig_sd*_*) sd*_* d_* `controls' copy_* outg_* if `s'legacy!=., partial(d_* sd*_* copy_* outg_*) cluster(sasid) savefirst first
			matrix temp = e(first)
			matrix persistence[`row'+2,`col']  = _b[years]
			matrix persistence[`row'+3,`col']  = _se[years]
			matrix persistence[`row'+4,`col']  = temp[7,1]
			matrix persistence[`row'+5,`col']  = _b[`s'legacy]
			matrix persistence[`row'+6,`col']  = _se[`s'legacy]
			matrix persistence[`row'+7,`col']  = temp[7,2]
			matrix persistence[`row'+10,`col'] = e(exexog_ct)-e(inexog_ct)
			matrix persistence[`row'+11,`col'] = e(N)	
		
		local ++col
		ivreg2 `s'gain (years = gf_elig_sd*_*) sd*_* d_* `controls' copy_* outg_* if `s'legacy!=., partial(d_* sd*_* copy_* outg_*) cluster(sasid) savefirst first
			matrix temp = e(first)
			matrix persistence[`row'+2,`col']   = _b[years]
			matrix persistence[`row'+3,`col'] = _se[years]
			matrix persistence[`row'+4,`col'] = temp[7,1]
			matrix persistence[`row'+10,`col'] = e(exexog_ct)-e(inexog_ct)
			matrix persistence[`row'+11,`col'] = e(N)

		preserve
		keep if `s'legacy!=.
		gen one = 1
		expandcl 2, cluster(one) gen(est)
		foreach var of varlist years `s'legacy gf_elig_sd*_* sd*_* d_* `controls' copy_* outg_* {
			forval est = 1/2 {
				gen est`est'_`var' = `var'*(est==`est')
			}
			drop `var'
		}
		gen Y = `s'score if est==1
		replace Y = `s'gain if est==2
		ivregress 2sls Y (est*_years est1_`s'legacy = est*_gf_elig_sd*_*) est1_sd*_* est2_sd*_* est*_d_* est*_foodfred_bl est*_lep_bl est*_c_state_mrawsc_bl est*_c_state_erawsc_bl est*_y_* est*_copy_* est*_outg_*
		lincom est1_years-est2_years
		matrix persistence[`row'+8,2] = r(estimate)
		matrix persistence[`row'+9,2] = r(se)
		restore
		
	local row = `row' + 12
	local col = 1
}

clear
svmat persistence
br

}
***********************************************
* Grandfathering IV with Matching on BL Score *
***********************************************
if `altgfing' == 1 {
use "${simsdata}simsoct10", clear
keep if (grade=="06" | grade=="07") & school=="00350435" 
keep sasid grade
destring grade, replace
destring sasid, replace
gen ingavin = 1
reshape wide ingavin, i(sasid) j(grade)
gen year6 = 2011
gen year7 = 2011
tempfile gavin_fa10
save `gavin_fa10', replace
	
use if boston5==1 & spedalt5==0 & year5>=2006 & hasbaselinedemos using "${cleandata}obsfile_2014", clear 
egen combscore_bl = rowmean(c_state_mrawsc5 c_state_erawsc5)
xtile combtile_bl = combscore_bl, nq(3)
egen cell = group(year5 masscode5 hisp5 black5 asian5 otherrace5 female5 sped5 foodfred5 combtile_bl)

merge 1:1 sasid using `gavin_fa10', keep(master match)
drop _merge

forval g = 6/7 {
	gen gf_elig`g' = 1 if ingavin`g'==1 & year`g'==2011
	replace gf_elig`g' = 0 if ingavin`g'==. & year`g'==2011
	bys cell: egen temp = max(gf_elig`g') if cell!=.
	gen samp`g' = (temp==1 & gf_elig`g'!=.) & bcharter5==0 & spedalt5==0
	drop temp
}

replace samp7 = 0 if samp6 == 1
gen matchedsample = (samp6==1 | samp7==1)
gen gf_elig = gf_elig6 if samp6
replace gf_elig = gf_elig7 if samp7
replace cell = . if matchedsample == 0

keep if matchedsample == 1

qui tab cell, gen(d_)

foreach g in 6 7 8 {
	gen attend_UP`g' = masscode`g'  == 4800405
	replace attend_UP`g' = . if masscode`g' == .
}

gen years1 = attend_UP7 if samp6==1
replace years1 = attend_UP8 if samp7==1
gen years2 = attend_UP8+attend_UP7*(1+repeat7) if samp6==1

forval out = 1/2 {	
	local outgr10 = 6+`out'
	local outgr09 = 6+`out'+1
	foreach s in m e {
		gen `s'score`out' = c_state_`s'rawsc`outgr10' if samp6==1
		replace `s'score`out' = c_state_`s'rawsc`outgr09' if samp7==1
	}
}

expand 2, gen(copy)
qui tab copy, gen(copy_)

gen years = years1 if copy==0
replace years = years2 if copy==1

gen grade = 7 if samp6==1 & copy==0
replace grade = 8 if (samp6==1 & copy==1) | (samp7==1 & copy==0)
qui tab grade, gen(outg_)

gen yearout = year7 if grade==7
replace yearout = year8 if grade==8
qui tab yearout, gen(y_)

local row = 1
local col = 1

foreach s in m e {

		gen `s'score  = `s'score1 if copy==0
		replace `s'score = `s'score2 if copy==1
		
		summ `s'score if gf_elig==0 & ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1))
			matrix gfingnoFD[`row',`col']   = r(mean)
			matrix gfingnoFD[`row'+2,`col'] = r(N)
		
		/* OLS */
		reg `s'score years d_* `controls' copy_* outg_* if ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)), cluster(sasid)
			matrix gfingnoFD[`row',`col'+1]   = _b[years]
			matrix gfingnoFD[`row'+1,`col'+1] = _se[years]
			matrix gfingnoFD[`row'+2,`col'+1] = e(N)	

		/* RF */
		areg `s'score gf_elig `controls' copy_* outg_* if years!=. & ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)), absorb(cell) cluster(sasid)
			matrix gfingnoFD[`row',`col'+2]   = _b[gf_elig]
			matrix gfingnoFD[`row'+1,`col'+2] = _se[gf_elig]
			matrix gfingnoFD[`row'+2,`col'+2] = e(N)
		
		/* FS */
		areg years gf_elig `controls' copy_* outg_* if `s'score!=. & ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)), absorb(cell) cluster(sasid)
			matrix gfingnoFD[`row',`col'+3]   = _b[gf_elig]
			matrix gfingnoFD[`row'+1,`col'+3] = _se[gf_elig]
			matrix gfingnoFD[`row'+2,`col'+3] = e(N)
			
		/* 2SLS */
		ivreg2 `s'score (years = gf_elig) d_* `controls' copy_* outg_* if ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)), partial(d_* copy_*) cluster(sasid)
			matrix gfingnoFD[`row',`col'+4]   = _b[years]
			matrix gfingnoFD[`row'+1,`col'+4] = _se[years]
			matrix gfingnoFD[`row'+2,`col'+4] = e(N)

	local row = `row' + 3
}

forval copy = 0/1 {

	foreach s in m e {
		
		summ `s'score if gf_elig==0 & ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)) & copy==`copy'
			matrix gfingnoFD[`row',`col']   = r(mean)
			matrix gfingnoFD[`row'+2,`col'] = r(N)
			
		/* OLS */
		reg `s'score years  d_* `controls' copy_* outg_* if ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)) & copy==`copy', cluster(sasid)
			matrix gfingnoFD[`row',`col'+1]   = _b[years]
			matrix gfingnoFD[`row'+1,`col'+1] = _se[years]
			matrix gfingnoFD[`row'+2,`col'+1] = e(N)
		
		/* RF */
		areg `s'score gf_elig `controls' copy_* outg_* if years!=. & ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)) & copy==`copy', absorb(cell) cluster(sasid)
			matrix gfingnoFD[`row',`col'+2]   = _b[gf_elig]
			matrix gfingnoFD[`row'+1,`col'+2] = _se[gf_elig]
			matrix gfingnoFD[`row'+2,`col'+2] = e(N)
		
		/* FS */
		areg years gf_elig `controls' copy_* outg_* if `s'score!=. & ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)) & copy==`copy', absorb(cell) cluster(sasid)
			matrix gfingnoFD[`row',`col'+3]   = _b[gf_elig]
			matrix gfingnoFD[`row'+1,`col'+3] = _se[gf_elig]
			matrix gfingnoFD[`row'+2,`col'+3] = e(N)
			
		/* 2SLS */
		ivreg2 `s'score (years = gf_elig) d_* `controls' copy_* outg_* if ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)) & copy==`copy', partial(d_* copy_*) cluster(sasid)
			matrix gfingnoFD[`row',`col'+4]   = _b[years]
			matrix gfingnoFD[`row'+1,`col'+4] = _se[years]
			matrix gfingnoFD[`row'+2,`col'+4] = e(N)

	local row = `row' + 3
	}
}

forval grade = 7/8 {

	foreach s in m e {
		
		summ `s'score if gf_elig==0 & ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)) & grade==`grade'
			matrix gfingnoFD[`row',`col']   = r(mean)
			matrix gfingnoFD[`row'+2,`col'] = r(N)
			
		/* OLS */
		reg `s'score years d_* `controls' copy_* outg_* if ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)) & grade==`grade',  cluster(sasid)
			matrix gfingnoFD[`row',`col'+1]   = _b[years]
			matrix gfingnoFD[`row'+1,`col'+1] = _se[years]
			matrix gfingnoFD[`row'+2,`col'+1] = e(N)	
		
		/* RF */
		areg `s'score gf_elig `controls' copy_* outg_* if years!=. & ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)) & grade==`grade', absorb(cell) cluster(sasid)
			matrix gfingnoFD[`row',`col'+2]   = _b[gf_elig]
			matrix gfingnoFD[`row'+1,`col'+2] = _se[gf_elig]
			matrix gfingnoFD[`row'+2,`col'+2] = e(N)
		
		/* FS */
		areg years gf_elig `controls' copy_* outg_* if `s'score!=. & ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)) & grade==`grade', absorb(cell) cluster(sasid)
			matrix gfingnoFD[`row',`col'+3]   = _b[gf_elig]
			matrix gfingnoFD[`row'+1,`col'+3] = _se[gf_elig]
			matrix gfingnoFD[`row'+2,`col'+3] = e(N)
			
		/* 2SLS */
		ivreg2 `s'score (years = gf_elig) d_* `controls' copy_* outg_* if ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)) & grade==`grade', partial(d_* copy_*) cluster(sasid)
			matrix gfingnoFD[`row',`col'+4]   = _b[years]
			matrix gfingnoFD[`row'+1,`col'+4] = _se[years]
			matrix gfingnoFD[`row'+2,`col'+4] = e(N)

	local row = `row' + 3
	}
}

/* linear exclusion restriction fix */			
local row = 1
foreach s in m e {
		
		gen `s'level = .
		replace `s'level = `s'score
		replace `s'score = `s'score - c_state_`s'rawsc6 if samp6==1
		replace `s'score = `s'score - c_state_`s'rawsc7 if samp7==1
		
		*MS added nonmissing restriction to keep consistent across tables
		summ `s'level if gf_elig==0 & ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1))
			matrix gfing[`row',`col']   = r(mean)
			matrix gfing[`row'+2,`col'] = r(N)
			
		/* OLS */
		reg `s'score years  d_* `controls' copy_* outg_* ,  cluster(sasid)
			matrix gfing[`row',`col'+1]   = _b[years]
			matrix gfing[`row'+1,`col'+1] = _se[years]
			matrix gfing[`row'+2,`col'+1] = e(N)	
		
		
		/* RF */
		areg `s'score gf_elig `controls' copy_* outg_* if years!=., absorb(cell) cluster(sasid)
			matrix gfing[`row',`col'+2]   = _b[gf_elig]
			matrix gfing[`row'+1,`col'+2] = _se[gf_elig]
			matrix gfing[`row'+2,`col'+2] = e(N)
		
		/* FS */
		areg years gf_elig `controls' copy_* outg_* if `s'score!=., absorb(cell) cluster(sasid)
			matrix gfing[`row',`col'+3]   = _b[gf_elig]
			matrix gfing[`row'+1,`col'+3] = _se[gf_elig]
			matrix gfing[`row'+2,`col'+3] = e(N)
			
		/* 2SLS */
		ivreg2 `s'score (years = gf_elig) d_* `controls' copy_* outg_* , partial(d_* copy_*) cluster(sasid)
			matrix gfing[`row',`col'+4]   = _b[years]
			matrix gfing[`row'+1,`col'+4] = _se[years]
			matrix gfing[`row'+2,`col'+4] = e(N)

	local row = `row' + 3
}

forval copy = 0/1 {
	
	foreach s in m e {
		
		*Summarize level instead of gain
		summ `s'level if gf_elig==0 & ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)) & copy==`copy'
			matrix gfing[`row',`col']   = r(mean)
			matrix gfing[`row'+2,`col'] = r(N)
			
		/* OLS */
		reg `s'score years  d_* `controls' copy_* outg_* if ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)) & copy==`copy', cluster(sasid)
			matrix gfing[`row',`col'+1]   = _b[years]
			matrix gfing[`row'+1,`col'+1] = _se[years]
			matrix gfing[`row'+2,`col'+1] = e(N)	
		
		/* RF */
		areg `s'score gf_elig `controls' copy_* outg_* if years!=. & ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)) & copy==`copy', absorb(cell) cluster(sasid)
			matrix gfing[`row',`col'+2]   = _b[gf_elig]
			matrix gfing[`row'+1,`col'+2] = _se[gf_elig]
			matrix gfing[`row'+2,`col'+2] = e(N)
		
		/* FS */
		areg years gf_elig `controls' copy_* outg_* if `s'score!=. & ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)) & copy==`copy', absorb(cell) cluster(sasid)
			matrix gfing[`row',`col'+3]   = _b[gf_elig]
			matrix gfing[`row'+1,`col'+3] = _se[gf_elig]
			matrix gfing[`row'+2,`col'+3] = e(N)
			
		/* 2SLS */
		ivreg2 `s'score (years = gf_elig) d_* `controls' copy_* outg_* if ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)) & copy==`copy', partial(d_* copy_*) cluster(sasid)
			matrix gfing[`row',`col'+4]   = _b[years]
			matrix gfing[`row'+1,`col'+4] = _se[years]
			matrix gfing[`row'+2,`col'+4] = e(N)

	local row = `row' + 3
	}
}

forval grade = 7/8 {
	
	foreach s in m e {
		
		*Summarize level instead of scores
		summ `s'level if gf_elig==0 & ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)) & grade==`grade'
			matrix gfing[`row',`col']   = r(mean)
			matrix gfing[`row'+2,`col'] = r(N)
			
		/* OLS */
		reg `s'score years d_* `controls' copy_* outg_* if ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)) & grade==`grade',  cluster(sasid)
			matrix gfing[`row',`col'+1]   = _b[years]
			matrix gfing[`row'+1,`col'+1] = _se[years]
			matrix gfing[`row'+2,`col'+1] = e(N)
		
		/* RF */
		areg `s'score gf_elig `controls' copy_* outg_* if years!=. & ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)) & grade==`grade', absorb(cell) cluster(sasid)
			matrix gfing[`row',`col'+2]   = _b[gf_elig]
			matrix gfing[`row'+1,`col'+2] = _se[gf_elig]
			matrix gfing[`row'+2,`col'+2] = e(N)
		
		/* FS */
		areg years gf_elig `controls' copy_* outg_* if `s'score!=. & ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)) & grade==`grade', absorb(cell) cluster(sasid)
			matrix gfing[`row',`col'+3]   = _b[gf_elig]
			matrix gfing[`row'+1,`col'+3] = _se[gf_elig]
			matrix gfing[`row'+2,`col'+3] = e(N)
			
		/* 2SLS */
		ivreg2 `s'score (years = gf_elig) d_* `controls' copy_* outg_* if ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1)) & grade==`grade', partial(d_* copy_*) cluster(sasid)
			matrix gfing[`row',`col'+4]   = _b[years]
			matrix gfing[`row'+1,`col'+4] = _se[years]
			matrix gfing[`row'+2,`col'+4] = e(N)

	local row = `row' + 3
	}
}

clear
svmat gfingnoFD
svmat gfing
br

}
***********************************************
* Peer composition                            *
***********************************************
if `switching' == 1 {
	use "${simsdata}simsoct10", clear
	keep if (grade=="06" | grade=="07") & school=="00350435" 
	keep sasid grade
	destring grade, replace
	destring sasid, replace
	gen ingavin = 1
	reshape wide ingavin, i(sasid) j(grade)
	gen year6 = 2011
	gen year7 = 2011
	tempfile gavin_fa10
	save `gavin_fa10', replace
			
	use if boston5==1 & spedalt5==0 & year5>=2006 & hasbaselinedemos using "${cleandata}obsfile_2014_wpeers", clear 
	egen cell = group(year5 masscode5 hisp5 black5 asian5 otherrace5 female5 sped5 foodfred5)

	merge 1:1 sasid using `gavin_fa10', keep(master match)
	drop _merge

	forval g = 6/7 {
		gen gf_elig`g' = 1 if ingavin`g'==1 & year`g'==2011
		replace gf_elig`g' = 0 if ingavin`g'==. & year`g'==2011
		bys cell: egen temp = max(gf_elig`g') if cell!=.
		gen samp`g' = (temp==1 & gf_elig`g'!=.) & bcharter5==0 & spedalt5==0
		drop temp
	}

	replace samp7 = 0 if samp6 == 1
	gen matchedsample = (samp6==1 | samp7==1)
	gen gf_elig = gf_elig6 if samp6
	replace gf_elig = gf_elig7 if samp7
	replace cell = . if matchedsample == 0

	keep if matchedsample == 1
	qui tab cell, gen(d_)
	gen grade = 7 if samp6==1
	replace grade = 8 if samp7==1 
	qui tab grade, gen(outg_)

	gen yearout = year7 if grade==7
	replace yearout = year8 if grade==8
	qui tab yearout, gen(y_)

	foreach t in mscore escore sped lep foodfred {
		gen peer_`t'1 = peer_`t'_bl5_for7 if samp6==1
		replace peer_`t'1 = peer_`t'_bl5_for8 if samp7==1
	}
	foreach g in 7 8 {
		gen attend_UP`g' = masscode`g'  == 4800405
		replace attend_UP`g' = . if masscode`g' == .
	}	
	gen attend_1 = attend_UP7 if samp6==1
	replace attend_1 = attend_UP8 if samp7==1
	
	local row = 1
	local col = 1	
	
	foreach t in sped foodfred lep mscore escore {
		summ peer_`t'1 if gf_elig==0
			matrix switching[`row',`col']   = r(mean)
			matrix switching[`row'+2,`col'] = r(N)
			
		ivreg2 peer_`t'1 (attend_1 = gf_elig) d_* outg_* `controls', r partial(d_* `controls')
			matrix switching[`row',`col'+1]   = _b[attend_1]
			matrix switching[`row'+1,`col'+1] = _se[attend_1]
			matrix switching[`row'+2,`col'+1] = e(N)
		local row = `row' + 3
	}
	
	clear
	svmat switching
	br
}
