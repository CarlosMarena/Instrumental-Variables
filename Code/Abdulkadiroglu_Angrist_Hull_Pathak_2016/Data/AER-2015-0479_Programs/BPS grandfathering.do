* P.HULL
* Dearborn, Harbor, and Orchard Gardens GF-IV

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

matrix cohorts             = J(200,10,.)
matrix dearharbbalance     = J(200,10,.)
matrix OGbalance           = J(200,10,.)
matrix OGgfing             = J(200,5,.)
matrix OGgfingnoFD         = J(200,5,.)
matrix harbdeargfing       = J(200,5,.)
matrix harbdeargfingnoFD   = J(200,5,.)
matrix switching           = J(200,10,.)

local cohorts   = 1 /* Table B2 (Panel B)     */
local balance   = 1 /* Table B11              */

local fdiv		= 1 /* Table 9                */
local nofdiv	= 1 
local switching = 1 /* Table 10 (Columns 3-4) */

local gfing = 0
if `fdiv' == 1 | `nofdiv' == 1 {
	local gfing = 1
}

***********************************************
* Counting match & follow-up                  *
***********************************************
if `cohorts' == 1 {
	use "${simsdata}simsoct09", clear
	keep if (grade=="06" | grade=="07") & (school=="00350257" | school=="00350426" | school=="00350074")
	gen inOG   = (school=="00350257")
	gen inharb = (school=="00350426")
	gen indear = (school=="00350074")
	keep sasid grade inOG inharb indear
	destring grade, replace
	destring sasid, replace
	reshape wide inOG inharb indear, i(sasid) j(grade)
	gen year6 = 2010
	gen year7 = 2010
	tempfile legacy_fa09
	save `legacy_fa09', replace

	use "${cleandata}obsfile_2014", clear 
	drop if (boston5 == 0 & boston6 == 0) | year5<2006
	
	merge 1:1 sasid using `legacy_fa09', keep(master match)
	drop _merge

	preserve
	egen harbcell = group(year5 masscode5 hisp5 black5 asian5 otherrace5 female5 sped5 foodfred5)
	gen dearcell = harbcell

	foreach sch in harb dear {
		forval g = 6/7 {
			gen `sch'_elig`g' = 1 if in`sch'`g'==1 & year`g'==2010
			replace `sch'_elig`g' = 0 if in`sch'`g'==. & year`g'==2010
			bys `sch'cell: egen temp = max(`sch'_elig`g') if `sch'cell!=.
			gen `sch'samp`g' = (temp==1 & `sch'_elig`g'!=.) & boston5==1 & bcharter5==0 & spedalt5==0
			drop temp
		}
		replace `sch'samp7 = 0 if `sch'samp6 == 1
		gen `sch'matchedsample = (`sch'samp6==1 | `sch'samp7==1)
	}

	keep if harb_elig6==1 | harb_elig7==1 | harbmatchedsample | dear_elig6==1 | dear_elig7==1 | dearmatchedsample
	drop *_bl
	foreach var in `covstems' white boston bcharter spedalt c_state_mrawsc c_state_erawsc year {
		gen `var'_bl  = `var'5 
		gen `var'_pbl = `var'4 
	}

	tempfile harbdear 
	save `harbdear', replace

	restore
	drop *_bl
	foreach var in `covstems' white boston bcharter spedalt c_state_mrawsc c_state_erawsc year {
		gen `var'_bl = .
		gen `var'_pbl = .
	}
	forval bl = 5/6 {

		gen sch_lettergrade_adj`bl' = .
		bys year`bl' masscode`bl': egen sch_boston`bl' = mean(boston`bl')
		gen temp = c_boston_mrawsc`bl' + c_boston_erawsc`bl' if sch_boston`bl' == 1
		bys year`bl' masscode`bl': egen sch_rawsc`bl' = mean(temp) if sch_boston`bl'==1
		drop temp
		bys year`bl' masscode`bl': gen temp = _n==1
		levelsof year`bl', local(years)
		foreach y in `years' {
			cap {
				pctile cut = sch_rawsc`bl' if temp==1 & year`bl'==`y', n(10)
				xtile bin = sch_rawsc`bl', cutpoint(cut)
				replace sch_lettergrade_adj`bl' = bin if year`bl'==`y'
				drop cut bin
			}
		}
		drop temp

		egen OGcell`bl' = group(year`bl' sch_lettergrade_adj`bl' hisp`bl' black`bl' asian`bl' otherrace`bl' female`bl' sped`bl' foodfred`bl')

		local g = `bl' + 1
		gen OG_elig`g' = 1 if inOG`g'==1 & year`g'==2010
		replace OG_elig`g' = 0 if inOG`g'==. & year`g'==2010
		bys OGcell`bl': egen temp = max(OG_elig`g') if OGcell`bl'!=.
		gen OGsamp`g' = (temp==1 & OG_elig`g'!=.) & boston`bl'==1 & bcharter`bl'==0 & spedalt`bl'==0
		replace OGcell`bl' = . if OGsamp`g'==0

		foreach var in `covstems' white boston bcharter spedalt c_state_mrawsc c_state_erawsc year {
			replace `var'_bl = `var'`bl' if OGsamp`g'==1
			local pbl = `bl' -1
			replace `var'_pbl = `var'`pbl' if OGsamp`g'==1
		}	

		drop sch_rawsc`bl' temp
	}

	egen OGcell = group(OGcell5 OGcell6), m
	gen OGmatchedsample = (OGsamp6==1 | OGsamp7==1)
	gen OG_elig = OG_elig6 if OGsamp6
	replace OG_elig = OG_elig7 if OGsamp7
	replace OGcell = . if OGmatchedsample == 0
	append using `harbdear'

	foreach sch in dear harb OG {
		egen `sch'_eligR  = rowmax(`sch'_elig6 `sch'_elig7)
		egen `sch'_eligwB = rowmax(`sch'_elig6 `sch'_elig7) if `sch'cell!=.
		egen `sch'_eligiD = rowmax(`sch'_elig6 `sch'_elig7) if `sch'cell!=. & boston_bl==1 & bcharter_bl==0 & spedalt_bl==0
		egen `sch'_eligiS = rowmax(`sch'_elig6 `sch'_elig7) if `sch'matchedsample==1
	}
	
	local row = 1 
	local col = 1
	
	foreach sch in dear harb OG {
		count if `sch'_eligR == 1
			matrix cohorts[`row',`col'] = r(N)
		local ++col
		count if `sch'_eligwB == 1
			matrix cohorts[`row',`col'] = r(N)
		local ++col	
		count if `sch'_eligiD == 1
			matrix cohorts[`row',`col'] = r(N)
		local ++col	
		count if `sch'_eligiS == 1
			matrix cohorts[`row',`col'] = r(N)
		local ++col
		count if `sch'matchedsample == 1 & `sch'_eligiS == 0
			matrix cohorts[`row',`col'] = r(N)
		local col = 1
		local ++row
	}
	clear
	svmat cohorts
	br
	
}

***********************************************
* Grandfathering Balance                      *
***********************************************
if `balance' == 1 {

use "${simsdata}simsoct09", clear
keep if (grade=="06" | grade=="07") & (school=="00350257" | school=="00350426" | school=="00350074")
gen inOG   = (school=="00350257")
gen indearharb = (school=="00350074" | school=="00350426")
keep sasid grade inOG indearharb
destring grade, replace
destring sasid, replace
reshape wide inOG indearharb, i(sasid) j(grade)
gen year6 = 2010
gen year7 = 2010
tempfile legacy_fa09
save `legacy_fa09', replace

use if (boston5==1 | boston6==1) & year5>=2006 using "${cleandata}obsfile_2014", clear 

merge 1:1 sasid using `legacy_fa09', keep(master match)
drop _merge

preserve
egen dearharbcell = group(year5 masscode5 hisp5 black5 asian5 otherrace5 female5 sped5 foodfred5)

foreach sch in dearharb {
	forval g = 6/7 {
		gen `sch'_elig`g' = 1 if in`sch'`g'==1 & year`g'==2010
		replace `sch'_elig`g' = 0 if in`sch'`g'==. & year`g'==2010
		bys `sch'cell: egen temp = max(`sch'_elig`g') if `sch'cell!=.
		gen `sch'samp`g' = (temp==1 & `sch'_elig`g'!=.) & boston5==1 & bcharter5==0 & spedalt5==0
		drop temp
	}

	replace `sch'samp7 = 0 if `sch'samp6 == 1
	gen `sch'matchedsample = (`sch'samp6==1 | `sch'samp7==1)
	gen `sch'_elig = `sch'_elig6 if `sch'samp6
	replace `sch'_elig = `sch'_elig7 if `sch'samp7
	replace `sch'cell = . if `sch'matchedsample == 0
}
keep if dearharbmatchedsample

drop *_bl
foreach var in `covstems' white boston bcharter spedalt c_state_mrawsc c_state_erawsc year {
	gen `var'_bl  = `var'5
	gen `var'_pbl = `var'4
	gen `var'_ppbl = `var'3 	
}

tempfile harbdear 
save `harbdear', replace

restore
drop *_bl
foreach var in `covstems' white boston bcharter spedalt c_state_mrawsc c_state_erawsc year {
	gen `var'_bl = .
	gen `var'_pbl = .
	gen `var'_ppbl = .
}
forval bl = 5/6 {

	gen sch_lettergrade_adj`bl' = .
	bys year`bl' masscode`bl': egen sch_boston`bl' = mean(boston`bl')
	gen temp = c_boston_mrawsc`bl' + c_boston_erawsc`bl' if sch_boston`bl' == 1
	bys year`bl' masscode`bl': egen sch_rawsc`bl' = mean(temp) if sch_boston`bl'==1
	drop temp
	bys year`bl' masscode`bl': gen temp = _n==1
	levelsof year`bl', local(years)
	foreach y in `years' {
		cap {
			pctile cut = sch_rawsc`bl' if temp==1 & year`bl'==`y', n(10)
			xtile bin = sch_rawsc`bl', cutpoint(cut)
			replace sch_lettergrade_adj`bl' = bin if year`bl'==`y'
			drop cut bin
		}
	}
	drop temp

	egen OGcell`bl' = group(year`bl' sch_lettergrade_adj`bl' hisp`bl' black`bl' asian`bl' otherrace`bl' female`bl' sped`bl' foodfred`bl')

	local g = `bl' + 1
	gen OG_elig`g' = 1 if inOG`g'==1 & year`g'==2010
	replace OG_elig`g' = 0 if inOG`g'==. & year`g'==2010
	bys OGcell`bl': egen temp = max(OG_elig`g') if OGcell`bl'!=.
	gen OGsamp`g' = (temp==1 & OG_elig`g'!=.) & boston`bl'==1 & bcharter`bl'==0 & spedalt`bl'==0
	replace OGcell`bl' = . if OGsamp`g'==0

	foreach var in `covstems' white boston bcharter spedalt c_state_mrawsc c_state_erawsc year {
		replace `var'_bl = `var'`bl' if OGsamp`g'==1
		local pbl = `bl' -1
		replace `var'_pbl = `var'`pbl' if OGsamp`g'==1
		local ppbl = `bl' - 2
		replace `var'_ppbl = `var'`ppbl' if OGsamp`g'==1		
	}	

	drop sch_rawsc`bl' temp
}

egen OGcell = group(OGcell5 OGcell6), m
gen OGmatchedsample = (OGsamp6==1 | OGsamp7==1)
gen OG_elig = OG_elig6 if OGsamp6
replace OG_elig = OG_elig7 if OGsamp7
replace OGcell = . if OGmatchedsample == 0
append using `harbdear'

foreach v in _elig matchedsample cell samp6 samp7 {
	replace dearharb`v' = 0 if dearharb`v' == .
	replace OG`v' = 0 if OG`v' == .
}

gen anymatchedsample = OGmatchedsample | dearharbmatchedsample
gen samp6 = OGsamp6 | dearharbsamp6
gen samp7 = OGsamp7 | dearharbsamp7 
gen any_elig = dearharb_elig | OG_elig==1
egen anycell = group(dearharbcell OGcell)

forval bl = 5/6 {
local g = `bl' + 1
gen compsamp`g' = 0
qui levelsof year`bl' if samp`g'==1, local(years)
qui count if anymatchedsample==1
local denom = r(N)
local prps ""
foreach y in `years' {
	qui count if year`bl'==`y' & samp`g'==1
	local prop`y' = r(N)/`denom'
	local prps "`prps', `prop`y''"
}
local prps = substr("`prps'",2,.)
foreach y in `years' {
	replace compsamp`g' = runiform()<`prop`y'' if year`bl'==`y' & compsamp`g' == 0 & boston`bl'==1 & spedalt`bl'==0
}
}
egen compsample = rowtotal(compsamp6 compsamp7)
replace compsample = 0 if compsample>1
foreach var in `covstems' white boston bcharter spedalt c_state_mrawsc c_state_erawsc year {
gen `var'_compbl = .
gen `var'_comppbl = .
}
forval bl = 5/6 {
local g = `bl' + 1
local pblg = `bl' - 1
foreach var in `covstems' white boston bcharter spedalt c_state_mrawsc c_state_erawsc year {
	replace `var'_compbl = `var'`bl' if compsamp`g' == 1
	replace `var'_comppbl = `var'`pbl' if compsamp`g' == 1
}
}

gen compsamp_charter = (charter6==1 | charter7==1 | charter8==1) if (charter6!=. | charter7!=. | charter8!=.) & compsamp6==1
replace compsamp_charter = (charter7==1 | charter8==1) if (charter7!=. | charter8!=.) & compsamp7==1

keep if anymatchedsample == 1 | compsample == 1

foreach g in 7 8 {
	gen attend_OG`g' = (masscode`g' == 350257)
	replace attend_OG`g' = . if masscode`g' == .
	gen attend_dearharb`g' = (masscode`g' == 350074 | masscode`g'  == 350426)
	replace attend_dearharb`g' = . if masscode`g' == .
	gen attend_any`g' = (masscode`g'  == 350257 | masscode`g'  == 350426 | masscode`g' == 350074)
	replace attend_any`g' = . if masscode`g' == .
}

gen any_attendOG   = (attend_OG7==1 | attend_OG8==1) if (attend_OG7!=. | attend_OG8!=.) & samp6==1
gen any_attenddearharb = (attend_dearharb7==1 | attend_dearharb8==1) if (attend_dearharb7!=. | attend_dearharb8!=.) & samp6==1
gen any_attendany  = (attend_any7==1 | attend_any8==1) if (attend_any7!=. | attend_any8!=.) & samp6==1
replace any_attendOG   = (attend_OG8==1) if (attend_OG8!=.) & samp7==1
replace any_attenddearharb = (attend_dearharb8==1) if (attend_dearharb8!=.) & samp7==1
replace any_attendany  = (attend_any8==1) if (attend_any8!=.) & samp7==1 

gen no_attrit0 = .
replace no_attrit0 = (c_state_mrawsc6!=. | c_state_erawsc6!=.) if samp6==1 | (compsamp6==1 & samp7==0)
replace no_attrit0 = (c_state_mrawsc7!=. | c_state_erawsc7!=.) if samp7==1 | (compsamp7==1 & samp6==0)
gen no_attrit1 = .
replace no_attrit1 = (c_state_mrawsc7!=. | c_state_erawsc7!=.) if samp6==1 | (compsamp6==1 & samp7==0)
replace no_attrit1 = (c_state_mrawsc8!=. | c_state_erawsc8!=.) if samp7==1 | (compsamp7==1 & samp6==0)
gen no_attrit2 = .
replace no_attrit2 = (c_state_mrawsc8!=. | c_state_erawsc8!=.) if samp6==1 | (compsamp6==1 & samp7==0)

gen gr_m_compbl = c_state_mrawsc_compbl - c_state_mrawsc_comppbl
gen gr_e_compbl = c_state_erawsc_compbl - c_state_erawsc_comppbl
gen gr_m_bl = c_state_mrawsc_bl - c_state_mrawsc_pbl
gen gr_e_bl = c_state_erawsc_bl - c_state_erawsc_pbl

foreach sch in dearharb OG {
	
	local row = 1
	local col = 1

	foreach dep in hisp black white asian female sped foodfred lep c_state_mrawsc gr_m c_state_erawsc gr_e {

	summ `dep'_compbl if compsample==1
		matrix `sch'balance[`row',`col']   = r(mean)
		matrix `sch'balance[`row'+2,`col'] = r(N)
	summ `dep'_compbl if compsample==1 & compsamp_charter==1
		matrix `sch'balance[`row',`col'+1]   = r(mean)
		matrix `sch'balance[`row'+2,`col'+1] = r(N)
	summ `dep'_bl if `sch'matchedsample==1 & any_attend`sch'==1
		matrix `sch'balance[`row',`col'+2]   = r(mean)
		matrix `sch'balance[`row'+2,`col'+2] = r(N)
	summ `dep'_bl if `sch'matchedsample==1 & `sch'_elig==1
		matrix `sch'balance[`row',`col'+3]   = r(mean)
		matrix `sch'balance[`row'+2,`col'+3] = r(N)
	summ `dep'_bl if `sch'matchedsample==1 & `sch'_elig==0
		matrix `sch'balance[`row',`col'+4]   = r(mean)
		matrix `sch'balance[`row'+2,`col'+4] = r(N)
		
	if "`dep'"=="lep" | "`dep'"=="c_state_mrawsc" | "`dep'"=="c_state_erawsc" {
	*No pre-baseline controls
		areg `dep'_bl `sch'_elig if `sch'matchedsample==1, r absorb(`sch'cell)
			matrix `sch'balance[`row',`col'+5]   = _b[`sch'_elig]
			matrix `sch'balance[`row'+1,`col'+5] = _se[`sch'_elig]
			matrix `sch'balance[`row'+2,`col'+5] = e(N)

	*Pre-baseline controls
		areg `dep'_bl `sch'_elig c_state_mrawsc_pbl c_state_erawsc_pbl if `sch'matchedsample==1, r absorb(`sch'cell)
			matrix `sch'balance[`row',`col'+6]   = _b[`sch'_elig]
			matrix `sch'balance[`row'+1,`col'+6] = _se[`sch'_elig]
			matrix `sch'balance[`row'+2,`col'+6] = e(N)
			
	*Pre-baseline controls - first exposure year
		areg `dep'_bl `sch'_elig c_state_mrawsc_pbl c_state_erawsc_pbl if `sch'matchedsample==1 & no_attrit0==1 & no_attrit1==1, r absorb(`sch'cell)
			matrix `sch'balance[`row',`col'+7]   = _b[`sch'_elig]
			matrix `sch'balance[`row'+1,`col'+7] = _se[`sch'_elig]
			matrix `sch'balance[`row'+2,`col'+7] = e(N)	
	}
	if "`dep'"=="gr_m" | "`dep'"=="gr_e" {
	*No pre-baseline controls
		areg `dep'_bl `sch'_elig if `sch'matchedsample==1, r absorb(`sch'cell)
			matrix `sch'balance[`row',`col'+5]   = _b[`sch'_elig]
			matrix `sch'balance[`row'+1,`col'+5] = _se[`sch'_elig]
			matrix `sch'balance[`row'+2,`col'+5] = e(N)

	*Pre-baseline controls
		areg `dep'_bl `sch'_elig c_state_mrawsc_ppbl c_state_erawsc_ppbl if `sch'matchedsample==1, r absorb(`sch'cell)
			matrix `sch'balance[`row',`col'+6]   = _b[`sch'_elig]
			matrix `sch'balance[`row'+1,`col'+6] = _se[`sch'_elig]
			matrix `sch'balance[`row'+2,`col'+6] = e(N)
			
	*Pre-baseline controls - first exposure year
		areg `dep'_bl `sch'_elig c_state_mrawsc_ppbl c_state_erawsc_ppbl if `sch'matchedsample==1 & no_attrit0==1 & no_attrit1==1, r absorb(`sch'cell)
			matrix `sch'balance[`row',`col'+7]   = _b[`sch'_elig]
			matrix `sch'balance[`row'+1,`col'+7] = _se[`sch'_elig]
			matrix `sch'balance[`row'+2,`col'+7] = e(N)	
	}	
			
	local row = `row' + 3

	} 
	foreach dep in no_attrit0 no_attrit1 no_attrit2 {

	summ `dep' if compsample==1
		matrix `sch'balance[`row',`col']   = r(mean)
		matrix `sch'balance[`row'+2,`col'] = r(N)
	summ `dep' if compsample==1 & compsamp_charter==1
		matrix `sch'balance[`row',`col'+1]   = r(mean)
		matrix `sch'balance[`row'+2,`col'+1] = r(N)
	summ `dep' if anymatchedsample==1 & any_attendany==1
		matrix `sch'balance[`row',`col'+2]   = r(mean)
		matrix `sch'balance[`row'+2,`col'+2] = r(N)
	summ `dep' if anymatchedsample==1 & any_elig==1
		matrix `sch'balance[`row',`col'+3]   = r(mean)
		matrix `sch'balance[`row'+2,`col'+3] = r(N)
	summ `dep' if anymatchedsample==1 & any_elig==0
		matrix `sch'balance[`row',`col'+4]   = r(mean)
		matrix `sch'balance[`row'+2,`col'+4] = r(N)

	areg `dep' any_elig if anymatchedsample==1, r absorb(anycell)
		matrix `sch'balance[`row',`col'+5]   = _b[any_elig]
		matrix `sch'balance[`row'+1,`col'+5] = _se[any_elig]
		matrix `sch'balance[`row'+2,`col'+5] = e(N)
	local row = `row' + 3

	} 
}
clear
svmat dearharbbalance
svmat OGbalance
br
}

***********************************************
* Grandfathering IV                           *
***********************************************
if `gfing' == 1 {

	use "${simsdata}simsoct09", clear
	keep if (grade=="06" | grade=="07") & (school=="00350257" | school=="00350426" | school=="00350074")
	gen inOG   = (school=="00350257")
	gen inharbdear = (school=="00350426" | school=="00350074")
	keep sasid grade inOG inharbdear
	destring grade, replace
	destring sasid, replace
	reshape wide inOG inharbdear, i(sasid) j(grade)
	gen year6 = 2010
	gen year7 = 2010
	tempfile legacy_fa09
	save `legacy_fa09', replace

	use if (boston5==1 | boston6==1) & year5>=2006 using "${cleandata}obsfile_2014", clear 

	merge 1:1 sasid using `legacy_fa09', keep(master match)
	drop _merge

	preserve
	egen harbdearcell = group(year5 masscode5 hisp5 black5 asian5 otherrace5 female5 sped5 foodfred5)
	forval g = 6/7 {
		gen harbdear_elig`g' = 1 if inharbdear`g'==1 & year`g'==2010
		replace harbdear_elig`g' = 0 if inharbdear`g'==. & year`g'==2010
		bys harbdearcell: egen temp = max(harbdear_elig`g') if harbdearcell!=.
		gen harbdearsamp`g' = (temp==1 & harbdear_elig`g'!=.) & boston5==1 & bcharter5==0 & spedalt5==0
		drop temp
	}

	replace harbdearsamp7 = 0 if harbdearsamp6 == 1
	gen harbdearmatchedsample = (harbdearsamp6==1 | harbdearsamp7==1)
	gen harbdear_elig = harbdear_elig6 if harbdearsamp6
	replace harbdear_elig = harbdear_elig7 if harbdearsamp7
	replace harbdearcell = . if harbdearmatchedsample == 0
	keep if harbdearmatchedsample

	drop *_bl
	foreach var in `covstems' white boston bcharter spedalt c_state_mrawsc c_state_erawsc year {
		gen `var'_bl   = `var'5 if harbdearmatchedsample == 1
		gen `var'_pbl  = `var'4 if harbdearmatchedsample == 1
		gen `var'_ppbl = `var'3 if harbdearmatchedsample == 1
	}

	tempfile harbdear 
	save `harbdear', replace

	restore
	drop *_bl
	foreach var in `covstems' white boston bcharter spedalt c_state_mrawsc c_state_erawsc year {
		gen `var'_bl = .
		gen `var'_pbl = .
		gen `var'_ppbl = .
	}
	forval bl = 5/6 {

		gen sch_lettergrade_adj`bl' = .
		bys year`bl' masscode`bl': egen sch_boston`bl' = mean(boston`bl')
		gen temp = c_boston_mrawsc`bl' + c_boston_erawsc`bl' if sch_boston`bl' == 1
		bys year`bl' masscode`bl': egen sch_rawsc`bl' = mean(temp) if sch_boston`bl'==1
		drop temp
		bys year`bl' masscode`bl': gen temp = _n==1
		levelsof year`bl', local(years)
		foreach y in `years' {
			cap {
				pctile cut = sch_rawsc`bl' if temp==1 & year`bl'==`y', n(10)
				xtile bin = sch_rawsc`bl', cutpoint(cut)
				replace sch_lettergrade_adj`bl' = bin if year`bl'==`y'
				drop cut bin
			}
		}
		drop temp

		egen OGcell`bl' = group(year`bl' sch_lettergrade_adj`bl' hisp`bl' black`bl' asian`bl' otherrace`bl' female`bl' sped`bl' foodfred`bl')

		local g = `bl' + 1
		gen OG_elig`g' = 1 if inOG`g'==1 & year`g'==2010
		replace OG_elig`g' = 0 if inOG`g'==. & year`g'==2010
		bys OGcell`bl': egen temp = max(OG_elig`g') if OGcell`bl'!=.
		gen OGsamp`g' = (temp==1 & OG_elig`g'!=.) & boston`bl'==1 & bcharter`bl'==0 & spedalt`bl'==0
		replace OGcell`bl' = . if OGsamp`g'==0

		foreach var in `covstems' white boston bcharter spedalt c_state_mrawsc c_state_erawsc year {
			replace `var'_bl = `var'`bl' if OGsamp`g'==1
			local pbl = `bl' -1
			replace `var'_pbl = `var'`pbl' if OGsamp`g'==1
			local ppbl = `bl' - 2
			replace `var'_ppbl = `var'`ppbl' if OGsamp`g'==1			
		}	

		drop sch_rawsc`bl' temp
	}

	egen OGcell = group(OGcell5 OGcell6), m
	gen OGmatchedsample = (OGsamp6==1 | OGsamp7==1)
	gen OG_elig = OG_elig6 if OGsamp6
	replace OG_elig = OG_elig7 if OGsamp7
	replace OGcell = . if OGmatchedsample == 0
	append using `harbdear'

	foreach v in _elig matchedsample cell samp6 samp7 {
		replace harbdear`v' = 0 if harbdear`v' == .
		replace OG`v' = 0 if OG`v' == .
	}

	gen anymatchedsample = OGmatchedsample | harbdearmatchedsample
	gen samp6 = OGsamp6 | harbdearsamp6
	gen samp7 = OGsamp7 | harbdearsamp7
	gen any_elig = harbdear_elig | OG_elig==1
	egen anycell = group(harbdearcell OGcell)

	keep if anymatchedsample == 1

	qui tab anycell, gen(d_any_)
	qui tab harbdearcell, gen(d_harbdear_)
	qui tab OGcell, gen(d_OG_)

	foreach g in 7 8 {
		gen attend_OG`g' = (masscode`g'  == 350257)
		replace attend_OG`g' = . if masscode`g' == .
		gen attend_harbdear`g' = (masscode`g'  == 350426 | masscode`g' == 350074)
		replace attend_harbdear`g' = . if masscode`g' == .
		gen attend_any`g' = (masscode`g'  == 350257 | masscode`g'  == 350426 | masscode`g' == 350074)
	}

	foreach sch in OG harbdear any {
		gen years`sch'1 = attend_`sch'7 if samp6==1
		replace years`sch'1 = attend_`sch'8 if samp7==1
		gen years`sch'2 = attend_`sch'8	+ attend_`sch'7*(1+repeat7) if samp6==1
	}

	foreach s in m e {
		forval out = 1/2 {
			gen `s'score`out' =.	
			forval g = 6/7{
				local `out'gr = `g' + `out'
				cap replace `s'score`out' = c_state_`s'rawsc``out'gr' if samp`g'==1
			}
		}
	}


	expand 2, gen(copy)
	qui tab copy, gen(copy_)
	
	foreach sch in OG harbdear any {
		gen years`sch' =.
		forval out = 1/2{
			local i = `out' -1
			replace years`sch' = years`sch'`out' if copy==`i'
		}
	}
	
	
	gen grade = 7 if (samp6==1 & copy==0)
	replace grade = 8 if (samp6==1 & copy==1) | (samp7==1 & copy==0)

	qui tab grade, gen(outg_)
	
	gen yearout=.
	forval g=7/8{
		replace yearout = year`g' if grade==`g'
	}
	
	qui tab yearout, gen(y_)

	tempfile usethis
	sa `usethis'
	
/************************ First-differenced IV ************************/		
	
	if `fdiv'==1 {

		foreach sch in harbdear OG {
			
			u `usethis', clear
			
			local row = 1
			local col = 1
			
			foreach s in m e {
			
					*Condition for being included in regression - have legacy year scores
					global cond "`sch'matchedsample==1 & ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1))"
					
					gen `s'score =.
					replace `s'score = `s'score1 if copy==0
					replace `s'score = `s'score2 if copy==1
					
					gen nongf`s'score = `s'score
						
					replace `s'score = `s'score - c_state_`s'rawsc6 if samp6==1
					replace `s'score = `s'score - c_state_`s'rawsc7 if samp7==1
						
					summ nongf`s'score if `sch'_elig==0 & ${cond}
						matrix `sch'gfing[`row',`col']   = r(mean)
						matrix `sch'gfing[`row'+2,`col'] = r(N)
						
					/* OLS */
					reg `s'score years`sch'  d_`sch'_* `controls' copy_* outg_* if ${cond}, cluster(sasid)
						matrix `sch'gfing[`row',`col'+1]   = _b[years`sch']
						matrix `sch'gfing[`row'+1,`col'+1] = _se[years`sch']
						matrix `sch'gfing[`row'+2,`col'+1] = e(N)	
					
					/* RF */
					areg `s'score `sch'_elig `controls' copy_* outg_* if years`sch'!=. & ${cond}, absorb(`sch'cell) cluster(sasid)
						matrix `sch'gfing[`row',`col'+2]   = _b[`sch'_elig]
						matrix `sch'gfing[`row'+1,`col'+2] = _se[`sch'_elig]
						matrix `sch'gfing[`row'+2,`col'+2] = e(N)
					
					/* FS */
					areg years`sch' `sch'_elig `controls' copy_* outg_* if `s'score!=. & ${cond}, absorb(`sch'cell) cluster(sasid)
						matrix `sch'gfing[`row',`col'+3]   = _b[`sch'_elig]
						matrix `sch'gfing[`row'+1,`col'+3] = _se[`sch'_elig]
						matrix `sch'gfing[`row'+2,`col'+3] = e(N)
						
					/* 2SLS */
					ivreg2 `s'score (years`sch' = `sch'_elig) d_`sch'_* `controls' copy_* outg_* if ${cond}, partial(d_`sch'_* copy_*) cluster(sasid)
						matrix `sch'gfing[`row',`col'+4]   = _b[years]
						matrix `sch'gfing[`row'+1,`col'+4] = _se[years]
						matrix `sch'gfing[`row'+2,`col'+4] = e(N)

				local row = `row' + 3
			}

			foreach copy in 0 1 {
				
				foreach s in m e {
					
					*Condition for being included in regression - have legacy year scores
					global cond "`sch'matchedsample==1 & ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1))"

					summ nongf`s'score if `sch'_elig==0 & ${cond} & copy==`copy'
						matrix `sch'gfing[`row',`col']   = r(mean)
						matrix `sch'gfing[`row'+2,`col'] = r(N)
						
						/* OLS */
						reg `s'score years`sch'  d_`sch'_* `controls' copy_* outg_* if ${cond} & copy==`copy', cluster(sasid)
							matrix `sch'gfing[`row',`col'+1]   = _b[years`sch']
							matrix `sch'gfing[`row'+1,`col'+1] = _se[years`sch']
							matrix `sch'gfing[`row'+2,`col'+1] = e(N)	
						
						/* RF */
						areg `s'score `sch'_elig `controls' copy_* outg_* if years`sch'!=. & ${cond} & copy==`copy', absorb(`sch'cell) cluster(sasid)
							matrix `sch'gfing[`row',`col'+2]   = _b[`sch'_elig]
							matrix `sch'gfing[`row'+1,`col'+2] = _se[`sch'_elig]
							matrix `sch'gfing[`row'+2,`col'+2] = e(N)
						
						/* FS */
						areg years`sch' `sch'_elig `controls' copy_* outg_* if `s'score!=. & ${cond} & copy==`copy', absorb(`sch'cell) cluster(sasid)
							matrix `sch'gfing[`row',`col'+3]   = _b[`sch'_elig]
							matrix `sch'gfing[`row'+1,`col'+3] = _se[`sch'_elig]
							matrix `sch'gfing[`row'+2,`col'+3] = e(N)
							
						/* 2SLS */
						ivreg2 `s'score (years`sch' = `sch'_elig) d_`sch'_* `controls' copy_* outg_* if ${cond} & copy==`copy', partial(d_`sch'_* copy_*) cluster(sasid)
							matrix `sch'gfing[`row',`col'+4]   = _b[years]
							matrix `sch'gfing[`row'+1,`col'+4] = _se[years]
							matrix `sch'gfing[`row'+2,`col'+4] = e(N)

				local row = `row' + 3
				}
			}

			foreach grade in 7 8 {
				
				foreach s in m e {
				
					*Condition for being included in regression - have legacy year scores
					global cond "`sch'matchedsample==1 & ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1))"
					
					summ nongf`s'score if `sch'_elig==0 & ${cond} & grade==`grade'
						matrix `sch'gfing[`row',`col']   = r(mean)
						matrix `sch'gfing[`row'+2,`col'] = r(N)
						
						/* OLS */
						reg `s'score years`sch'  d_`sch'_* `controls' copy_* outg_* if ${cond} & grade==`grade', cluster(sasid)
							matrix `sch'gfing[`row',`col'+1]   = _b[years`sch']
							matrix `sch'gfing[`row'+1,`col'+1] = _se[years`sch']
							matrix `sch'gfing[`row'+2,`col'+1] = e(N)	
						
						/* RF */
						areg `s'score `sch'_elig `controls' copy_* outg_* if years`sch'!=. & ${cond} & grade==`grade', absorb(`sch'cell) cluster(sasid)
							matrix `sch'gfing[`row',`col'+2]   = _b[`sch'_elig]
							matrix `sch'gfing[`row'+1,`col'+2] = _se[`sch'_elig]
							matrix `sch'gfing[`row'+2,`col'+2] = e(N)
						
						/* FS */
						areg years`sch' `sch'_elig `controls' copy_* outg_* if `s'score!=. & ${cond} & grade==`grade', absorb(`sch'cell) cluster(sasid)
							matrix `sch'gfing[`row',`col'+3]   = _b[`sch'_elig]
							matrix `sch'gfing[`row'+1,`col'+3] = _se[`sch'_elig]
							matrix `sch'gfing[`row'+2,`col'+3] = e(N)
							
						/* 2SLS */
						ivreg2 `s'score (years`sch' = `sch'_elig) d_`sch'_* `controls' copy_* outg_* if ${cond} & grade==`grade', partial(d_`sch'_* copy_*) cluster(sasid)
							matrix `sch'gfing[`row',`col'+4]   = _b[years]
							matrix `sch'gfing[`row'+1,`col'+4] = _se[years]
							matrix `sch'gfing[`row'+2,`col'+4] = e(N)

					local row = `row' + 3
				}
			}

	
		
		}
	}
	/********************************* Levels IV *********************************/	
	if `nofdiv'==1 {
		foreach sch in harbdear OG {
			
			u `usethis', clear
			
			local row = 1
			local col = 1
			
			foreach s in m e {
			
					*Condition for being included in regression - have legacy year scores
					global cond "`sch'matchedsample==1 & ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1))"
					
					gen `s'score =.
					replace `s'score = `s'score1 if copy==0
					replace `s'score = `s'score2 if copy==1
					
					gen nongf`s'score = `s'score
						
					summ nongf`s'score if `sch'_elig==0 & ${cond}
						matrix `sch'gfingnoFD[`row',`col']   = r(mean)
						matrix `sch'gfingnoFD[`row'+2,`col'] = r(N)
						
					/* OLS */
					reg `s'score years`sch'  d_`sch'_* `controls' copy_* outg_* if ${cond}, cluster(sasid)
						matrix `sch'gfingnoFD[`row',`col'+1]   = _b[years`sch']
						matrix `sch'gfingnoFD[`row'+1,`col'+1] = _se[years`sch']
						matrix `sch'gfingnoFD[`row'+2,`col'+1] = e(N)	
					
					/* RF */
					areg `s'score `sch'_elig `controls' copy_* outg_* if years`sch'!=. & ${cond}, absorb(`sch'cell) cluster(sasid)
						matrix `sch'gfingnoFD[`row',`col'+2]   = _b[`sch'_elig]
						matrix `sch'gfingnoFD[`row'+1,`col'+2] = _se[`sch'_elig]
						matrix `sch'gfingnoFD[`row'+2,`col'+2] = e(N)
					
					/* FS */
					areg years`sch' `sch'_elig `controls' copy_* outg_* if `s'score!=. & ${cond}, absorb(`sch'cell) cluster(sasid)
						matrix `sch'gfingnoFD[`row',`col'+3]   = _b[`sch'_elig]
						matrix `sch'gfingnoFD[`row'+1,`col'+3] = _se[`sch'_elig]
						matrix `sch'gfingnoFD[`row'+2,`col'+3] = e(N)
						
					/* 2SLS */
					ivreg2 `s'score (years`sch' = `sch'_elig) d_`sch'_* `controls' copy_* outg_* if ${cond}, partial(d_`sch'_* copy_*) cluster(sasid)
						matrix `sch'gfingnoFD[`row',`col'+4]   = _b[years]
						matrix `sch'gfingnoFD[`row'+1,`col'+4] = _se[years]
						matrix `sch'gfingnoFD[`row'+2,`col'+4] = e(N)

				local row = `row' + 3
			}

			foreach copy in 0 1 {
				
				foreach s in m e {
					
					*Condition for being included in regression - have legacy year scores
					global cond "`sch'matchedsample==1 & ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1))"

					summ nongf`s'score if `sch'_elig==0 & ${cond} & copy==`copy'
						matrix `sch'gfingnoFD[`row',`col']   = r(mean)
						matrix `sch'gfingnoFD[`row'+2,`col'] = r(N)
						
						/* OLS */
						reg `s'score years`sch'  d_`sch'_* `controls' copy_* outg_* if ${cond} & copy==`copy', cluster(sasid)
							matrix `sch'gfingnoFD[`row',`col'+1]   = _b[years`sch']
							matrix `sch'gfingnoFD[`row'+1,`col'+1] = _se[years`sch']
							matrix `sch'gfingnoFD[`row'+2,`col'+1] = e(N)	
						
						/* RF */
						areg `s'score `sch'_elig `controls' copy_* outg_* if years`sch'!=. & ${cond} & copy==`copy', absorb(`sch'cell) cluster(sasid)
							matrix `sch'gfingnoFD[`row',`col'+2]   = _b[`sch'_elig]
							matrix `sch'gfingnoFD[`row'+1,`col'+2] = _se[`sch'_elig]
							matrix `sch'gfingnoFD[`row'+2,`col'+2] = e(N)
						
						/* FS */
						areg years`sch' `sch'_elig `controls' copy_* outg_* if `s'score!=. & ${cond} & copy==`copy', absorb(`sch'cell) cluster(sasid)
							matrix `sch'gfingnoFD[`row',`col'+3]   = _b[`sch'_elig]
							matrix `sch'gfingnoFD[`row'+1,`col'+3] = _se[`sch'_elig]
							matrix `sch'gfingnoFD[`row'+2,`col'+3] = e(N)
							
						/* 2SLS */
						ivreg2 `s'score (years`sch' = `sch'_elig) d_`sch'_* `controls' copy_* outg_* if ${cond} & copy==`copy', partial(d_`sch'_* copy_*) cluster(sasid)
							matrix `sch'gfingnoFD[`row',`col'+4]   = _b[years]
							matrix `sch'gfingnoFD[`row'+1,`col'+4] = _se[years]
							matrix `sch'gfingnoFD[`row'+2,`col'+4] = e(N)

				local row = `row' + 3
				}
			}

			foreach grade in 7 8 {
				
				foreach s in m e {
				
					*Condition for being included in regression - have legacy year scores
					global cond "`sch'matchedsample==1 & ((c_state_`s'rawsc6!=. & samp6==1) | (c_state_`s'rawsc7!=. & samp7==1))"
					
					summ nongf`s'score if `sch'_elig==0 & ${cond} & grade==`grade'
						matrix `sch'gfingnoFD[`row',`col']   = r(mean)
						matrix `sch'gfingnoFD[`row'+2,`col'] = r(N)
						
						/* OLS */
						reg `s'score years`sch'  d_`sch'_* `controls' copy_* outg_* if ${cond} & grade==`grade', cluster(sasid)
							matrix `sch'gfingnoFD[`row',`col'+1]   = _b[years`sch']
							matrix `sch'gfingnoFD[`row'+1,`col'+1] = _se[years`sch']
							matrix `sch'gfingnoFD[`row'+2,`col'+1] = e(N)	
						
						/* RF */
						areg `s'score `sch'_elig `controls' copy_* outg_* if years`sch'!=. & ${cond} & grade==`grade', absorb(`sch'cell) cluster(sasid)
							matrix `sch'gfingnoFD[`row',`col'+2]   = _b[`sch'_elig]
							matrix `sch'gfingnoFD[`row'+1,`col'+2] = _se[`sch'_elig]
							matrix `sch'gfingnoFD[`row'+2,`col'+2] = e(N)
						
						/* FS */
						areg years`sch' `sch'_elig `controls' copy_* outg_* if `s'score!=. & ${cond} & grade==`grade', absorb(`sch'cell) cluster(sasid)
							matrix `sch'gfingnoFD[`row',`col'+3]   = _b[`sch'_elig]
							matrix `sch'gfingnoFD[`row'+1,`col'+3] = _se[`sch'_elig]
							matrix `sch'gfingnoFD[`row'+2,`col'+3] = e(N)
							
						/* 2SLS */
						ivreg2 `s'score (years`sch' = `sch'_elig) d_`sch'_* `controls' copy_* outg_* if ${cond} & grade==`grade', partial(d_`sch'_* copy_*) cluster(sasid)
							matrix `sch'gfingnoFD[`row',`col'+4]   = _b[years]
							matrix `sch'gfingnoFD[`row'+1,`col'+4] = _se[years]
							matrix `sch'gfingnoFD[`row'+2,`col'+4] = e(N)

					local row = `row' + 3
				}
			}

	
		
		}
		
	}
	clear
	svmat harbdeargfing 
	svmat OGgfing 
	svmat harbdeargfingnoFD 
	svmat OGgfingnoFD
}
***********************************************
* Peer composition                            *
***********************************************
if `switching' == 1 {

	use "${simsdata}simsoct09", clear
	keep if (grade=="06" | grade=="07") & (school=="00350257" | school=="00350426" | school=="00350074")
	gen inOG   = (school=="00350257")
	gen inharbdear = (school=="00350426" | school=="00350074")
	keep sasid grade inOG inharbdear
	destring grade, replace
	destring sasid, replace
	reshape wide inOG inharbdear, i(sasid) j(grade)
	gen year6 = 2010
	gen year7 = 2010
	tempfile legacy_fa09
	save `legacy_fa09', replace

	use if (boston5==1 | boston6==1) & year5>=2006 using "${cleandata}obsfile_2014_wpeers", clear 

	merge 1:1 sasid using `legacy_fa09', keep(master match)
	drop _merge

	preserve
	egen harbdearcell = group(year5 masscode5 hisp5 black5 asian5 otherrace5 female5 sped5 foodfred5)
	forval g = 6/7 {
		gen harbdear_elig`g' = 1 if inharbdear`g'==1 & year`g'==2010
		replace harbdear_elig`g' = 0 if inharbdear`g'==. & year`g'==2010
		bys harbdearcell: egen temp = max(harbdear_elig`g') if harbdearcell!=.
		gen harbdearsamp`g' = (temp==1 & harbdear_elig`g'!=.) & boston5==1 & bcharter5==0 & spedalt5==0
		drop temp
	}

	replace harbdearsamp7 = 0 if harbdearsamp6 == 1
	gen harbdearmatchedsample = (harbdearsamp6==1 | harbdearsamp7==1)
	gen harbdear_elig = harbdear_elig6 if harbdearsamp6
	replace harbdear_elig = harbdear_elig7 if harbdearsamp7
	replace harbdearcell = . if harbdearmatchedsample == 0
	keep if harbdearmatchedsample

	drop *_bl
	foreach var in `covstems' white boston bcharter spedalt c_state_mrawsc c_state_erawsc year {
		gen `var'_bl   = `var'5 if harbdearmatchedsample == 1
		gen `var'_pbl  = `var'4 if harbdearmatchedsample == 1
		gen `var'_ppbl = `var'3 if harbdearmatchedsample == 1
	}

	tempfile harbdear 
	save `harbdear', replace

	restore
	drop *_bl
	foreach var in `covstems' white boston bcharter spedalt c_state_mrawsc c_state_erawsc year {
		gen `var'_bl = .
		gen `var'_pbl = .
		gen `var'_ppbl = .
	}
	forval bl = 5/6 {

		gen sch_lettergrade_adj`bl' = .
		bys year`bl' masscode`bl': egen sch_boston`bl' = mean(boston`bl')
		gen temp = c_boston_mrawsc`bl' + c_boston_erawsc`bl' if sch_boston`bl' == 1
		bys year`bl' masscode`bl': egen sch_rawsc`bl' = mean(temp) if sch_boston`bl'==1
		drop temp
		bys year`bl' masscode`bl': gen temp = _n==1
		levelsof year`bl', local(years)
		foreach y in `years' {
			cap {
				pctile cut = sch_rawsc`bl' if temp==1 & year`bl'==`y', n(10)
				xtile bin = sch_rawsc`bl', cutpoint(cut)
				replace sch_lettergrade_adj`bl' = bin if year`bl'==`y'
				drop cut bin
			}
		}
		drop temp

		egen OGcell`bl' = group(year`bl' sch_lettergrade_adj`bl' hisp`bl' black`bl' asian`bl' otherrace`bl' female`bl' sped`bl' foodfred`bl')

		local g = `bl' + 1
		gen OG_elig`g' = 1 if inOG`g'==1 & year`g'==2010
		replace OG_elig`g' = 0 if inOG`g'==. & year`g'==2010
		bys OGcell`bl': egen temp = max(OG_elig`g') if OGcell`bl'!=.
		gen OGsamp`g' = (temp==1 & OG_elig`g'!=.) & boston`bl'==1 & bcharter`bl'==0 & spedalt`bl'==0
		replace OGcell`bl' = . if OGsamp`g'==0

		foreach var in `covstems' white boston bcharter spedalt c_state_mrawsc c_state_erawsc year {
			replace `var'_bl = `var'`bl' if OGsamp`g'==1
			local pbl = `bl' -1
			replace `var'_pbl = `var'`pbl' if OGsamp`g'==1
			local ppbl = `bl' - 2
			replace `var'_ppbl = `var'`ppbl' if OGsamp`g'==1			
		}	

		drop sch_rawsc`bl' temp
	}

	egen OGcell = group(OGcell5 OGcell6), m
	gen OGmatchedsample = (OGsamp6==1 | OGsamp7==1)
	gen OG_elig = OG_elig6 if OGsamp6
	replace OG_elig = OG_elig7 if OGsamp7
	replace OGcell = . if OGmatchedsample == 0
	append using `harbdear'

	foreach v in _elig matchedsample cell samp6 samp7 {
		replace harbdear`v' = 0 if harbdear`v' == .
		replace OG`v' = 0 if OG`v' == .
	}

	gen anymatchedsample = OGmatchedsample | harbdearmatchedsample
	gen samp6 = OGsamp6 | harbdearsamp6
	gen samp7 = OGsamp7 | harbdearsamp7
	gen any_elig = harbdear_elig | OG_elig==1
	egen anycell = group(harbdearcell OGcell)

	keep if anymatchedsample == 1

	qui tab anycell, gen(d_any_)
	qui tab harbdearcell, gen(d_harbdear_)
	qui tab OGcell, gen(d_OG_)

	foreach g in 7 8 {
		gen attend_OG`g' = (masscode`g'  == 350257)
		replace attend_OG`g' = . if masscode`g' == .
		gen attend_harbdear`g' = (masscode`g'  == 350426 | masscode`g' == 350074)
		replace attend_harbdear`g' = . if masscode`g' == .
		gen attend_any`g' = (masscode`g'  == 350257 | masscode`g'  == 350426 | masscode`g' == 350074)
	}

	foreach sch in OG harbdear any {
		gen years`sch'1 = attend_`sch'7 if samp6==1
		replace years`sch'1 = attend_`sch'8 if samp7==1
		gen years`sch'2 = attend_`sch'8	+ attend_`sch'7*(1+repeat7) if samp6==1
	}

	foreach s in m e {
		forval out = 1/2 {
			gen `s'score`out' =.	
			forval g = 6/7{
				local `out'gr = `g' + `out'
				cap replace `s'score`out' = c_state_`s'rawsc``out'gr' if samp`g'==1
			}
		}
	}


	
	gen grade = 7 if samp6==1
	replace grade = 8 if samp7==1
	qui tab grade, gen(outg_)
	
	gen yearout=.
	forval g=7/8{
		replace yearout = year`g' if grade==`g'
	}
	qui tab yearout, gen(y_)

	foreach t in mscore escore sped lep foodfred {
		gen peer_`t'1 = peer_`t'_bl5_for7 if samp6==1
		replace peer_`t'1 = peer_`t'_bl5_for8 if harbdearsamp7==1
		replace peer_`t'1 = peer_`t'_bl6_for8 if OGsamp7==1
	}
	local row = 1
	local col = 1	
	foreach sch in harbdear OG {
		foreach t in sped foodfred lep mscore escore {
			summ peer_`t'1 if `sch'_elig==0 & `sch'matchedsample==1
				matrix switching[`row',`col']   = r(mean)
				matrix switching[`row'+2,`col'] = r(N)
				
			ivreg2 peer_`t'1 (years`sch'1 = `sch'_elig) d_`sch'_* outg_* `controls' if `sch'matchedsample==1, r partial(d_`sch'_* `controls')
				matrix switching[`row',`col'+1]   = _b[years`sch'1]
				matrix switching[`row'+1,`col'+1] = _se[years`sch'1]
				matrix switching[`row'+2,`col'+1] = e(N)
			local row = `row' + 3
		}
		local col = `col' + 3
		local row = 1
	}
	clear
	svmat switching
	br
}
