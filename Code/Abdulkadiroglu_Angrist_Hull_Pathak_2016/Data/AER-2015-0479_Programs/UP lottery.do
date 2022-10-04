* P.HULL
* UP Lottery Analysis

clear all
cap log close
set more off
set trace off
set matsize 10000
set maxvar 10000
set seed 42
adopath + "Z:\ado\plus\"

global sims        "Z:/charters/MA DESE data/SIMS/"
global lotterydata "Z:/charters/MA Lotteries/all lotto data/"
global combdata    "Z:/charters/Combined data/"
global cleandata   "X:/indistrict/Data/"
global bpsdata     "Z:/BPS/Peter/BPSVA/Data/"

global programs "Z:/charters/indistrict/Programs/"

local simsdate "10 30 2014"

local covstems "female black hisp asian otherrace sped lep foodfred fem_min"
local controls "hisp_bl black_bl white_bl asian_bl female_bl foodfred_bl sped_bl lep_bl c_state_mrawsc_bl c_state_erawsc_bl y_*"			

matrix balance       = J(200,200,.)
matrix lottery       = J(200,200,.)
matrix choice        = J(20,20,.)
matrix sample        = J(20,20,.)

local setup          = 1 /* Table B9             */
local switchingsetup = 1 

local balance        = 1 /* Table B8             */
local lottery        = 1 /* Table 8              */
local choice         = 1 /* Table B10 (Panel B)  */
local switching      = 1 /* Table 10 (Column 5)  */

***********************************************
* Construct analysis sample                   *
***********************************************
if `setup'==1 {

use "${lotterydata}up academy/matched/UP2011_Matched", clear
replace offeredadmissionever = . if lotterynumber>=258 /* end of september cut for ever offer */
qui count
matrix sample[1,1] = r(N)
drop if lateapplicant == 1
qui count
matrix sample[2,1] = r(N)
drop if prioritygroup == 3
qui count
matrix sample[3,1] = r(N)
drop if prioritygroup == 1
qui count
matrix sample[4,1] = r(N)
keep applicationyear sasid gradeapplying offeredadmissioninitially offeredadmissionever 
rename offeredadmissioninitially initial_offerUP
rename offeredadmissionever offerUP
gen applyUP = 1
tempfile up2011
save `up2011', replace

insheet using "${lotterydata}up academy/formatted/Spring 2011/Up Academy 2011 Gr 7 Formatted_PHEDIT.csv", clear
replace offeredadmissionever = . if offeredadmissionever == 0
qui count
matrix sample[4,2] = r(N)
keep applicationyear sasid gradeapplying offeredadmissioninitially offeredadmissionever 
rename offeredadmissioninitially initial_offerUP
rename offeredadmissionever offerUP
gen applyUP = 1
tempfile up2011gr7
save `up2011gr7', replace

use "${lotterydata}up academy/matched/UP2012_Matched", clear
qui count
matrix sample[1,3] = r(N)
drop if lateapplicant == 1
qui count
matrix sample[2,3] = r(N)
drop if prioritygroup == 3
qui count
matrix sample[3,3] = r(N)
drop if prioritygroup == 1
qui count
matrix sample[4,3] = r(N)
keep applicationyear sasid gradeapplying offeredadmissioninitially offeredadmissionever 
rename offeredadmissioninitially initial_offerUP
rename offeredadmissionever offerUP
gen applyUP = 1
tempfile up2012
save `up2012', replace

use `up2011'
append using `up2011gr7'
append using `up2012' 
drop if sasid==.

rename applicationyear yearapp
rename gradeapplying gradeapp

replace applyUP = 0 if applyUP==.

foreach s in UP {
	egen riskset_`s' = group(yearapp gradeapp apply`s')
	replace initial_offer`s' = 0 if initial_offer`s'==. 
	replace offer`s' = 0 if offer`s'==. 
	replace offer`s' = 1 if initial_offer`s' == 1
}
egen riskset = group(riskset_*)
qui tab riskset, gen(d_)

merge m:1 sasid using "${cleandata}simsmcaswide `simsdate'"
forval g = 3/8 {
		foreach t in m e {
			qui replace c_state_`t'rawsc`g' = c_boston_`t'rawsc`g'
		}
}
gen temp1 = year5==2011 | year5==2012
gen temp2 = year6==2011
keep if temp1 | temp2
gen temp3 = runiform()>0.5 if (temp1 & temp2)
replace temp2 = 0 if temp3==1
replace temp1 = 0 if temp3==0
gen baselinegrade = gradeapp-1
replace baselinegrade = 5 if temp1==1
replace baselinegrade = 6 if temp2==1
levelsof baselinegrade, local(b_grade)
drop temp*
foreach v in boston spedalt charter year `covstems' c_state_mrawsc c_state_erawsc {
	gen `v'_bl = .
	foreach g in `b_grade' {
		replace `v'_bl = `v'`g' if baselinegrade==`g'
	}
}
keep if (year_bl == yearapp) | yearapp == .

qui count if gradeapp==6 & yearapp==2011
matrix sample[5,1] = r(N)
qui count if gradeapp==7 & yearapp==2011
matrix sample[5,2] = r(N)
qui count if gradeapp==6 & yearapp==2012
matrix sample[5,3] = r(N)
qui count if gradeapp==6 & yearapp==2011 & boston_bl==1 & spedalt_bl==0 & charter_bl==0
matrix sample[6,1] = r(N)
qui count if gradeapp==7 & yearapp==2011 & boston_bl==1 & spedalt_bl==0 & charter_bl==0
matrix sample[6,2] = r(N)
qui count if gradeapp==6 & yearapp==2012 & boston_bl==1 & spedalt_bl==0 & charter_bl==0
matrix sample[6,3] = r(N)

forval g = 3/8 {
	foreach t in m e {
		qui replace c_state_`t'rawsc`g' = c_boston_`t'rawsc`g'
	}
}

gen white_bl = 1-(hisp_bl + black_bl + asian_bl + otherrace_bl)

gen hasbaselinedemos = !missing(female_bl, black_bl, hisp_bl, asian_bl, otherrace_bl, sped_bl, lep_bl, foodfred_bl)

foreach g in 6 7 8 9 10 11 12 {
	gen attend_UP`g' = masscode`g'  == 4800405
	replace attend_UP`g' = . if masscode`g' == .
}

gen has_1yr_out = .
replace has_1yr_out = (c_state_mrawsc6!=. & c_state_erawsc6!=.) if baselinegrade==5 & year5+1<=2014
replace has_1yr_out = (c_state_mrawsc7!=. & c_state_erawsc7!=.) if baselinegrade==6 & year6+1<=2014
gen has_2yr_out = .
replace has_2yr_out = (c_state_mrawsc7!=. & c_state_erawsc7!=.) if baselinegrade==5 & year5+2<=2014
replace has_2yr_out = (c_state_mrawsc8!=. & c_state_erawsc8!=.) if baselinegrade==6 & year6+2<=2014
gen has_3yr_out = (c_state_mrawsc8!=. & c_state_erawsc8!=.) if baselinegrade==5 & year5+3<=2014

compress
save "${cleandata}analysisfile_2014", replace		
clear
svmat sample
br
}
if `switchingsetup'==1 {
	use "${cleandata}simsmcaswide `simsdate'", clear
	forval g = 3/8 {
			foreach t in m e {
				qui replace c_state_`t'rawsc`g' = c_boston_`t'rawsc`g'
			}
	}
	keep sasid year6 masscode6 year7 masscode7 year8 masscode8 c_state_mrawsc5 c_state_erawsc5 sped5 lep5 foodfred5 c_state_mrawsc6 c_state_erawsc6 sped6 lep6 foodfred6 
	forval outgr = 6/8 {
		preserve
		collapse (mean) c_state_mrawsc5 c_state_erawsc5 sped5 lep5 foodfred5 c_state_mrawsc6 c_state_erawsc6 sped6 lep6 foodfred6, by(year`outgr' masscode`outgr')
		forval blgr = 5/6 {
			gen peer_mscore_bl`blgr'_for`outgr' = c_state_mrawsc`blgr'
			gen peer_escore_bl`blgr'_for`outgr' = c_state_erawsc`blgr'
			gen peer_sped_bl`blgr'_for`outgr' = sped`blgr'
			gen peer_lep_bl`blgr'_for`outgr' = lep`blgr'
			gen peer_foodfred_bl`blgr'_for`outgr' = foodfred`blgr'
		}
		drop if year`outgr' == . | masscode`outgr' == .
		keep masscode`outgr' year`outgr' peer_*
		tempfile peer_for`outgr'
		save `peer_for`outgr'', replace
		restore
	}	
	
	use "${cleandata}analysisfile_2014", clear
		drop _merge
	forval outgr = 6/8 {
		merge m:1 year`outgr' masscode`outgr' using `peer_for`outgr''
		drop _merge
	}	
	
	save "${cleandata}analysisfile_2014_wpeers", replace
}

***********************************************
* Balance regressions                         *
***********************************************
if `balance' == 1 {

local row = 1
local col = 1

use if boston_bl==1 & spedalt_bl==0 using "${cleandata}analysisfile_2014", clear 
replace applyUP = 0 if applyUP == 1 & gradeapp == 6 & charter5 == 1
replace applyUP = 0 if applyUP == 1 & gradeapp == 7 & charter6 == 1
gen compsample = 0
qui levelsof year_bl, local(years)
qui count if applyUP==1
local denom = r(N)
local prps ""
foreach y in `years' {
	qui count if year_bl==`y' & applyUP==1 & gradeapp==6
	local prop`y' = r(N)/`denom'
	local prps "`prps', `prop`y''"
}
qui count if year_bl==2011 & applyUP==1 & gradeapp==7
local prop62011 = r(N)/`denom'
local prps = substr("`prps'",2,.)
disp "`prps'"
local maxprop = max(`prps')
disp `maxprop'

foreach y in `years' {
	replace compsample = runiform()<(`prop`y''/`maxprop') if year_bl == `y' & baselinegr == 5 & compsample == 0 & spedalt5==0
}
replace compsample = runiform()<(`prop62011'/`maxprop') if year_bl == 2011 & baselinegr == 6 & compsample == 0 & spedalt6 == 0

replace offerUP = offerUP - initial_offerUP

local samp1  "compsample == 1"
local samp2  "applyUP == 1 & gradeapp==6"
local samp3  "applyUP == 1 & gradeapp==7"
local samp4  "applyUP == 1"
local samp5  "applyUP == 1"

forval s = 1/5 {
	local row = 1
	foreach dep in hisp black white asian female sped foodfred lep {
		summ `dep'_bl if `samp`s''
			matrix balance[`row',`col']   = r(mean)
			matrix balance[`row'+2,`col'] = r(N)
		local row = `row' +2
	}	

	local row = `row'+2
	foreach dep in c_state_mrawsc c_state_erawsc {
		summ `dep'_bl if `samp`s''
			matrix balance[`row',`col']   = r(mean)
			matrix balance[`row'+2,`col'] = r(N)
		local row = `row'+3
	}
	
	local row = `row'+1
	foreach dep in has_1yr_out has_2yr_out has_3yr_out  {
		summ `dep' if `samp`s''
			matrix balance[`row',`col']   = r(mean)
			matrix balance[`row'+2,`col'] = r(N)
		local row = `row'+3
	}

	local ++col

	if "`s'"=="5" {
		foreach z in initial_offer offer {
			local row = 1
			foreach dep in hisp black white asian female sped foodfred lep {
				reg `dep'_bl `z'UP d_* if `samp`s'' , r
					matrix balance[`row',`col']   = _b[`z'UP]
					matrix balance[`row'+1,`col'] = _se[`z'UP]				
					matrix balance[`row'+2,`col'] = e(N)
				reg `dep'_bl `z'UP d_* if `samp`s'' & has_1yr_out==1, r
					matrix balance[`row',`col'+2]   = _b[`z'UP]
					matrix balance[`row'+1,`col'+2] = _se[`z'UP]				
					matrix balance[`row'+2,`col'+2] = e(N)
				local row = `row' + 2
			}

			local row = `row'+2
			foreach dep in c_state_mrawsc c_state_erawsc {
				reg `dep'_bl `z'UP d_* if `samp`s'' , r
					matrix balance[`row',`col']   = _b[`z'UP]
					matrix balance[`row'+1,`col'] = _se[`z'UP]				
					matrix balance[`row'+2,`col'] = e(N)
				reg `dep'_bl `z'UP d_* if `samp`s'' & has_1yr_out==1 , r
					matrix balance[`row',`col'+2]   = _b[`z'UP]
					matrix balance[`row'+1,`col'+2] = _se[`z'UP]				
					matrix balance[`row'+2,`col'+2] = e(N)
				local row = `row'+3
			}

			local row = `row'+1
			foreach dep in has_1yr_out has_2yr_out has_3yr_out {
				cap {
				reg `dep' `z'UP d_* if `samp`s'' , r
					matrix balance[`row',`col']   = _b[`z'UP]
					matrix balance[`row'+1,`col'] = _se[`z'UP]				
					matrix balance[`row'+2,`col'] = e(N)
				}
				local row = `row'+3
			}
			local ++col
		}
		local ++col
	}
	
}

replace offerUP = offerUP + initial_offerUP
gen years_UP1 = attend_UP6 if gradeapp==6
replace years_UP1 = attend_UP7 if gradeapp==7

local row = 1
foreach dep in hisp_bl black_bl white_bl asian_bl female_bl sped_bl foodfred_bl lep_bl c_state_mrawsc_bl c_state_erawsc_bl has_1yr_out has_2yr_out has_3yr_out {
	
	qui reg offerUP d_* if applyUP == 1 & `dep'!=. & years_UP1!=.
	predict pZ, xb
	replace pZ = 0.000001 if pZ<=0
	replace pZ = 0.999999 if pZ>=1
	gen kappa = 1 - offerUP*(1-years_UP1)/pZ - (1-offerUP)*years_UP1/(1-pZ)
	summ `dep' if applyUP == 1 & `dep'!=. & years_UP1!=. [iw=kappa]
		matrix balance[`row',4]   = r(mean)
	summ kappa if applyUP == 1 & `dep'!=. & years_UP1!=.
		matrix balance[`row'+2,4] = `=`=r(N)'*`=r(mean)''
	drop pZ kappa
	if "`dep'" == "lep_bl" | "`dep'" == "c_state_erawsc_bl" {
		local row = `row' + 4	
	}
	else if "`dep'" == "c_state_mrawsc_bl" | "`dep'" == "has_1yr_out" | "`dep'" == "has_2yr_out"{
		local row = `row' + 3
	}
	else {
		local row = `row' + 2
	}
}

clear
svmat balance
br
}

***********************************************
* Lottery 2SLS                                *
***********************************************
if `lottery' == 1 {

use if boston_bl==1 & spedalt_bl==0 & charter_bl==0 using "${cleandata}analysisfile_2014", clear 
replace applyUP = 0 if applyUP == 1 & gradeapp == 6 & charter5 == 1
replace applyUP = 0 if applyUP == 1 & gradeapp == 7 & charter6 == 1
keep if applyUP == 1
gen years_UP1 = attend_UP6 if gradeapp==6
replace years_UP1 = attend_UP7 if gradeapp==7
gen years_UP2 = attend_UP7 + attend_UP6*(1+repeat6) if gradeapp==6
replace years_UP2 = attend_UP8 + attend_UP7*(1+repeat7) if gradeapp==7
gen years_UP3 = attend_UP8 + attend_UP7*(1+repeat7) + attend_UP6*(1+repeat6) if gradeapp==6

/* make offer waitlist offer (mutually exclusive instruments) */
replace offerUP = offerUP - initial_offerUP

expand 3
bys sasid: gen copy = _n
qui tab copy, gen(copy_)
gen years_UP = years_UP1 if copy==1
replace years_UP = years_UP2 if copy==2
replace years_UP = years_UP3 if copy==3

gen gradeout = 6 if copy==1 & gradeapp==6
replace gradeout = 7 if copy==2 & gradeapp==6
replace gradeout = 8 if copy==3 & gradeapp==6
replace gradeout = 7 if copy==1 & gradeapp==7
replace gradeout = 8 if copy==2 & gradeapp==7

qui tab gradeout, gen(og_)

gen yearout = year6 if copy==1 & gradeapp==6
replace yearout = year7 if copy==2 & gradeapp==6
replace yearout = year8 if copy==3 & gradeapp==6
replace yearout = year7 if copy==1 & gradeapp==7
replace yearout = year8 if copy==2 & gradeapp==7

qui tab yearout, gen(y_)

replace gradeout = 78 if gradeout==7 | gradeout==8
gen copyout = copy
replace copyout = 23 if copy==2 | copy==3
local row = 1
local col = 1

foreach s in m e {
		
		gen `s'level = c_state_`s'rawsc6 if copy==1 & gradeapp==6
		replace `s'level = c_state_`s'rawsc7 if copy==2 & gradeapp==6
		replace `s'level = c_state_`s'rawsc8 if copy==3 & gradeapp==6
		replace `s'level = c_state_`s'rawsc7 if copy==1 & gradeapp==7
		replace `s'level = c_state_`s'rawsc8 if copy==2 & gradeapp==7
		
		gen `s'score = c_state_`s'rawsc6 - c_state_`s'rawsc_bl if copy==1 & gradeapp==6
		replace `s'score = c_state_`s'rawsc7 - c_state_`s'rawsc_bl if copy==2 & gradeapp==6
		replace `s'score = c_state_`s'rawsc8 - c_state_`s'rawsc_bl if copy==3 & gradeapp==6
		replace `s'score = c_state_`s'rawsc7 - c_state_`s'rawsc_bl if copy==1 & gradeapp==7
		replace `s'score = c_state_`s'rawsc8 - c_state_`s'rawsc_bl if copy==2 & gradeapp==7
		
		*Summarize levels instead of gains
		qui summ `s'level if offerUP==0 & initial_offerUP==0 & `s'score!=.
			matrix lottery[`row',`col'] = r(mean)
			matrix lottery[`row'+1,`col'] = r(sd)
			matrix lottery[`row'+2,`col'] = r(N)
		
		/*OLS*/
		local col = `col'+2
		reg `s'score years_UP `controls' d_* copy_* og_*, cluster(sasid)
		matrix lottery[`row',`col'] = _b[years_UP]
			matrix lottery[`row'+1,`col'] = _se[years_UP]
			matrix lottery[`row'+2,`col'] = e(N)
		
		/* RF */
		local col = `col'+2
		reg `s'score initial_offerUP offerUP `controls' d_* copy_* og_* if years_UP!=., cluster(sasid)
			matrix lottery[`row',`col'] = _b[initial_offerUP]
			matrix lottery[`row'+1,`col'] = _se[initial_offerUP]
			matrix lottery[`row',`col'+1] = _b[offerUP]
			matrix lottery[`row'+1,`col'+1] = _se[offerUP]		
				
		/* FS */
		local col = `col'+2
		reg years_UP initial_offerUP offerUP `controls' d_* copy_* og_* if `s'score!=., cluster(sasid)
			matrix lottery[`row',`col'] = _b[initial_offerUP]
			matrix lottery[`row'+1,`col'] = _se[initial_offerUP]
			matrix lottery[`row',`col'+1] = _b[offerUP]
			matrix lottery[`row'+1,`col'+1] = _se[offerUP]		

		/* 2SLS */
		local col = `col'+2
		ivreg2 `s'score (years_UP = initial_offerUP offerUP) `controls' d_* copy_* og_*, cluster(sasid) partial(d_* `controls' copy_*)
		matrix lottery[`row',`col'] = _b[years_UP]
			matrix lottery[`row'+1,`col'] = _se[years_UP]
			matrix lottery[`row'+2,`col'] = e(N)
			matrix lottery[`row',`col'+1] = e(jp)

		local row = `row' + 3
		local col = 1
}

foreach copy in 1 23 {
	foreach s in m e {
	
		*Summarize levels instead of gains
		qui summ `s'level if offerUP==0 & initial_offerUP==0 & copyout==`copy'  & `s'score!=.
			matrix lottery[`row',`col'] = r(mean)
			matrix lottery[`row'+1,`col'] = r(sd)
			matrix lottery[`row'+2,`col'] = r(N)
				
		/* OLS */
		local col = `col'+2
		reg `s'score years_UP `controls' d_* copy_* og_* if copyout==`copy', cluster(sasid)
			matrix lottery[`row',`col'] = _b[years_UP]
			matrix lottery[`row'+1,`col'] = _se[years_UP]
			matrix lottery[`row'+2,`col'] = e(N)
		
		/* RF */
		local col = `col'+2
		reg `s'score initial_offerUP offerUP `controls' d_* copy_* og_* if years_UP!=. & copyout==`copy', cluster(sasid)
			matrix lottery[`row',`col'] = _b[initial_offerUP]
			matrix lottery[`row'+1,`col'] = _se[initial_offerUP]
			matrix lottery[`row',`col'+1] = _b[offerUP]
			matrix lottery[`row'+1,`col'+1] = _se[offerUP]		
				
		/* FS */
		local col = `col'+2
		reg years_UP initial_offerUP offerUP `controls' d_* copy_* og_* if `s'score!=. & copyout==`copy', cluster(sasid)
			matrix lottery[`row',`col'] = _b[initial_offerUP]
			matrix lottery[`row'+1,`col'] = _se[initial_offerUP]
			matrix lottery[`row',`col'+1] = _b[offerUP]
			matrix lottery[`row'+1,`col'+1] = _se[offerUP]		

		/* 2SLS */
		local col = `col'+2
		ivreg2 `s'score (years_UP = initial_offerUP offerUP) `controls' d_* copy_* og_* if copyout==`copy', cluster(sasid) partial(d_* `controls' copy_*)
			matrix lottery[`row',`col'] = _b[years_UP]
			matrix lottery[`row'+1,`col'] = _se[years_UP]
			matrix lottery[`row'+2,`col'] = e(N)
			matrix lottery[`row',`col'+1] = e(jp)

		local row = `row' + 3
		local col = 1
	}
}

foreach gr in 6 78 {
	foreach s in m e {
		
		*Summarize levels instead of gains
		qui summ `s'level if offerUP==0 & initial_offerUP==0 & gradeout==`gr'  & `s'score!=.
			matrix lottery[`row',`col'] = r(mean)
			matrix lottery[`row'+1,`col'] = r(sd)
			matrix lottery[`row'+2,`col'] = r(N)
				
		/* OLS */
		local col = `col'+2
		reg `s'score years_UP `controls' d_* copy_* og_* if gradeout==`gr', cluster(sasid)
			matrix lottery[`row',`col'] = _b[years_UP]
			matrix lottery[`row'+1,`col'] = _se[years_UP]
			matrix lottery[`row'+2,`col'] = e(N)
		
		/* RF */
		local col = `col'+2
		reg `s'score initial_offerUP offerUP `controls' d_* copy_* og_* if years_UP!=. & gradeout==`gr', cluster(sasid)
			matrix lottery[`row',`col'] = _b[initial_offerUP]
			matrix lottery[`row'+1,`col'] = _se[initial_offerUP]
			matrix lottery[`row',`col'+1] = _b[offerUP]
			matrix lottery[`row'+1,`col'+1] = _se[offerUP]		
				
		/* FS */
		local col = `col'+2
		reg years_UP initial_offerUP offerUP `controls' d_* copy_* og_* if `s'score!=. & gradeout==`gr', cluster(sasid)
			matrix lottery[`row',`col'] = _b[initial_offerUP]
			matrix lottery[`row'+1,`col'] = _se[initial_offerUP]
			matrix lottery[`row',`col'+1] = _b[offerUP]
			matrix lottery[`row'+1,`col'+1] = _se[offerUP]		

		/* 2SLS */
		local col = `col'+2
		ivreg2 `s'score (years_UP = initial_offerUP offerUP) `controls' d_* copy_* og_* if gradeout==`gr', cluster(sasid) partial(d_* `controls' copy_*)
			matrix lottery[`row',`col'] = _b[years_UP]
			matrix lottery[`row'+1,`col'] = _se[years_UP]
			matrix lottery[`row'+2,`col'] = e(N)
			matrix lottery[`row',`col'+1] = e(jp)

		local row = `row' + 3
		local col = 1
	}
}

clear
svmat lottery
br

}


***********************************************
* Lottery complier school choice              *
***********************************************
if `choice' == 1 {

use if boston_bl==1 & spedalt_bl==0 & charter_bl==0 using "${cleandata}analysisfile_2014", clear 
replace applyUP = 0 if applyUP == 1 & gradeapp == 6 & charter5 == 1
replace applyUP = 0 if applyUP == 1 & gradeapp == 7 & charter6 == 1
keep if applyUP == 1
gen up           = attend_UP6==1 & attend_UP6!=. & has_1yr_out==1 if gradeapp == 6
replace up       = attend_UP7==1 & attend_UP7!=. & has_1yr_out==1 if gradeapp == 7
gen bcharter     = bcharter6==1 & attend_UP6==0 & has_1yr_out==1 if gradeapp == 6
replace bcharter = bcharter7==1 & attend_UP7==0 & has_1yr_out==1 if gradeapp == 7
gen bpublic      = bcharter6==0 & boston6==1 & has_1yr_out==1 if gradeapp == 6
replace bpublic  = bcharter7==0 & boston7==1 & has_1yr_out==1 if gradeapp == 7
gen massoth      = bcharter6==0 & boston6==0 & has_1yr_out==1 if gradeapp == 6
replace massoth  = bcharter7==0 & boston7==0 & has_1yr_out==1 if gradeapp == 7
gen leftmass     = bcharter6==. & has_1yr_out==0 if gradeapp == 6
replace leftmass = bcharter7==. & has_1yr_out==0 if gradeapp == 7

keep if has_1yr_out==1

gen yearout = year6
qui tab yearout, gen(y_)

local con1 "if Z==0"
local con2 "if Z==1"
local con3 "[iweight=kappa0]"
local con4 "[iweight=kappa1]"
gen Z = offerUP
gen D = attend_UP6 if gradeapp == 6
replace D = attend_UP7 if gradeapp == 7
gen oneminusD = 1-D
foreach var of varlist d_* `controls' {
	drop if `var' == .
}
qui probit Z d_* `controls'
predict pZ
gen kappa_0 = (1-D)*((1-Z)-(1-pZ))/(pZ*(1-pZ))

local row = 1
foreach s in up bcharter bpublic massoth leftmass {
	summ `s' if Z==0
		matrix choice[`row',1] = r(mean)
		matrix choice[`row'+1,1] = r(N)
	summ `s' if Z==1
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
* Peer composition                            *
***********************************************
if `switching' == 1 {

	use if boston_bl==1 & spedalt_bl==0 & charter_bl==0 using "${cleandata}analysisfile_2014_wpeers", clear 
	replace applyUP = 0 if applyUP == 1 & gradeapp == 6 & charter5 == 1
	replace applyUP = 0 if applyUP == 1 & gradeapp == 7 & charter6 == 1
	keep if applyUP == 1
	gen years_UP1 = attend_UP6 if gradeapp==6
	replace years_UP1 = attend_UP7 if gradeapp==7

	foreach t in mscore escore sped lep foodfred {
		gen peer_`t'1 = peer_`t'_bl5_for6 if applyUP == 1 & gradeapp == 6
		replace peer_`t'1 = peer_`t'_bl6_for7 if applyUP == 1 & gradeapp == 7
	}
	
	gen gradeout = 6 if gradeapp==6
	replace gradeout = 7 if gradeapp==7
	qui tab gradeout, gen(og_)

	gen yearout = year6 if gradeapp==6
	replace yearout = year7 if gradeapp==7
	qui tab yearout, gen(y_)
	
	local row = 1
	local col = 1

	foreach t in sped foodfred lep mscore escore {
		
		*Summarize levels instead of gains
		qui summ peer_`t'1 if offerUP==0 & initial_offerUP==0 & peer_`t'1!=.
			matrix switching[`row',`col'] = r(mean)
			matrix switching[`row'+1,`col'] = r(sd)		

		/* 2SLS */
		local col = `col'+2
		ivreg2 peer_`t'1 (years_UP = initial_offerUP offerUP) `controls' d_* og_*, cluster(sasid) partial(d_* `controls')
			matrix switching[`row',`col'] = _b[years_UP]
			matrix switching[`row'+1,`col'] = _se[years_UP]
			matrix switching[`row'+2,`col'] = e(N)
			matrix switching[`row',`col'+1] = e(jp)

		local row = `row' + 3
		local col = 1
	}
	
	clear
	svmat switching
	br
}
