* P.HULL
* Dearborn, Harbor, and Orchard Gardens Diff-in-Diff

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
global graphs      "X:/indistrict/Graphs/"

matrix harbDiD     = J(200,12,.)
matrix dearDiD     = J(200,12,.)
matrix dearharbDiD = J(200,12,.)
matrix OGDiD       = J(200,12,.)

local dearharb = 1
local OG       = 1

***********************************************
* Harbor+Dearborn
***********************************************
if `dearharb' == 1 {

use "${simsdata}simsoct09", clear
keep if (grade=="06" | grade=="07") & (school=="00350426" | school=="00350074")
gen inharb = (school=="00350426")
gen indear = (school=="00350074")
keep sasid grade inharb indear
destring grade, replace
destring sasid, replace
reshape wide inharb indear, i(sasid) j(grade)
gen year6 = 2010
gen year7 = 2010
tempfile legacy_fa09
save `legacy_fa09', replace

use if (boston5==1 | boston6==1) & year5>=2006 using "${cleandata}dindfile_2014", clear 
egen harbcell = group(year5 masscode5 hisp5 black5 asian5 otherrace5 female5 sped5 foodfred5)
gen dearcell = harbcell
merge 1:1 sasid using `legacy_fa09', keep(master match)
drop _merge
	
foreach sch in harb dear {
	forval g = 6/7 {
		gen `sch'_elig`g' = 1 if in`sch'`g'==1 & year`g'==2010
		replace `sch'_elig`g' = 0 if in`sch'`g'==. & year`g'==2010
		bys `sch'cell: egen temp = max(`sch'_elig`g') if `sch'cell!=.
		gen `sch'samp`g' = (temp==1 & `sch'_elig`g'!=.) & bcharter5==0 & spedalt5==0
		drop temp
	}

	replace `sch'samp7 = 0 if `sch'samp6 == 1
	gen `sch'matchedsample = (`sch'samp6==1 | `sch'samp7==1)
	gen `sch'_elig = `sch'_elig6 if `sch'samp6
	replace `sch'_elig = `sch'_elig7 if `sch'samp7
	replace `sch'cell = . if `sch'matchedsample == 0
}

keep if harbmatchedsample | dearmatchedsample
gen dearharbmatchedsample = 1
egen dearharb_elig = rowmax(dear_elig harb_elig)
egen dearharbcell = group(dearcell harbcell), m
gen dearharbsamp6 = dearsamp6 | harbsamp6
gen dearharbsamp7 = dearsamp7 | harbsamp7

foreach sch in dearharb {
	preserve
	local row  = 1
	local rrow = 1
	local col  = 1

	forval out = 0/5 {	
		local outgr10 = 3+`out'
		local outgr09 = 3+`out'+1
		foreach s in m e {
			gen `s'score`out' = c_state_`s'rawsc`outgr10' if `sch'samp6==1
			replace `s'score`out' = c_state_`s'rawsc`outgr09' if `sch'samp7==1
		}
	}
	reshape long mscore escore, i(sasid) j(grade)
	replace grade = grade+3
	drop if grade==5 & `sch'samp7==1
	replace grade=5 if grade==4 & `sch'samp7==1
	replace grade=4 if grade==3 & `sch'samp7==1

	forvalues i=3/8 {
		gen dgrade`i'=(grade==`i')
		gen dgrade`i'_elig=dgrade`i'*`sch'_elig
		gen dgrade`i'_nelig=dgrade`i'*(1-`sch'_elig)
	}

	foreach s in m e {

		foreach g in 3 4 5 6 7 8 {
			areg `s'score dgrade`g' `sch'_elig dgrade`g'_elig if (grade==5 | grade==`g') & `sch'matchedsample==1, cluster(sasid) absorb(`sch'cell)
			matrix `sch'DiD[`row',`col']   = _b[dgrade`g'_elig]
			matrix `sch'DiD[`row'+1,`col'] = _se[dgrade`g'_elig]
			matrix `sch'DiD[`row'+2,`col'] = e(N)
			lincom _b[dgrade`g'_elig] + _b[`sch'_elig] + _b[dgrade`g'] + _b[_cons]
			matrix `sch'DiD[`rrow',`col'+2] = r(estimate)
			matrix `sch'DiD[`rrow',`col'+3] = r(estimate) + 1.96*r(se)
			matrix `sch'DiD[`rrow',`col'+4] = r(estimate) - 1.96*r(se)
			lincom _b[dgrade`g'] + _b[_cons]
			matrix `sch'DiD[`rrow',`col'+5] = r(estimate)
			matrix `sch'DiD[`rrow',`col'+6] = r(estimate) + 1.96*r(se)
			matrix `sch'DiD[`rrow',`col'+7] = r(estimate) - 1.96*r(se)
			matrix `sch'DiD[`rrow',`col'+8] = _b[dgrade`g'_elig]
			matrix `sch'DiD[`rrow',`col'+9] = _b[dgrade`g'_elig]+1.96*_se[dgrade`g'_elig]
			matrix `sch'DiD[`rrow',`col'+10] = _b[dgrade`g'_elig]-1.96*_se[dgrade`g'_elig]	
			local row = `row' + 3
			local ++rrow
		}
		local rrow = `rrow' + 12
	}

	clear
	svmat `sch'DiD
	keep if `sch'DiD3!=.
	replace `sch'DiD1 = 1-int((_n-1)/6)
	replace `sch'DiD2 = mod(_n-1,6)+3
	rename `sch'DiD1 math
	rename `sch'DiD2 grade

	replace grade = grade - 6	 

	forval math = 0/1 {
		if `math' == 0 {
			local stitle "ELA"
			local yclr "bg"
		}
		if `math' == 1 {
			local stitle "Math"
			local yclr "black"
		}
		twoway (line `sch'DiD6 grade if math==`math', lcolor(gs8) lwidth(thick)) (line `sch'DiD3 grade if math==`math', lcolor(black) lwidth(thick)), ///
			   legend(label(1 "Grandfathering-eligible") label(2 "Grandfathering-ineligible") cols(2) size(small)) ///
			   xline(0, lcolor(black) lpattern(dash)) scheme(s1color) ytitle(,size(small)) xtitle(" ") b2title("Grade relative to turnaround",size(small)) xlabel(,labsize(small)) ///
			   ylabel(-0.4(0.1)0.2,labsize(small)) ytitle("Standardized MCAS score", color(`yclr')) text(0.2 -1.6 "Pre-turnaround", size(small)) text(0.2 1.05   "Post-turnaround", size(small)) ///
			   text(-0.47 -1 "(Baseline)", size(small)) text(-0.47 0 "(Last legacy)", size(small)) legend(off) subtitle("`stitle'") xsize(4.5) ysize(4.25)
		graph save "${graphs}`sch'_DinD1_math`math'", replace
		graph export "${graphs}`sch'_DinD1_math`math'.wmf", replace
		
		twoway (line `sch'DiD10 grade if math==`math', lcolor(gs8) lpattern(dash)) (line `sch'DiD9 grade if math==`math', lcolor(black) lwidth(thick)) ///
			   (line `sch'DiD11 grade if math==`math', lcolor(gs8) lpattern(dash)), ///
			   legend(label(2 "GF-eligible - GF-ineligible, relative to matching grade") label(1 "95% CI") order(2 1) cols(2) size(small)) ///
			   xline(0, lcolor(black) lpattern(dash)) yline(0,lcolor(black) lpattern(shortdash)) ///
			   scheme(s1color) ytitle(,size(small)) xtitle(" ") b2title("Grade relative to turnaround",size(small)) xlabel(,labsize(small)) ///
			   ylabel(-0.4(0.1)0.4,labsize(small)) ytitle("Standardized MCAS score", color(`yclr')) text(0.4 -1.6 "Pre-turnaround", size(small)) text(0.4 1.05   "Post-turnaround", size(small)) ///
			   text(-0.495 -1 "(Baseline)", size(small)) text(-0.495 0 "(Last legacy)", size(small)) legend(off) subtitle("`stitle'") xsize(4.5) ysize(4.25)
		graph save "${graphs}`sch'_DinD2_math`math'", replace
		graph export "${graphs}`sch'_DinD2_math`math'.wmf", replace	

	}
	restore
}
}

***********************************************
* Orchard Gardens
***********************************************
if `OG' == 1 {

use "${simsdata}simsoct09", clear
keep if (grade=="06" | grade=="07") & school=="00350257"
gen inOG = 1
keep sasid grade inOG
destring grade, replace
destring sasid, replace
reshape wide inOG, i(sasid) j(grade)
gen year6 = 2010
gen year7 = 2010
tempfile legacy_fa09
save `legacy_fa09', replace

use if (boston5==1 | boston6==1) & year5>=2006 using "${cleandata}dindfile_2014", clear 
merge 1:1 sasid using `legacy_fa09', keep(master match)
drop _merge


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

	drop sch_*rawsc`bl' temp
}

egen OGcell = group(OGcell5 OGcell6), m
gen OGmatchedsample = (OGsamp6==1 | OGsamp7==1)
gen OG_elig = OG_elig6 if OGsamp6
replace OG_elig = OG_elig7 if OGsamp7
replace OGcell = . if OGmatchedsample == 0

keep if OGmatchedsample == 1

local row  = 1
local rrow = 1
local col  = 1

local sch "OG"
forval out = 0/7 {	
	local outgr10 = 1+`out'
	local outgr09 = 1+`out'+1
	foreach s in m e {
		gen `s'score`out' = c_state_`s'rawsc`outgr10' if `sch'samp6==1
		replace `s'score`out' = c_state_`s'rawsc`outgr09' if `sch'samp7==1
	}
}
reshape long mscore escore, i(sasid) j(grade)

replace grade = grade+1

forvalues i=2/8 {
	gen dgrade`i'=(grade==`i')
	gen dgrade`i'_elig=dgrade`i'*`sch'_elig
	gen dgrade`i'_nelig=dgrade`i'*(1-`sch'_elig)
}

foreach s in m e {

	foreach g in 2 3 4 5 6 7 8 {
		areg `s'score dgrade`g' `sch'_elig dgrade`g'_elig if (grade==5 | grade==`g') & `sch'matchedsample==1, cluster(sasid) absorb(`sch'cell)
		matrix `sch'DiD[`row',`col']   = _b[dgrade`g'_elig]
		matrix `sch'DiD[`row'+1,`col'] = _se[dgrade`g'_elig]
		matrix `sch'DiD[`row'+2,`col'] = e(N)
		lincom _b[dgrade`g'_elig] + _b[`sch'_elig] + _b[dgrade`g'] + _b[_cons]
		matrix `sch'DiD[`rrow',`col'+2] = r(estimate)
		matrix `sch'DiD[`rrow',`col'+3] = r(estimate) + 1.96*r(se)
		matrix `sch'DiD[`rrow',`col'+4] = r(estimate) - 1.96*r(se)
		lincom _b[dgrade`g'] + _b[_cons]
		matrix `sch'DiD[`rrow',`col'+5] = r(estimate)
		matrix `sch'DiD[`rrow',`col'+6] = r(estimate) + 1.96*r(se)
		matrix `sch'DiD[`rrow',`col'+7] = r(estimate) - 1.96*r(se)
		matrix `sch'DiD[`rrow',`col'+8] = _b[dgrade`g'_elig]
		matrix `sch'DiD[`rrow',`col'+9] = _b[dgrade`g'_elig]+1.96*_se[dgrade`g'_elig]
		matrix `sch'DiD[`rrow',`col'+10] = _b[dgrade`g'_elig]-1.96*_se[dgrade`g'_elig]	
		local row = `row' + 3
		local ++rrow
	}
	local rrow = `rrow' + 12
}

clear
svmat `sch'DiD
keep if `sch'DiD3!=.
replace `sch'DiD1 = 1-int((_n-1)/7)
replace `sch'DiD2 = mod(_n-1,7)+2
rename `sch'DiD1 math
rename `sch'DiD2 grade

replace grade = grade - 6	 

forval math = 0/1 {
	if `math' == 0 {
		local stitle "ELA"
		local yclr "bg"
	}
	if `math' == 1 {
		local stitle "Math"
		local yclr "black"
	}
	twoway (line `sch'DiD3 grade if math==`math', lcolor(black) lwidth(thick)) (line `sch'DiD6 grade if math==`math', lcolor(gs8) lwidth(thick)), ///
		   legend(label(1 "Grandfathering-eligible") label(2 "Grandfathering-ineligible") cols(2) size(small)) ///
		   xline(0, lcolor(black) lpattern(dash)) scheme(s1color) ytitle(,size(small)) xtitle(" ") b2title("Grade relative to turnaround",size(small)) xlabel(-4(1)2,labsize(small)) ///
		   ylabel(-0.6(0.2)0.4,labsize(small)) ytitle("Standardized MCAS score", color(`yclr')) text(0.4 -2 "Pre-turnaround", size(small)) text(0.4 1.05   "Post-turnaround", size(small)) ///
		   text(-0.72 -1 "(Baseline)", size(small)) text(-0.72 0 "(Last legacy)", size(small)) legend(off) subtitle("`stitle'") xsize(4.5) ysize(4.25)
	graph save "${graphs}`sch'_DinD1_math`math'", replace
	graph export "${graphs}`sch'_DinD1_math`math'.wmf", replace

	twoway (line `sch'DiD10 grade if math==`math', lcolor(gs8) lpattern(dash)) (line `sch'DiD9 grade if math==`math', lcolor(black) lwidth(thick)) ///
		   (line `sch'DiD11 grade if math==`math', lcolor(gs8) lpattern(dash)), ///
		   legend(label(2 "GF-eligible - GF-ineligible, relative to matching grade") label(1 "95% CI") order(2 1) cols(2) size(small)) ///
		   xline(0, lcolor(black) lpattern(dash)) yline(0,lcolor(black) lpattern(shortdash)) ///
		   scheme(s1color) ytitle(,size(small)) xtitle(" ") b2title("Grade relative to turnaround",size(small)) xlabel(-4(1)2,labsize(small)) ///
		   ylabel(-0.4(0.2)1,labsize(small)) ytitle("Standardized MCAS score", color(`yclr')) text(1 -2 "Pre-turnaround", size(small)) text(1 1.05   "Post-turnaround", size(small)) ///
		   text(-0.58 -1 "(Baseline)", size(small)) text(-0.58 0 "(Last legacy)", size(small)) legend(off) subtitle("`stitle'") xsize(4.5) ysize(4.25)
	graph save "${graphs}`sch'_DinD2_math`math'", replace
	graph export "${graphs}`sch'_DinD2_math`math'.wmf", replace	

}

}
