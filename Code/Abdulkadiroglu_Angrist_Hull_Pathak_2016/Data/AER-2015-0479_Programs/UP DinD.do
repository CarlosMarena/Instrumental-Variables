* P.HULL
* UP Grandfathering Diff-in-Diffs Analysis

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

global programs    "X:/indistrict/Programs/"
global graphs      "X:/indistrict/Graphs/"

local simsdate "10 30 2014"

local covstems "female black hisp asian otherrace sped lep foodfred fem_min"

matrix DiD = J(200,12,.)

local setup = 1
local DiD   = 1

***********************************************
* Match applicants to State data
***********************************************
if `setup'==1 {
use "${cleandata}simsmcaswide `simsdate'", clear

forval g = 3/8 {
		foreach t in m e {
			qui replace c_state_`t'rawsc`g' = c_bps_`t'rawsc`g'
		}
}


gen baselinegrade = 5

foreach v in boston year `covstems' c_state_mrawsc c_state_erawsc {
	gen `v'_bl = `v'5
}

gen white_bl = 1-(hisp_bl + black_bl + asian_bl + otherrace_bl)

gen hasbaselinedemos = !missing(female_bl, black_bl, hisp_bl, asian_bl, otherrace_bl, sped_bl, lep_bl, foodfred_bl)

compress
save "${cleandata}dindfile_2014", replace	
}

***********************************************
* Differences-in-differences
***********************************************
if `DiD' == 1 {

local row  = 1
local rrow = 1
local col  = 1

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

use if boston5==1 & spedalt5==0 & year5>=2006 & hasbaselinedemos using "${cleandata}dindfile_2014", clear 
egen cell = group(masscode5 year5 hisp5 black5 asian5 otherrace5 female5 sped5 foodfred5)

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

forval out = 0/5 {	
	local outgr10 = 3+`out'
	local outgr09 = 3+`out'+1
	foreach s in m e {
		gen `s'score`out' = c_state_`s'rawsc`outgr10' if samp6==1
		replace `s'score`out' = c_state_`s'rawsc`outgr09' if samp7==1
	}
}
reshape long mscore escore, i(sasid) j(grade)
replace grade = grade+3
drop if grade==5 & samp7==1
replace grade=5 if grade==4 & samp7==1
replace grade=4 if grade==3 & samp7==1

forvalues i=3/8 {
	gen dgrade`i'=(grade==`i')
	gen dgrade`i'_elig=dgrade`i'*gf_elig
	gen dgrade`i'_nelig=dgrade`i'*(1-gf_elig)
}

foreach s in m e {

	foreach g in 3 4 5 6 7 8 {
		areg `s'score dgrade`g' gf_elig dgrade`g'_elig if grade==5 | grade==`g', cluster(sasid) absorb(cell)
		matrix DiD[`row',`col']   = _b[dgrade`g'_elig]
		matrix DiD[`row'+1,`col'] = _se[dgrade`g'_elig]
		matrix DiD[`row'+2,`col'] = e(N)
		lincom _b[dgrade`g'_elig] + _b[gf_elig] + _b[dgrade`g'] + _b[_cons]
		matrix DiD[`rrow',`col'+2] = r(estimate)
		matrix DiD[`rrow',`col'+3] = r(estimate) + 1.96*r(se)
		matrix DiD[`rrow',`col'+4] = r(estimate) - 1.96*r(se)
		lincom _b[dgrade`g'] + _b[_cons]
		matrix DiD[`rrow',`col'+5] = r(estimate)
		matrix DiD[`rrow',`col'+6] = r(estimate) + 1.96*r(se)
		matrix DiD[`rrow',`col'+7] = r(estimate) - 1.96*r(se)
		matrix DiD[`rrow',`col'+8] = _b[dgrade`g'_elig]
		matrix DiD[`rrow',`col'+9] = _b[dgrade`g'_elig]+1.96*_se[dgrade`g'_elig]
		matrix DiD[`rrow',`col'+10] = _b[dgrade`g'_elig]-1.96*_se[dgrade`g'_elig]	
		local row = `row' + 3
		local ++rrow
	}
	local rrow = `rrow' + 12
}

clear
svmat DiD
keep if DiD3!=.
replace DiD1 = 1-int((_n-1)/6)
replace DiD2 = mod(_n-1,6)+3
rename DiD1 math
rename DiD2 grade

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

	twoway (line DiD6 grade if math==`math', lcolor(gs8) lwidth(thick)) (line DiD3 grade if math==`math', lcolor(black) lwidth(thick)), ///
		   legend(label(1 "Grandfathering-eligible") label(2 "Grandfathering-ineligible") cols(2) size(small)) ///
		   xline(0, lcolor(black) lpattern(dash)) scheme(s1color) ytitle(,size(small)) xtitle(" ") b2title("Grade relative to takeover",size(small)) xlabel(,labsize(small)) ///
		   ylabel(-0.4(0.1)0.6,labsize(small)) ytitle("Standardized MCAS score", color(`yclr')) text(0.6 -1.6 "Pre-UP takeover", size(small)) text(0.6 1.05   "Post-UP takeover", size(small)) ///
		   text(-0.52 -1 "(Baseline)", size(small)) text(-0.52 0 "(Last legacy)", size(small)) legend(off) subtitle("`stitle'") xsize(4.5) ysize(4.25)
	graph save "${graphs}UP_DinD1_math`math'", replace
	graph export "${graphs}UP_DinD1_math`math'.wmf", replace
	
	twoway (line DiD10 grade if math==`math', lcolor(gs8) lpattern(dash)) (line DiD9 grade if math==`math', lcolor(black) lwidth(thick)) ///
		   (line DiD11 grade if math==`math', lcolor(gs8) lpattern(dash)), ///
		   legend(label(2 "GF-eligible - GF-ineligible, relative to baseline") label(1 "95% CI") order(2 1) cols(2) size(small)) ///
		   xline(0, lcolor(black) lpattern(dash)) yline(0,lcolor(black) lpattern(shortdash)) ///
		   scheme(s1color) ytitle(,size(small)) xtitle(" ") b2title("Grade relative to takeover",size(small)) xlabel(,labsize(small)) ///
		   ylabel(-0.4(0.2)1,labsize(small)) ytitle("Standardized MCAS score", color(`yclr')) text(1 -1.6 "Pre-UP takeover", size(small)) text(1 1.05   "Post-UP takeover", size(small)) ///
		   text(-0.58 -1 "(Baseline)", size(small)) text(-0.58 0 "(Last legacy)", size(small)) legend(off) subtitle("`stitle'") xsize(4.5) ysize(4.25)
	graph save "${graphs}UP_DinD2_math`math'", replace
	graph export "${graphs}UP_DinD2_math`math'.wmf", replace	

}

clear
svmat DiD
br
}
