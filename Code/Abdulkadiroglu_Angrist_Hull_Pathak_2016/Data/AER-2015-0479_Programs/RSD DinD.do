* P.HULL
* RSD Grandfathering Diff-in-Diffs Analysis

clear all
cap log close
set more off
set trace off
set matsize 10000
set seed 42
adopath + "Z:\ado\plus\"

global cleandata   "X:/indistrict/Data/"

global programs "X:/indistrict/Programs/"
global graphs   "X:/indistrict/Graphs/"

local matchingdemos "female hisp black white asian otherrace"

matrix dind = J(200,200,.)

local setup = 1
local dind  = 1

local ch_to_ch = 0 // flag for including charter-to-charter takeovers

***********************************************
* STEP 1 Find grandfathering matches
***********************************************
if `setup'==1 {

	use "${cleandata}SIMS_LEAP_20082014_wide 20141030", clear


	/* restandardize to direct-runs */
	forval g = 3/8 {
			foreach t in math ela science social {
				replace leap_crsdnola_`t'sc`g' = leap_crsdnola_d_`t'sc`g'
		}
	}
	
	forval gr = 4/7 { 
		local bl_gr = `gr' - 1

		gen sch_lettergrade_adj`bl_gr' = int(sps`bl_gr'/5)
	
		egen cell`gr'   = group(spyear`bl_gr' sch_lettergrade_adj`bl_gr' `matchingdemos' sped`bl_gr' foodfred`bl_gr')

		gen gf_FL_elig`gr' = (spyear`gr' == 2010 & (fallschoolcode`gr' == "396005") & fallgrade`gr' == `gr') // Eligible for FirstLine
		gen gf_NB_elig`gr' = (spyear`gr' == 2010 & (fallschoolcode`gr' == "396040") & fallgrade`gr' == `gr') // Eligible for New Beginnings			   
		gen gf_RN_elig`gr' = (spyear`gr' == 2010 & (fallschoolcode`gr' == "396014" | fallschoolcode`gr' == "396013") & fallgrade`gr' == `gr') ///
						   | (spyear`gr' == 2011 & (fallschoolcode`gr' == "396021") & fallgrade`gr' == `gr')  ///
						   | (spyear`gr' == 2013 & (fallschoolcode`gr' == "396034" | fallschoolcode`gr' == "396203") & fallgrade`gr' == `gr') // Eligible for ReNEW		   
		gen gf_SP_elig`gr' = (spyear`gr' == 2010 & (fallschoolcode`gr' == "396030") & fallgrade`gr' == `gr') // Eligible for Spirit			   
		gen gf_FR_elig`gr' = (spyear`gr' == 2012 & (fallschoolcode`gr' == "396001") & fallgrade`gr' == `gr') // Eligible for Friends		   
		
		if `ch_to_ch' == 1 {
			gen gf_CH_elig`gr' = (spyear`gr' == 2010 & (fallschoolcode`gr' == "386001") & fallgrade`gr' == `gr') ///
							   | (spyear`gr' == 2012 & (fallschoolcode`gr' == "394003") & fallgrade`gr' == `gr') // Eligible for Choice		   
			gen gf_CC_elig`gr' = (spyear`gr' == 2011 & (fallschoolcode`gr' == "395006") & fallgrade`gr' == `gr') ///
							   | (spyear`gr' == 2013 & (fallschoolcode`gr' == "396009" | fallschoolcode`gr' == "396010") & fallgrade`gr' == `gr') // Eligible for Crescent City  
			gen gf_AR_elig`gr' = (spyear`gr' == 2013 & (fallschoolcode`gr' == "376001") & fallgrade`gr' == `gr') // Eligible for Arise		
			gen gf_CP_elig`gr' = (spyear`gr' == 2013 & (fallschoolcode`gr' == "379001") & fallgrade`gr' == `gr') // NOLA College Prep			
		}
		else {
			gen gf_CC_elig`gr' = (spyear`gr' == 2013 & (fallschoolcode`gr' == "396009" | fallschoolcode`gr' == "396010") & fallgrade`gr' == `gr') // Eligible for Crescent City  
		}
		
		egen gf_any_elig`gr' = rowmax(gf_*_elig`gr')
		replace gf_any_elig`gr' = . if (fallschoolcode`gr' == "") | (spyear`gr' == .)
		bys cell`gr': egen temp = max(gf_any_elig`gr') if cell`gr'!=.
		if `ch_to_ch' == 1 {
			gen samp`gr' = (temp == 1) & (gf_any_elig`gr' != .) & (rsd_nola`bl_gr'==1) & (alternative`bl_gr'==0)
		}
		else {
			gen samp`gr' = (temp == 1) & (gf_any_elig`gr' != .) & (charter`bl_gr'==0) & (rsd_nola`bl_gr'==1) & (alternative`bl_gr'==0)
		}
		drop temp

		foreach v of varlist cell`gr' gf_*_elig`gr' {
			replace `v' = . if samp`gr' == 0
		}
		
	}

	/* Keep the first grade sample that a student is in */
	forval gr = 4/7 {
		local grp = `gr' + 1
		forval j = `grp'/7 {
			replace samp`j' = 0 if samp`gr' == 1
			foreach v of varlist cell`j' gf_*_elig`j' {
				replace `v' = 0 if samp`j' == 0
			}
		}
	}

	egen matchedsample =  rowtotal(samp*)
	assert matchedsample <= 1
	keep if matchedsample == 1

	if `ch_to_ch' == 1 {
	foreach v in FL NB RN SP FR CH CC AR CP any {
		egen gf_`v'_elig = rowmax(gf_`v'_elig*)
	}
	}
	else {
	foreach v in FL NB RN SP FR CC any {
		egen gf_`v'_elig = rowmax(gf_`v'_elig*)
	}	
	}
	
	egen match_cell = group(cell*), missing

	forval i = 1/3 {
		gen mscore`i' = .
		gen escore`i' = .
		gen cscore`i' = .
		gen sscore`i' = .
		
		forval gr = 4/7 { 
			local outgr = `gr' + `i'
			foreach t in m e c s {
				if "`t'" == "m" {
					local test "math"
				}
				if "`t'" == "e" {
					local test "ela"
				}
				if "`t'" == "c" {
					local test "science"
				}
				if "`t'" == "s" {
					local test "social"
				}
				replace `t'score`i' = leap_crsdnola_`test'sc`outgr' if samp`gr' == 1
			}
		}
	}
	forval i = 0/4 {
		gen mscore9`i' = .
		gen escore9`i' = .
		gen cscore9`i' = .
		gen sscore9`i' = .
	
		forval gr = 4/7 { 
			local outgr = `gr' - `i'
			foreach t in m e c s {
				if "`t'" == "m" {
					local test "math"
				}
				if "`t'" == "e" {
					local test "ela"
				}
				if "`t'" == "c" {
					local test "science"
				}
				if "`t'" == "s" {
					local test "social"
				}
				cap replace `t'score9`i' = leap_crsdnola_`test'sc`outgr' if samp`gr' == 1
			}
		}
	}
	
	gen bl_sped = .
	gen bl_lep = .
	gen bl_foodfred = .
	gen bl_m = .
	gen bl_e = .
	gen bl_c = .
	gen bl_s = .
	gen bl_charter = .
	forval gr = 4/7 {
		local bl_gr = `gr' - 1
		replace bl_sped = sped`bl_gr' if samp`gr' == 1
		replace bl_lep = lep`bl_gr' if samp`gr' == 1
		replace bl_foodfred = foodfred`bl_gr' if samp`gr' == 1
		replace bl_m = leap_crsdnola_mathsc`bl_gr' if samp`gr' == 1
		replace bl_e = leap_crsdnola_elasc`bl_gr' if samp`gr' == 1
		replace bl_c = leap_crsdnola_sciencesc`bl_gr' if samp`gr' == 1
		replace bl_s = leap_crsdnola_socialsc`bl_gr' if samp`gr' == 1		
		replace bl_charter = charter`bl_gr' if samp`gr' == 1
	}

	keep sasid mscore* escore* cscore* sscore* gf_*elig* *cell* bl*_* *samp* spyear3 spyear4 spyear5 spyear6 spyear7 spyear8 schoolcode3 schoolcode4 schoolcode5 schoolcode6 schoolcode7 schoolcode8 leap_crsdnola_mathsc3 leap_crsdnola_elasc3 leap_crsdnola_sciencesc3 leap_crsdnola_socialsc3 hisp* black* white* asian* female* foodfred3 foodfred4 foodfred5 foodfred6 foodfred7 foodfred8 sped3 sped4 sped5 sped6 sped7 sped8 lep3 lep4 lep5 lep6 lep7 lep8
	compress

	save "${cleandata}NOLA_all_dd2014", replace
}

if `dind' == 1 {
use "${cleandata}NOLA_all_dd2014", clear

qui tab match_cell, gen(dg_)
gen samp = .
forval gr = 4/7 {
	replace samp = `gr' if samp`gr'==1
}

reshape long mscore escore cscore sscore, i(sasid samp) j(time)
replace time = 90 - time if time >= 90
replace time = time + 1

local row  = 1
local rrow = 1
local col  = 1

foreach s in m e c s {

local f_row = `row'

	forvalues i = -3/4 {
		gen dtime = (time==`i')
		gen dtime_elig   = dtime * gf_any_elig
		local col = 1
		cap {
		areg `s'score dtime gf_any_elig dtime_elig if time==0 | time==`i', cluster(sasid) absorb(match_cell)
			matrix dind[`row',`col']   = _b[dtime_elig]
			matrix dind[`row'+1,`col'] = _se[dtime_elig]
			matrix dind[`row'+2,`col'] = e(N)
			lincom _b[dtime_elig] + _b[gf_any_elig] + _b[dtime] + _b[_cons]
			matrix dind[`rrow',`col'+2]   = r(estimate)
			matrix dind[`rrow',`col'+3] = r(estimate) + 1.96*r(se)
			matrix dind[`rrow',`col'+4] = r(estimate) - 1.96*r(se)
			lincom _b[dtime] + _b[_cons]
			matrix dind[`rrow',`col'+5]   = r(estimate)
			matrix dind[`rrow',`col'+6] = r(estimate) + 1.96*r(se)
			matrix dind[`rrow',`col'+7] = r(estimate) - 1.96*r(se)
			matrix dind[`rrow',`col'+8] = _b[dtime_elig]
			matrix dind[`rrow',`col'+9] = _b[dtime_elig]+1.96*_se[dtime_elig]
			matrix dind[`rrow',`col'+10] = _b[dtime_elig]-1.96*_se[dtime_elig]	
		}
		local row = `row' + 3
		local ++rrow
		drop dtime*
		}
		local rrow = `rrow' + 18
	}

clear
svmat dind
keep dind1-dind11
keep if dind3!=.
replace dind1 = int((_n-1)/8)
replace dind2 = mod(_n-1,8)-4
drop if dind2 >= 3 | dind2 <= -4

rename dind1 test
rename dind2 grade


forval test = 0/3 {
	if `test' == 0 {
		local stitle "Math"
		local yclr "black"
	}
	if `test' == 1 {
		local stitle "ELA"
		local yclr "bg"
	}
	if `test' == 2 {
		local stitle "Science"
		local yclr "black"
	}
	if `test' == 3 {
		local stitle "Social science"
		local yclr "bg"
	}
	twoway (line dind6 grade if test==`test', lcolor(gs8) lwidth(thick)) (line dind3 grade if test==`test', lcolor(black) lwidth(thick)) , ///
		   legend(label(1 "Grandfathering-eligible") label(2 "Grandfathering-ineligible") cols(2) size(small)) ///
		   xline(0, lcolor(black) lpattern(dash)) scheme(s1color) ytitle(,size(small)) xtitle(" ") b2title("Grade relative to takeover",size(small)) xlabel(,labsize(small)) /// 
		   ylabel(-0.1(0.1)0.5,labsize(small)) ytitle("Standardized LEAP/iLEAP score", color(`yclr')) text(0.5 -1.6 "Pre-charter takeover", size(small)) text(0.5 1.05   "Post-charter takeover", size(small)) ///
		   text(-0.17 -1 "(Baseline)", size(small)) text(-0.17 0 "(Last legacy)", size(small)) legend(off) subtitle("`stitle'") xsize(4.5) ysize(4.25)

	graph save "${graphs}RSD_DinD1_math`test'", replace
	graph export "${graphs}RSD_DinD1_math`test'.wmf", replace

replace dind11 = -0.4 if dind11<-0.4
replace dind10 = 0.4 if dind10>0.4
	twoway (line dind10 grade if test==`test', lcolor(gs8) lpattern(dash)) (line dind9 grade if test==`test', lcolor(black) lwidth(thick)) ///
		   (line dind11 grade if test==`test', lcolor(gs8) lpattern(dash)), ///
		   legend(label(2 "GF-eligible - GF-ineligible, relative to baseline") label(1 "95% CI") order(2 1) cols(2) size(small)) ///
		   xline(0, lcolor(black) lpattern(dash)) yline(0,lcolor(black) lpattern(shortdash)) ///
		   scheme(s1color) ytitle(,size(small)) xtitle(,size(small)) xtitle(" ") b2title("Grade relative to takeover",size(small)) xlabel(,labsize(small)) ///
		   ylabel(-0.4(0.1)0.4,labsize(small)) ytitle("Standardized LEAP/iLEAP score", color(`yclr')) text(0.4 -1.6 "Pre-charter takeover", size(small)) text(0.4 1.05   "Post-charter takeover", size(small)) ///
		   text(-0.495 -1 "(Baseline)", size(small)) text(-0.495 0 "(Last legacy)", size(small)) legend(off) subtitle("`stitle'") xsize(4.5) ysize(4.25)
	graph save "${graphs}RSD_DinD2_math`test'", replace
	graph export "${graphs}RSD_DinD2_math`test'.wmf", replace	
	
}

clear
svmat dind
br
}
