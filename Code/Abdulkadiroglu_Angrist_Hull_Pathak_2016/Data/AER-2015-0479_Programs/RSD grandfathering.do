* P.HULL
* RSD Grandfathering Analysis

clear all
cap log close
set more off
set trace off
set matsize 10000
set seed 42
adopath + "Z:\ado\plus\"

global nolasims    "X:/NOLA/Combined data/sims_clean/"
global cleandata   "X:/indistrict/Data/"
global closures    "X:/indistrict/Tables/stata output/output/"
global types       "X:/indistrict/Tables/stata output/"

global programs "X:/indistrict/Programs/"
global graphs   "X:/indistrict/Graphs/"

local matchingdemos "female hisp black white asian otherrace"
local covs          "bl_sped bl_foodfred bl_lep bl_m bl_e bl_c bl_s y_*"

local simsdate 20141030

matrix statusgrowth   = J(200,100,.)
matrix balance        = J(200,100,.)
matrix gfing          = J(200,5,.)
matrix gfingnoFD      = J(200,5,.)
matrix choice         = J(100,10,.)
matrix fallbackchange = J(200,100,.)
matrix persistence    = J(200,100,.)
matrix cohorts        = J(200,100,.)
matrix switching      = J(200,100,.)

local setup           = 1
local switchingsetup  = 1

local statusgrowth    = 1 /* Fig 1                    */
local balance         = 1 /* Tab 2; Tab B3 (Panel A)  */
local gfing           = 1 /* Tab 3; Tab B4            */
local choice          = 1 /* Tab 4                    */
local fallbackchange  = 1 /* Tab 5                    */
local switching       = 1 /* Tab 10 (Column 1)        */
local persistence     = 1 /* Tab A3                   */	
local fallbackgraph   = 1 /* Fig B3                   */
local cohorts         = 1 /* Tab B2 (Panel A)         */
local subgroup        = 1 /* Tab B5                   */
local altfbchange     = 1 /* Tab B6                   */

***********************************************
* Construct main RSD GF analysis file         *
***********************************************
if `setup'==1 {

	use "${cleandata}SIMS_LEAP_20082014_wide `simsdate'", clear

	forval gr = 4/7 { 
		local bl_gr = `gr' - 1

		gen sch_lettergrade_adj`bl_gr' = int(sps`bl_gr'/5)
		
		egen cell`gr'   = group(spyear`bl_gr' sch_lettergrade_adj`bl_gr' `matchingdemos' sped`bl_gr' foodfred`bl_gr')   

		gen gf_FL_elig`gr' = (spyear`gr' == 2010 & (fallschoolcode`gr' == "396005") & fallgrade`gr' == `gr') // Eligible for FirstLine
		gen gf_NB_elig`gr' = (spyear`gr' == 2010 & (fallschoolcode`gr' == "396040") & fallgrade`gr' == `gr') // Eligible for New Beginnings						   
		
		gen gf_RN_elig`gr' = (spyear`gr' == 2010 & (fallschoolcode`gr' == "396014" /*| fallschoolcode`gr' == "396013"*/) & fallgrade`gr' == `gr') ///
						   | (spyear`gr' == 2011 & (fallschoolcode`gr' == "396021") & fallgrade`gr' == `gr')  ///
						   | (spyear`gr' == 2013 & (fallschoolcode`gr' == "396034" | fallschoolcode`gr' == "396203") & fallgrade`gr' == `gr') // Eligible for ReNEW
	
		gen gf_SP_elig`gr' = (spyear`gr' == 2010 & (fallschoolcode`gr' == "396030") & fallgrade`gr' == `gr') // Eligible for Spirit			   
		gen gf_FR_elig`gr' = (spyear`gr' == 2012 & (fallschoolcode`gr' == "396001") & fallgrade`gr' == `gr') // Eligible for Friends  	
		gen gf_CC_elig`gr' = (spyear`gr' == 2013 & (fallschoolcode`gr' == "396009" | fallschoolcode`gr' == "396010") & fallgrade`gr' == `gr') // Eligible for Crescent City  
				
		egen gf_any_elig`gr' = rowmax(gf_*_elig`gr')
		replace gf_any_elig`gr' = . if (fallschoolcode`gr' == "") | (spyear`gr' == .)
		bys cell`gr': egen temp = max(gf_any_elig`gr') if cell`gr'!=.
		gen samp`gr' = (temp == 1) & (gf_any_elig`gr' != .) & (charter`bl_gr'==0) & (rsd_nola`bl_gr'==1) & (alternative`bl_gr'==0)
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

	foreach v in FL NB /*RN*/ SP FR CC any {
		egen gf_`v'_elig = rowmax(gf_`v'_elig*)
	}
	
	egen match_cell = group(cell*), missing

	forval i = 0/4 {
		foreach t in m e c s {
			gen has_`t'score`i' = .
			gen `t'score`i' = .
			forval j = 1/5 {
				gen `t'status`j'_`i' = .
			}
			gen `t'gfscore`i' = .
		}
		
		foreach s in FL NB RN SP FR CC {
			gen attend_`s'`i' = .
		}
		
		gen repeatout`i' = .
		gen incharter`i' = .
		forval gr = 4/7 { 
			local outgr = `gr' + `i'
			local lgr = `gr' - 1
			
			if `outgr' <= 8 {
				
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
					replace has_`t'score`i' = (leap_crsdnola_`test'sc`outgr'!=.) if samp`gr' == 1 & spyear`lgr'+1+`i'<=2014
					replace `t'gfscore`i' = leap_crsdnola_`test'sc`gr' if samp`gr' == 1	
					forval j = 1/5 {
						replace `t'status`j'_`i' = leap_`test'_achieve`j'`outgr' if samp`gr' == 1
						replace `t'status`j'_`i' = leap_`test'_achieve`j'`outgr'  if samp`gr' == 1
					}
				}

				replace attend_FL`i' = 1 if (schoolcode`outgr' == "399004") & samp`gr' == 1 & schoolcode`outgr'!=""
				replace attend_FL`i' = 0 if (schoolcode`outgr' != "399004") & samp`gr' == 1 & schoolcode`outgr'!=""

				replace attend_NB`i' = 1 if (schoolcode`outgr' == "300004") & samp`gr' == 1 & schoolcode`outgr'!=""
				replace attend_NB`i' = 0 if (schoolcode`outgr' != "300004") & samp`gr' == 1 & schoolcode`outgr'!=""				

				replace attend_RN`i' = 1 if (schoolcode`outgr' == "369001" | schoolcode`outgr' == "369002" | schoolcode`outgr' == "369003" | schoolcode`outgr' == "369006") & samp`gr' == 1 & schoolcode`outgr'!=""
				replace attend_RN`i' = 0 if (schoolcode`outgr' != "369001" & schoolcode`outgr' != "369002" & schoolcode`outgr' != "369003" & schoolcode`outgr' != "369006") & samp`gr' == 1 & schoolcode`outgr'!=""

				replace attend_SP`i' = 1 if (schoolcode`outgr' == "367001") & samp`gr' == 1 & schoolcode`outgr'!=""
				replace attend_SP`i' = 0 if (schoolcode`outgr' != "367001") & samp`gr' == 1 & schoolcode`outgr'!=""

				replace attend_FR`i' = 1 if (schoolcode`outgr' == "391002") & samp`gr' == 1 & schoolcode`outgr'!=""
				replace attend_FR`i' = 0 if (schoolcode`outgr' != "391002") & samp`gr' == 1 & schoolcode`outgr'!=""

				replace attend_CC`i' = 1 if (schoolcode`outgr' == "363002") & samp`gr' == 1 & schoolcode`outgr'!=""
				replace attend_CC`i' = 0 if (schoolcode`outgr' != "363002") & samp`gr' == 1 & schoolcode`outgr'!=""
					
				replace repeatout`i' = repeat`outgr' if samp`gr' == 1
				replace incharter`i' = charter`outgr' if samp`gr' == 1
				
			}
		}
		egen attend_any`i' = rowmax(attend_*`i')	
	}

	foreach v in FL NB RN SP FR CC any {
		gen years_`v'1 = attend_`v'1
		gen temp = attend_`v'1*(1+repeatout1)
		egen years_`v'2 = rowtotal(temp attend_`v'2)
		replace years_`v'2 = . if attend_`v'1==. & attend_`v'2==.	
		gen temp2 = attend_`v'2*(1+repeatout2)
		egen years_`v'3 = rowtotal(temp temp2 attend_`v'3)
		replace years_`v'3 = . if attend_`v'1==. & attend_`v'2==. & attend_`v'3==.	
		gen temp3 = attend_`v'3*(1+repeatout3)
		egen years_`v'4 = rowtotal(temp temp2 temp3 attend_`v'4)
		replace years_`v'4 = . if attend_`v'1==. & attend_`v'2==. & attend_`v'3==. & attend_`v'4==.			
		drop temp temp2 temp3
	}	
	
	foreach v in sped foodfred lep {
		gen bl_`v' = .
		forval gr = 4/7 {
			local bl_gr = `gr' - 1
			replace bl_`v' = `v'`bl_gr' if samp`gr' == 1
		}
	}
	gen bl_m = .
	gen bl_e = .
	gen bl_c = .
	gen bl_s = .
	gen bl_charter = .
	gen bl_year = .
	forval gr = 4/7 {
		local bl_gr = `gr' - 1
		replace bl_m = leap_crsdnola_mathsc`bl_gr' if samp`gr' == 1
		replace bl_e = leap_crsdnola_elasc`bl_gr' if samp`gr' == 1	
		replace bl_c = leap_crsdnola_sciencesc`bl_gr' if samp`gr' == 1
		replace bl_s = leap_crsdnola_socialsc`bl_gr' if samp`gr' == 1	
		replace bl_charter = charter`bl_gr' if samp`gr' == 1
		replace bl_year = spyear`bl_gr' if samp`gr' == 1
	}
	egen temp = rowtotal(samp*)

	forval out = 1/4 {
		gen year`out' = .
		gen school`out' = ""
		forval gr = 4/7 {
			local outgr = `gr' + `out'
			replace year`out'   = spyear`outgr' if samp`gr'==1
			replace school`out' = schoolcode`outgr' if samp`gr'==1
		}
	}

	qui tab match_cell, gen(d_)

	keep sasid mscore* mgfscore* *has_mscore* leap_math_achieve* escore* egfscore* leap_ela_achieve* *has_escore* cscore* cgfscore* *has_cscore* sscore* sgfscore* *has_sscore* incharter* mstatus* estatus* cstatus* sstatus* *attend_* years_* gf_*elig* *cell* d_* bl_* *samp* spyear3 spyear4 spyear5 spyear6 spyear7 spyear8 rsd_nola3 rsd_nola4 rsd_nola5 rsd_nola6 rsd_nola7 rsd_nola8 alternative3 alternative4 alternative5 alternative6 alternative7 alternative8 schoolcode3 *schoolcode4 *schoolcode5 *schoolcode6 *schoolcode7 *schoolcode8 *schoolname4 *schoolname5 *schoolname6 *schoolname7 *schoolname8 fallgrade* leap_crsdnola_mathsc* leap_crsdnola_elasc* leap_crsdnola_sciencesc* leap_crsdnola_socialsc* hisp* black* white* otherrace* asian* female* foodfred3 foodfred4 foodfred5 foodfred6 foodfred7 foodfred8 sped3 sped4 sped5 sped6 sped7 sped8 lep3 lep4 lep5 lep6 lep7 lep8 repeat3 repeat4 repeat5 repeat6 repeat7 repeat8 charter3 charter4 charter5 charter6 charter7 charter8 years_till_take* repeatout* sps*
	compress

	save "${cleandata}NOLA_rsdbl_gf2014", replace
}

***********************************************
* Compute peer quality                        *
***********************************************
if `switchingsetup'==1 {
	use "${cleandata}SIMS_LEAP_20082014_wide `simsdate'", clear
	keep sasid spyear5 schoolcode5 spyear6 schoolcode6 spyear7 schoolcode7 spyear8 schoolcode8 ///
		leap_crsdnola_*sc3 leap_crsdnola_*sc4 leap_crsdnola_*sc5 leap_crsdnola_*sc6 ///
		sped3 sped4 sped5 sped6 lep3 lep4 lep5 lep6 foodfred3 foodfred4 foodfred5 foodfred6
	forval outgr = 5/8 {
		preserve
		collapse (mean) leap_crsdnola_*sc3 leap_crsdnola_*sc4 leap_crsdnola_*sc5 leap_crsdnola_*sc6 sped3 sped4 sped5 sped6 lep3 lep4 lep5 lep6 foodfred3 foodfred4 foodfred5 foodfred6, by(spyear`outgr' schoolcode`outgr')
		forval blgr = 3/6 {
			gen peer_mscore`blgr'_for`outgr' = leap_crsdnola_mathsc`blgr'
			gen peer_escore`blgr'_for`outgr' = leap_crsdnola_elasc`blgr' 
			gen peer_cscore`blgr'_for`outgr' = leap_crsdnola_sciencesc`blgr'
			gen peer_sscore`blgr'_for`outgr' = leap_crsdnola_socialsc`blgr'
			gen peer_sped`blgr'_for`outgr' = sped`blgr'
			gen peer_lep`blgr'_for`outgr' = lep`blgr'
			gen peer_foodfred`blgr'_for`outgr' = foodfred`blgr'
		}
		drop if spyear`outgr' == . | schoolcode`outgr' == ""
		keep schoolcode`outgr' spyear`outgr' peer_*
		tempfile peer_for`outgr'
		save `peer_for`outgr'', replace
		restore
	}
	use "${cleandata}NOLA_rsdbl_gf2014", clear
	forval outgr = 5/8 {
		merge m:1 spyear`outgr' schoolcode`outgr' using `peer_for`outgr''
		drop _merge
	}
	foreach t in mscore escore cscore sscore sped lep foodfred {
		forval out = 1/4 {
			gen peer_`t'`out' = .
			forval samp = 4/7 {
				local blgr  = `samp' - 1
				local outgr = `samp' + `out'
				cap replace peer_`t'`out' = peer_`t'`blgr'_for`outgr' if samp`samp' == 1
			}
		}
	}
	drop peer_*_*
	save "${cleandata}NOLA_rsdbl_gf2014_wpeers", replace
}
***********************************************
* Describe grandfathering cohorts             *
***********************************************
if `cohorts' == 1 {
	use "${cleandata}NOLA_rsdbl_gf2014", clear

	foreach stub in R wB iD iS {
		foreach sch in dil lau liv har gen sar jos hab1 hab2 sch1 sch2 {
			gen gf_`sch'_elig`stub' = .
		}
	}
	
	forval g = 4/7 {
		  
		local bl_gr = `g' - 1
		gen sch_lettergrade_adj`bl_gr' = int(sps`bl_gr'/5)
		egen acell`g'   = group(spyear`bl_gr' sch_lettergrade_adj`bl_gr' `matchingdemos' sped`bl_gr')

		replace gf_dil_eligR = 1 if (spyear`g' == 2010 & fallschoolcode`g' == "396005" & fallgrade`g' == `g')
		replace gf_lau_eligR = 1 if (spyear`g' == 2010 & fallschoolcode`g' == "396013" & fallgrade`g' == `g')
		replace gf_liv_eligR = 1 if (spyear`g' == 2010 & fallschoolcode`g' == "396014" & fallgrade`g' == `g')
		replace gf_har_eligR = 1 if (spyear`g' == 2010 & fallschoolcode`g' == "396030" & fallgrade`g' == `g')
		replace gf_gen_eligR = 1 if (spyear`g' == 2010 & fallschoolcode`g' == "396040" & fallgrade`g' == `g')
		replace gf_sar_eligR = 1 if (spyear`g' == 2011 & fallschoolcode`g' == "396021" & fallgrade`g' == `g')
		replace gf_jos_eligR = 1 if (spyear`g' == 2012 & fallschoolcode`g' == "396001" & fallgrade`g' == `g')	
		replace gf_hab1_eligR= 1 if (spyear`g' == 2013 & fallschoolcode`g' == "396034" & fallgrade`g' == `g')
		replace gf_hab2_eligR= 1 if (spyear`g' == 2013 & fallschoolcode`g' == "396203" & fallgrade`g' == `g')
		replace gf_sch1_eligR= 1 if (spyear`g' == 2013 & fallschoolcode`g' == "396009" & fallgrade`g' == `g')
		replace gf_sch2_eligR= 1 if (spyear`g' == 2013 & fallschoolcode`g' == "396010" & fallgrade`g' == `g')
		
		replace gf_dil_eligwB = 1 if (spyear`g' == 2010 & fallschoolcode`g' == "396005" & fallgrade`g' == `g' & acell`g'!=.)
		replace gf_lau_eligwB = 1 if (spyear`g' == 2010 & fallschoolcode`g' == "396013" & fallgrade`g' == `g' & acell`g'!=.)
		replace gf_liv_eligwB = 1 if (spyear`g' == 2010 & fallschoolcode`g' == "396014" & fallgrade`g' == `g' & acell`g'!=.)
		replace gf_har_eligwB = 1 if (spyear`g' == 2010 & fallschoolcode`g' == "396030" & fallgrade`g' == `g' & acell`g'!=.)
		replace gf_gen_eligwB = 1 if (spyear`g' == 2010 & fallschoolcode`g' == "396040" & fallgrade`g' == `g' & acell`g'!=.)
		replace gf_sar_eligwB = 1 if (spyear`g' == 2011 & fallschoolcode`g' == "396021" & fallgrade`g' == `g' & acell`g'!=.)
		replace gf_jos_eligwB = 1 if (spyear`g' == 2012 & fallschoolcode`g' == "396001" & fallgrade`g' == `g' & acell`g'!=.)	
		replace gf_hab1_eligwB= 1 if (spyear`g' == 2013 & fallschoolcode`g' == "396034" & fallgrade`g' == `g' & acell`g'!=.)
		replace gf_hab2_eligwB= 1 if (spyear`g' == 2013 & fallschoolcode`g' == "396203" & fallgrade`g' == `g' & acell`g'!=.)
		replace gf_sch1_eligwB= 1 if (spyear`g' == 2013 & fallschoolcode`g' == "396009" & fallgrade`g' == `g' & acell`g'!=.)
		replace gf_sch2_eligwB= 1 if (spyear`g' == 2013 & fallschoolcode`g' == "396010" & fallgrade`g' == `g' & acell`g'!=.)

		replace gf_dil_eligiD = 1 if (spyear`g' == 2010 & fallschoolcode`g' == "396005" & fallgrade`g' == `g' & charter`bl_gr'==0 & rsd_nola`bl_gr'==1 & alternative`bl_gr'==0 & acell`g'!=.)
		replace gf_lau_eligiD = 1 if (spyear`g' == 2010 & fallschoolcode`g' == "396013" & fallgrade`g' == `g' & charter`bl_gr'==0 & rsd_nola`bl_gr'==1 & alternative`bl_gr'==0 & acell`g'!=.)
		replace gf_liv_eligiD = 1 if (spyear`g' == 2010 & fallschoolcode`g' == "396014" & fallgrade`g' == `g' & charter`bl_gr'==0 & rsd_nola`bl_gr'==1 & alternative`bl_gr'==0 & acell`g'!=.)
		replace gf_har_eligiD = 1 if (spyear`g' == 2010 & fallschoolcode`g' == "396030" & fallgrade`g' == `g' & charter`bl_gr'==0 & rsd_nola`bl_gr'==1 & alternative`bl_gr'==0 & acell`g'!=.)
		replace gf_gen_eligiD = 1 if (spyear`g' == 2010 & fallschoolcode`g' == "396040" & fallgrade`g' == `g' & charter`bl_gr'==0 & rsd_nola`bl_gr'==1 & alternative`bl_gr'==0 & acell`g'!=.)
		replace gf_sar_eligiD = 1 if (spyear`g' == 2011 & fallschoolcode`g' == "396021" & fallgrade`g' == `g' & charter`bl_gr'==0 & rsd_nola`bl_gr'==1 & alternative`bl_gr'==0 & acell`g'!=.)
		replace gf_jos_eligiD = 1 if (spyear`g' == 2012 & fallschoolcode`g' == "396001" & fallgrade`g' == `g' & charter`bl_gr'==0 & rsd_nola`bl_gr'==1 & alternative`bl_gr'==0 & acell`g'!=.)	
		replace gf_hab1_eligiD= 1 if (spyear`g' == 2013 & fallschoolcode`g' == "396034" & fallgrade`g' == `g' & charter`bl_gr'==0 & rsd_nola`bl_gr'==1 & alternative`bl_gr'==0 & acell`g'!=.)
		replace gf_hab2_eligiD= 1 if (spyear`g' == 2013 & fallschoolcode`g' == "396203" & fallgrade`g' == `g' & charter`bl_gr'==0 & rsd_nola`bl_gr'==1 & alternative`bl_gr'==0 & acell`g'!=.)
		replace gf_sch1_eligiD= 1 if (spyear`g' == 2013 & fallschoolcode`g' == "396009" & fallgrade`g' == `g' & charter`bl_gr'==0 & rsd_nola`bl_gr'==1 & alternative`bl_gr'==0 & acell`g'!=.)
		replace gf_sch2_eligiD= 1 if (spyear`g' == 2013 & fallschoolcode`g' == "396010" & fallgrade`g' == `g' & charter`bl_gr'==0 & rsd_nola`bl_gr'==1 & alternative`bl_gr'==0 & acell`g'!=.)

		replace gf_dil_eligiS = (spyear`g' == 2010 & fallschoolcode`g' == "396005" & fallgrade`g' == `g') if samp`g'==1
		replace gf_lau_eligiS = (spyear`g' == 2010 & fallschoolcode`g' == "396013" & fallgrade`g' == `g') if samp`g'==1
		replace gf_liv_eligiS = (spyear`g' == 2010 & fallschoolcode`g' == "396014" & fallgrade`g' == `g') if samp`g'==1
		replace gf_har_eligiS = (spyear`g' == 2010 & fallschoolcode`g' == "396030" & fallgrade`g' == `g') if samp`g'==1
		replace gf_gen_eligiS = (spyear`g' == 2010 & fallschoolcode`g' == "396040" & fallgrade`g' == `g') if samp`g'==1
		replace gf_sar_eligiS = (spyear`g' == 2011 & fallschoolcode`g' == "396021" & fallgrade`g' == `g') if samp`g'==1
		replace gf_jos_eligiS = (spyear`g' == 2012 & fallschoolcode`g' == "396001" & fallgrade`g' == `g') if samp`g'==1	
		replace gf_hab1_eligiS= (spyear`g' == 2013 & fallschoolcode`g' == "396034" & fallgrade`g' == `g') if samp`g'==1
		replace gf_hab2_eligiS= (spyear`g' == 2013 & fallschoolcode`g' == "396203" & fallgrade`g' == `g') if samp`g'==1
		replace gf_sch1_eligiS= (spyear`g' == 2013 & fallschoolcode`g' == "396009" & fallgrade`g' == `g') if samp`g'==1
		replace gf_sch2_eligiS= (spyear`g' == 2013 & fallschoolcode`g' == "396010" & fallgrade`g' == `g') if samp`g'==1
		foreach s in dil lau liv har gen sar jos hab1 hab2 sch1 sch2 {
			replace gf_`s'_eligiS = . if spyear`g' == . & samp`g'==1
		}
		
	}
	
	foreach stub in R wB iD iS {
		egen gf_any_elig`stub' = rowmax(gf_*_elig`stub')
	}
	count if gf_any_eligiS != gf_any_elig & matchedsample==1
	if r(N)!=0 {
		stop
	}

	local row = 1
	local col = 1

	forval i = 1/4 {
		egen has_score`i' = rowmax(has_mscore`i' has_escore`i')
	}

	local c = 1
	foreach sch in dil lau liv har gen sar jos hab1 hab2 sch1 sch2 any {
			matrix cohorts[`row',`col'] = `c'
			local ++c
		local ++col
		bys match_cell: egen temp = max(gf_`sch'_eligiS)
		count if gf_`sch'_eligR == 1
			matrix cohorts[`row',`col'] = r(N)
		local ++col
		count if gf_`sch'_eligwB == 1
			matrix cohorts[`row',`col'] = r(N)
		local ++col
		count if gf_`sch'_eligiD == 1
			matrix cohorts[`row',`col'] = r(N)
		local ++col
		count if gf_`sch'_eligiS == 1
			matrix cohorts[`row',`col'] = r(N)
		local ++col
		count if temp == 1 & gf_`sch'_eligiS == 0
			matrix cohorts[`row',`col'] = r(N)
		drop temp
		local col = 1
		local ++row
	}
	clear
	svmat cohorts
	br
	
}
***********************************************
* Check GF balance                            *
***********************************************
if `balance' == 1 {
	use "${cleandata}NOLA_rsdbl_gf2014", clear
	egen temp = rowtotal(samp*)

	foreach var in m e c s lep {
		gen pbl_`var' = .
		gen ppbl_`var' = .
	}
	forval samp = 4/7 {
		local pbl  = `samp' - 2
		local ppbl = `samp' - 3
		foreach p in pbl ppbl {
			cap replace `p'_m = leap_crsdnola_mathsc``p'' if samp`samp' == 1
			cap replace `p'_e = leap_crsdnola_elasc``p'' if samp`samp' == 1	
			cap replace `p'_c = leap_crsdnola_sciencesc``p'' if samp`samp' == 1	
			cap replace `p'_s = leap_crsdnola_socialsc``p'' if samp`samp' == 1	
			cap replace `p'_lep = lep``p'' if samp`samp' == 1
		}
	}

	gen compsamp = 0
	forval g = 4/7 {
		local bl_gr = `g' - 1
		gen compsamp`g' = 0
		qui levelsof spyear`bl_gr' if samp`g'==1 & spyear`bl_gr'!=., local(years)
		qui count if samp`g'==1
		local denom = r(N)
		local prps ""
		foreach y in `years' {
			qui count if samp`g'==1 & spyear`bl_gr'==`y'
			local prop`y' = r(N)/`denom'
			local prps "`prps',`prop`y''"
		}
		local prps = substr("`prps'",2,.)
		disp "`prps'"
		local maxprop = max(`prps')
		disp `maxprop'
		foreach y in `years' {
			replace compsamp`g' = runiform()<(`prop`y''/`maxprop') if spyear`bl_gr'==`y' /*& schoolcode`g'!=""*/ & compsamp==0 & rsd_nola`bl_gr'==1 & alternative`bl_gr'==0
		}
		replace compsamp = 1 if compsamp`g' == 1
	}

	keep if matchedsamp==1 | compsamp==1
	gen any_years = (years_any1>0 & years_any1!=.) | (years_any2>0 & years_any2!=.) | (years_any3>0 & years_any3!=.) if years_any1!=. | years_any2!=. | years_any3!=.
		
	gen alt_sped_bl = .
	gen alt_foodfred_bl = .
	gen alt_lep_bl = .
	foreach t in m e c s {
		gen alt_bl_`t' = .
		gen alt_pbl_`t' = .
	}
	gen alt_charter = .
	forval t = 0/4 {
		gen alt_has_score`t' = .
	}	
	forval g = 4/7 {
		local plg = `g' - 2
		local lg  = `g'-1
		local gp0 = `g'
		local gp1 = `g' + 1
		local gp2 = `g' + 2
		local gp3 = `g' + 3
		local gp4 = `g' + 4
		
		replace alt_sped_bl = sped`lg' if compsamp`g' == 1
		replace alt_foodfred_bl = foodfred`lg' if compsamp`g' == 1
		replace alt_lep_bl = lep`lg' if compsamp`g' == 1
		replace alt_bl_m = leap_crsdnola_mathsc`lg' if compsamp`g' == 1
		replace alt_bl_e = leap_crsdnola_elasc`lg' if compsamp`g' == 1
		replace alt_bl_c = leap_crsdnola_sciencesc`lg' if compsamp`g' == 1
		replace alt_bl_s = leap_crsdnola_socialsc`lg' if compsamp`g' == 1
		replace alt_pbl_m = leap_crsdnola_mathsc`plg' if compsamp`g' == 1
		replace alt_pbl_e = leap_crsdnola_elasc`plg' if compsamp`g' == 1
		replace alt_pbl_c = leap_crsdnola_sciencesc`plg' if compsamp`g' == 1
		replace alt_pbl_s = leap_crsdnola_socialsc`plg' if compsamp`g' == 1
		forval out = 0/4 {
			if `gp`out'' <= 8 {
				replace alt_has_score`out' = (leap_crsdnola_mathsc`gp`out''!=. | leap_crsdnola_elasc`gp`out''!=. | leap_crsdnola_sciencesc`gp`out''!=. | leap_crsdnola_socialsc`gp`out''!=.) if compsamp`g' == 1 & spyear`lg'+1+`out'<=2014
			}
		}
		
	}

	replace alt_charter = (charter4==1 | charter5==1 | charter6==1 | charter7==1 | charter8==1) if (charter4!=. | charter5!=. | charter6!=. | charter7!=. | charter8!=.) & compsamp4==1
	replace alt_charter = (charter5==1 | charter6==1 | charter7==1 | charter8==1) if (charter5!=. | charter6!=. | charter7!=. | charter8!=.) & compsamp5==1
	replace alt_charter = (charter6==1 | charter7==1 | charter8==1) if (charter6!=. | charter7!=. | charter8!=.) & compsamp6==1
	replace alt_charter = (charter7==1 | charter8==1) if (charter7!=. | charter8!=.) & compsamp7==1
	
	local row = 1
	local col = 1

	forval i = 0/4 {
		egen has_score`i' = rowmax(has_mscore`i' has_escore`i')
	}
	
	gen alt_g_m = alt_bl_m - alt_pbl_m
	gen alt_g_e = alt_bl_e - alt_pbl_e
	foreach v in hisp black white asian female alt_sped_bl alt_foodfred_bl alt_lep alt_bl_m alt_g_m alt_bl_e alt_g_e alt_has_score0 alt_has_score1 alt_has_score2 alt_has_score3 alt_has_score4 {
		summ `v' if compsamp==1
			matrix balance[`row',`col']   = r(mean)
			matrix balance[`row'+1,`col'] = r(sd)
			matrix balance[`row'+2,`col'] = r(N)
		local ++col
		summ `v' if compsamp==1 & alt_charter==1
			matrix balance[`row',`col']   = r(mean)
			matrix balance[`row'+1,`col'] = r(sd)
			matrix balance[`row'+2,`col'] = r(N)
		local col = 1
		local row = `row' + 3
	}
	
	local row = 1
	local col = 3
	gen g_m = bl_m - pbl_m
	gen g_e = bl_e - pbl_e
	foreach v in hisp black white asian female bl_sped bl_foodfred bl_lep bl_m g_m bl_e g_e has_score0 has_score1 has_score2 has_score3 has_score4 {
		summ `v' if matchedsample==1 & any_years==1
			matrix balance[`row',`col']   = r(mean)
			matrix balance[`row'+1,`col'] = r(sd)
			matrix balance[`row'+2,`col'] = r(N)	
		local ++col
		summ `v' if matchedsample==1 & gf_any_elig==1
			matrix balance[`row',`col']   = r(mean)
			matrix balance[`row'+1,`col'] = r(sd)
			matrix balance[`row'+2,`col'] = r(N)
		local ++col
		if substr("`v'",1,9) != "has_score" {
			qui areg gf_any_elig if matchedsample==1 & `v'!=. & years_any1!=., absorb(match_cell)
			predict pZ, xbd
			replace pZ = 0.000001 if pZ<=0
			replace pZ = 0.999999 if pZ>=1
			gen kappa = 1 - gf_any_elig*(1-years_any1)/pZ - (1-gf_any_elig)*years_any1/(1-pZ)
			summ `v' if matchedsample==1 & `v'!=. & years_any1!=. [iw=kappa]
				matrix balance[`row',`col']   = r(mean)
				matrix balance[`row'+1,`col'] = r(sd)
			summ kappa if matchedsample==1 & `v'!=. & years_any1!=.
				matrix balance[`row'+2,`col'] = `=`=r(N)'*`=r(mean)''
			drop pZ kappa
		}
		local ++col
		summ `v' if matchedsample==1 & gf_any_elig==0
			matrix balance[`row',`col']   = r(mean)
			matrix balance[`row'+1,`col'] = r(sd)
			matrix balance[`row'+2,`col'] = r(N)			
		local col = 3
		local row = `row' + 3
	}
	
	local row = 1
	local col = 7
	foreach v in hisp black white asian female bl_sped bl_foodfred bl_lep bl_m g_m bl_e g_e {
		if "`v'" == "bl_lep" | "`v'" == "bl_m" | "`v'" == "g_m" | "`v'" == "bl_e" | "`v'" == "g_e" {
			areg `v' gf_any_elig if matchedsample==1, absorb(match_cell) r
				matrix balance[`row',`col']   = _b[gf_any_elig]
				matrix balance[`row'+1,`col'] = _se[gf_any_elig]
				matrix balance[`row'+2,`col'] = e(N)
			areg `v' gf_any_elig if matchedsample==1 & has_score1==1 & has_score0==1, absorb(match_cell) r
				matrix balance[`row',`col'+1]   = _b[gf_any_elig]
				matrix balance[`row'+1,`col'+1] = _se[gf_any_elig]
				matrix balance[`row'+2,`col'+1] = e(N)
		}
		local row = `row' + 3
	}
	
	foreach v in has_score0 has_score1 has_score2 has_score3 has_score4 {
		areg `v' gf_any_elig if matchedsample==1, absorb(match_cell) r
			matrix balance[`row',`col']   = _b[gf_any_elig]
			matrix balance[`row'+1,`col'] = _se[gf_any_elig]
			matrix balance[`row'+2,`col'] = e(N)	
		local row = `row' + 3
	}
	
clear
svmat balance
br
}

***********************************************
* Run basic GF spec (w/ and w/o gains)        *
***********************************************
if `gfing' == 1 {
	use "${cleandata}NOLA_rsdbl_gf2014", clear
	keep if matchedsample == 1
	
	expand 4
	bys sasid: gen copy = _n
	qui tab copy, gen(copy_)
	gen years_any = years_any1 if copy==1
	replace years_any = years_any2 if copy==2
	replace years_any = years_any3 if copy==3
	replace years_any = years_any4 if copy==4
	
	gen grade     = 5 if (samp4==1 & copy==1)
	replace grade = 6 if (samp4==1 & copy==2) | (samp5==1 & copy==1)
	replace grade = 7 if (samp4==1 & copy==3) | (samp5==1 & copy==2) | (samp6==1 & copy==1)
	replace grade = 8 if (samp4==1 & copy==4) | (samp5==1 & copy==3) | (samp6==1 & copy==2) | (samp7==1 & copy==1)
	qui tab grade, gen(outg_)
	
	gen yearout = .
	forval grade = 5/8 {
		replace yearout = spyear`grade' if grade==`grade'
	}
	qui tab yearout, gen(y_)

	local row = 1
	local col = 1
	
	gen copyout = copy
	replace copyout = 34 if copy==3 | copy==4
	replace grade = 78 if grade==7 | grade==8
	replace grade = 56 if grade==5 | grade==6
	
	foreach t in m e c s {
	
		gen `t'score  = `t'score1 if copy==1
		replace `t'score = `t'score2 if copy==2
		replace `t'score = `t'score3 if copy==3
		replace `t'score = `t'score4 if copy==4

		summ `t'score if gf_any_elig==0 & `t'gfscore0!=.
				matrix gfingnoFD[`row',`col']   = r(mean)
				matrix gfingnoFD[`row'+2,`col'] = r(N)
			
			/* OLS */
			reg `t'score years_any  d_* `covs' copy_* outg_* if `t'gfscore0!=., cluster(sasid)
				matrix gfingnoFD[`row',`col'+1]   = _b[years_any]
				matrix gfingnoFD[`row'+1,`col'+1] = _se[years_any]
				matrix gfingnoFD[`row'+2,`col'+1] = e(N)	
				
			/* RF */
			areg `t'score gf_any_elig `covs' copy_* outg_* if years_any!=. & `t'gfscore0!=., absorb(match_cell) cluster(sasid)
				matrix gfingnoFD[`row',`col'+2]   = _b[gf_any_elig]
				matrix gfingnoFD[`row'+1,`col'+2] = _se[gf_any_elig]
				matrix gfingnoFD[`row'+2,`col'+2] = e(N)
			
			/* FS */
			areg years_any gf_any_elig `covs' copy_* outg_* if `t'score!=. & `t'gfscore0!=., absorb(match_cell) cluster(sasid)
				matrix gfingnoFD[`row',`col'+3]   = _b[gf_any_elig]
				matrix gfingnoFD[`row'+1,`col'+3] = _se[gf_any_elig]
				matrix gfingnoFD[`row'+2,`col'+3] = e(N)
				
			/* 2SLS */
			ivregress 2sls `t'score (years_any = gf_any_elig) d_* `covs' copy_* outg_* if `t'gfscore0!=., cluster(sasid)
				matrix gfingnoFD[`row',`col'+4]   = _b[years_any]
				matrix gfingnoFD[`row'+1,`col'+4] = _se[years_any]
				matrix gfingnoFD[`row'+2,`col'+4] = e(N)
			
		local row = `row' + 3		
		
	}
	
	foreach copy in 1 2 34 {
		foreach t in m e c s {

			summ `t'score if gf_any_elig==0 & `t'gfscore0!=. & copyout==`copy'
				matrix gfingnoFD[`row',`col']   = r(mean)
				matrix gfingnoFD[`row'+2,`col'] = r(N)
				
			/* OLS */
			reg `t'score years_any d_* `covs' copy_* outg_* if  `t'gfscore0!=. & copyout==`copy', cluster(sasid)
				matrix gfingnoFD[`row',`col'+1]   = _b[years_any]
				matrix gfingnoFD[`row'+1,`col'+1] = _se[years_any]
				matrix gfingnoFD[`row'+2,`col'+1] = e(N)
				
				
			/* RF */
			areg `t'score gf_any_elig `covs' copy_* outg_* if years_any!=. & `t'gfscore0!=. & copyout==`copy', absorb(match_cell) cluster(sasid)
				matrix gfingnoFD[`row',`col'+2]   = _b[gf_any_elig]
				matrix gfingnoFD[`row'+1,`col'+2] = _se[gf_any_elig]
				matrix gfingnoFD[`row'+2,`col'+2] = e(N)
			
			/* FS */
			areg years_any gf_any_elig `covs' copy_* outg_* if `t'score!=. & `t'gfscore0!=. & copyout==`copy', absorb(match_cell) cluster(sasid)
				matrix gfingnoFD[`row',`col'+3]   = _b[gf_any_elig]
				matrix gfingnoFD[`row'+1,`col'+3] = _se[gf_any_elig]
				matrix gfingnoFD[`row'+2,`col'+3] = e(N)
				
			/* 2SLS */
			ivregress 2sls `t'score (years_any = gf_any_elig) d_* `covs' copy_* outg_* if  `t'gfscore0!=. & copyout==`copy', cluster(sasid)
				matrix gfingnoFD[`row',`col'+4]   = _b[years_any]
				matrix gfingnoFD[`row'+1,`col'+4] = _se[years_any]
				matrix gfingnoFD[`row'+2,`col'+4] = e(N)
			

		local row = `row' + 3
		}
	}	

	foreach grade in 56 78 {
		foreach t in m e c s{

			summ `t'score if gf_any_elig==0 & `t'gfscore0!=. & grade==`grade'
				matrix gfingnoFD[`row',`col']   = r(mean)
				matrix gfingnoFD[`row'+2,`col'] = r(N)
				
			/* OLS */
			reg `t'score years_any d_* `covs' copy_* outg_* if  `t'gfscore0!=. & grade==`grade', cluster(sasid)
				matrix gfingnoFD[`row',`col'+1]   = _b[years_any]
				matrix gfingnoFD[`row'+1,`col'+1] = _se[years_any]
				matrix gfingnoFD[`row'+2,`col'+1] = e(N)		
				
			/* RF */
			areg `t'score gf_any_elig `covs' copy_* outg_* if years_any!=. & `t'gfscore0!=. & grade==`grade', absorb(match_cell) cluster(sasid)
				matrix gfingnoFD[`row',`col'+2]   = _b[gf_any_elig]
				matrix gfingnoFD[`row'+1,`col'+2] = _se[gf_any_elig]
				matrix gfingnoFD[`row'+2,`col'+2] = e(N)
			
			/* FS */
			areg years_any gf_any_elig `covs' copy_* outg_* if `t'score!=. & `t'gfscore0!=. & grade==`grade', absorb(match_cell) cluster(sasid)
				matrix gfingnoFD[`row',`col'+3]   = _b[gf_any_elig]
				matrix gfingnoFD[`row'+1,`col'+3] = _se[gf_any_elig]
				matrix gfingnoFD[`row'+2,`col'+3] = e(N)
				
			/* 2SLS */
			ivregress 2sls `t'score (years_any = gf_any_elig) d_* `covs' copy_* outg_* if  `t'gfscore0!=. & grade==`grade', cluster(sasid)
				matrix gfingnoFD[`row',`col'+4]   = _b[years_any]
				matrix gfingnoFD[`row'+1,`col'+4] = _se[years_any]
				matrix gfingnoFD[`row'+2,`col'+4] = e(N)
			
		local row = `row' + 3
		}

	}

	local row = 1
	/* linear exclusion restriction fix */	
	foreach t in m e c  s {
	
			gen `t'level = .
			replace `t'level= `t'score
			
			replace `t'score = `t'score1-`t'gfscore0 if copy==1
			replace `t'score = `t'score2-`t'gfscore0 if copy==2
			replace `t'score = `t'score3-`t'gfscore0 if copy==3
			replace `t'score = `t'score4-`t'gfscore0 if copy==4
			
			*Show level for comparison group mean
			summ `t'level if gf_any_elig==0 & `t'gfscore0!=.
				matrix gfing[`row',`col']   = r(mean)
				matrix gfing[`row'+2,`col'] = r(N)
				
			/* OLS */
			reg `t'score years_any d_* `covs' copy_* outg_*, cluster(sasid)
				matrix gfing[`row',`col'+1]   = _b[years_any]
				matrix gfing[`row'+1,`col'+1] = _se[years_any]
				matrix gfing[`row'+2,`col'+1] = e(N)		
				
			/* RF */
			areg `t'score gf_any_elig `covs' copy_* outg_* if years_any!=., absorb(match_cell) cluster(sasid)
				matrix gfing[`row',`col'+2]   = _b[gf_any_elig]
				matrix gfing[`row'+1,`col'+2] = _se[gf_any_elig]
				matrix gfing[`row'+2,`col'+2] = e(N)
			
			/* FS */
			areg years_any gf_any_elig `covs' copy_* outg_* if `t'score!=., absorb(match_cell) cluster(sasid)
				matrix gfing[`row',`col'+3]   = _b[gf_any_elig]
				matrix gfing[`row'+1,`col'+3] = _se[gf_any_elig]
				matrix gfing[`row'+2,`col'+3] = e(N)
				
			/* 2SLS */
			ivregress 2sls `t'score (years_any = gf_any_elig) d_* `covs' copy_* outg_*, cluster(sasid)
				matrix gfing[`row',`col'+4]   = _b[years_any]
				matrix gfing[`row'+1,`col'+4] = _se[years_any]
				matrix gfing[`row'+2,`col'+4] = e(N)	
			
		local row = `row' + 3
		local col = 1
	}

	foreach copy in 1 2 34 {
		foreach t in m e c s {

			*Show level for comparison group mean
			summ `t'level if gf_any_elig==0 & `t'gfscore0!=. & copyout==`copy'
				matrix gfing[`row',`col']   = r(mean)
				matrix gfing[`row'+2,`col'] = r(N)
			
			/* OLS */
			reg `t'score years_any d_* `covs' copy_* outg_* if copyout==`copy', cluster(sasid)
				matrix gfing[`row',`col'+1]   = _b[years_any]
				matrix gfing[`row'+1,`col'+1] = _se[years_any]
				matrix gfing[`row'+2,`col'+1] = e(N)	
			
			/* RF */
			areg `t'score gf_any_elig `covs' copy_* outg_* if years_any!=. & copyout==`copy', absorb(match_cell) cluster(sasid)
				matrix gfing[`row',`col'+2]   = _b[gf_any_elig]
				matrix gfing[`row'+1,`col'+2] = _se[gf_any_elig]
				matrix gfing[`row'+2,`col'+2] = e(N)
			
			/* FS */
			areg years_any gf_any_elig `covs' copy_* outg_* if `t'score!=. & copyout==`copy', absorb(match_cell) cluster(sasid)
				matrix gfing[`row',`col'+3]   = _b[gf_any_elig]
				matrix gfing[`row'+1,`col'+3] = _se[gf_any_elig]
				matrix gfing[`row'+2,`col'+3] = e(N)
				
			/* 2SLS */
			ivregress 2sls `t'score (years_any = gf_any_elig) d_* `covs' copy_* outg_* if copyout==`copy', cluster(sasid)
				matrix gfing[`row',`col'+4]   = _b[years_any]
				matrix gfing[`row'+1,`col'+4] = _se[years_any]
				matrix gfing[`row'+2,`col'+4] = e(N)
			
		local row = `row' + 3
		local col = 1
		}
	}

	foreach grade in 56 78 {
		foreach t in m e c s {

			*Show level for comparison group mean
			summ `t'level if gf_any_elig==0 & `t'gfscore0!=. & grade==`grade'
				matrix gfing[`row',`col']   = r(mean)
				matrix gfing[`row'+2,`col'] = r(N)
			
			/* OLS */
			reg `t'score years_any d_* `covs' copy_* outg_* if grade==`grade', cluster(sasid)
				matrix gfing[`row',`col'+1]   = _b[years_any]
				matrix gfing[`row'+1,`col'+1] = _se[years_any]
				matrix gfing[`row'+2,`col'+1] = e(N)	
			
			
			/* RF */
			areg `t'score gf_any_elig `covs' copy_* outg_* if years_any!=. & grade==`grade', absorb(match_cell) cluster(sasid)
				matrix gfing[`row',`col'+2]   = _b[gf_any_elig]
				matrix gfing[`row'+1,`col'+2] = _se[gf_any_elig]
				matrix gfing[`row'+2,`col'+2] = e(N)
			
			/* FS */
			areg years_any gf_any_elig `covs' copy_* outg_* if `t'score!=. & grade==`grade', absorb(match_cell) cluster(sasid)
				matrix gfing[`row',`col'+3]   = _b[gf_any_elig]
				matrix gfing[`row'+1,`col'+3] = _se[gf_any_elig]
				matrix gfing[`row'+2,`col'+3] = e(N)
				
			/* 2SLS */
			ivregress 2sls `t'score (years_any = gf_any_elig) d_* `covs' copy_* outg_* if grade==`grade', cluster(sasid)
				matrix gfing[`row',`col'+4]   = _b[years_any]
				matrix gfing[`row'+1,`col'+4] = _se[years_any]
				matrix gfing[`row'+2,`col'+4] = e(N)
			
		local row = `row' + 3
		local col = 1
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
	use "${cleandata}NOLA_rsdbl_gf2014", clear
	keep if matchedsample == 1
	
	local row = 1
	local col = 1
		
	forval out = 1/4 {
	
		gen year_`out' = .
		gen grade_`out' = .
		gen ms_samp_`out' = 0
		gen takeover_full_`out' = .
		gen takeover_other_`out' = .
		gen startup_`out' = .
		gen unclear_`out' = .
		forval gr = 4/7 {
			local grout = `gr' + `out'
			local lg = `gr' - 1
			if `grout' <= 8 {
				replace year_`out' = spyear`grout' if samp`gr'==1
				replace grade_`out' = `grout' if samp`gr'==1
				gen schoolcode = schoolcode`grout' if samp`gr'==1 & spyear`lg'+1+`out'<=2014 & ((has_mscore`out'==1 & mgfscore0!=.) | (has_escore`out'==1 & egfscore0!=.)) & bl_m!=. & bl_e!=. & bl_c!=. & bl_s!=.
				gen schoolname = schoolname`grout' if samp`gr'==1 & spyear`lg'+1+`out'<=2014 & ((has_mscore`out'==1 & mgfscore0!=.) | (has_escore`out'==1 & egfscore0!=.)) & bl_m!=. & bl_e!=. & bl_c!=. & bl_s!=.
				merge m:1 schoolcode schoolname using "${types}rsd_charter_types", keepusing(takeover_full takeover_other startup unclear) keep(master match)
				replace takeover_full_`out' = takeover_full if _merge==3
				replace takeover_other_`out' = takeover_other if _merge==3
				replace startup_`out' = startup if _merge==3
				replace unclear_`out' = unclear if _merge==3
				drop schoolcode schoolname takeover_full takeover_other startup unclear _merge
			}
			replace ms_samp_`out' = (`grout'<=8) if samp`gr'==1 & spyear`lg'+1+`out'<=2014 & ((has_mscore`out'==1 & mgfscore0!=.) | (has_escore`out'==1 & egfscore0!=.)) & bl_m!=. & bl_e!=. & bl_c!=. & bl_s!=.
		}

		qui tab year_`out', gen(y_)
		qui tab grade_`out', gen(gout_)
		
		gen study_`out'      = incharter`out'==1 & attend_any`out'==1 & (has_mscore`out'==1 | has_escore`out'==1) if ms_samp_`out'==1
		gen othtake_`out'    = incharter`out'==1 & attend_any`out'==0 & (takeover_full_`out'==1 | takeover_other_`out'==1) & (has_mscore`out'==1 | has_escore`out'==1) if ms_samp_`out'==1
		gen othch_`out'      = incharter`out'==1 & attend_any`out'==0 & (startup_`out'==1 | unclear_`out'==1) & (has_mscore`out'==1 | has_escore`out'==1) if ms_samp_`out'==1
		gen public_`out'     = incharter`out'==0 & attend_any`out'!=. & (has_mscore`out'==1 | has_escore`out'==1) if ms_samp_`out'==1
		gen leftrsd_`out'    = attend_any`out' == . | (has_mscore`out'==0 & has_escore`out'==0) if ms_samp_`out'==1
		
		gen D_`out' = study_`out' if ms_samp_`out'==1
		gen oneminusD_`out' = 1-D_`out'
		gen Z = gf_any_elig
		
		gen samp = ms_samp_`out' 
		foreach var of varlist gf_any_elig D_`out' study_`out' othtake_`out' othch_`out' public_`out' leftrsd_`out' d_* gout_* `covs' {
			replace samp = 0 if `var' == .
		}
		
		qui probit Z d_* gout_* `covs'
		predict pZ if samp == 1
		gen kappa_0 = (1-D_`out')*((1-Z)-(1-pZ))/(pZ*(1-pZ)) if samp==1

		foreach s in study othtake othch public leftrsd {
			summ `s'_`out' if gf_any_elig==0 & ms_samp_`out'==1
				matrix choice[`row',1] = r(mean)
				matrix choice[`row'+1,1] = r(N)
			summ `s'_`out' if gf_any_elig==1 & ms_samp_`out'==1
				matrix choice[`row',2] = r(mean)	
				matrix choice[`row'+1,2] = r(N)

			gen Y = `s'_`out'*oneminusD_`out'

			cap noisily {
				summ `s'_`out' [iweight=kappa_0]
				matrix choice[`row',3] = r(mean)
				
			}
			drop Y
			local ++row
			if "`s'" == "leftrsd" {
				local ++row
			}
		}		

		local ++row
		drop y_* gout_* samp pZ kappa_0 Z
	}
	clear
	svmat choice
	br

}

***********************************************
* RSD growth in LEAP/iLEAP status             *
***********************************************
if `statusgrowth' == 1 {
	use "${cleandata}SIMS_LEAP_20082014_wide `simsdate'", clear

	forval g = 5/8 {
		egen b_or_a_m`g' = rowmax(leap_math_achieve1`g' leap_math_achieve2`g' leap_math_achieve3`g') if rsd_nola`g'==1 & alternative`g'==0
		egen b_or_a_e`g' = rowmax(leap_ela_achieve1`g' leap_ela_achieve2`g' leap_ela_achieve3`g') if rsd_nola`g'==1 & alternative`g'==0		
	}
	keep b_or_a_m* b_or_a_e* spyear5 spyear6 spyear7 spyear8 charter5 charter6 charter7 charter8 sasid
	
	reshape long spyear b_or_a_m b_or_a_e charter, i(sasid) j(grade)
	collapse (mean) b_or_a_m b_or_a_e, by(spyear charter grade)
	
	collapse (mean) b_or_a_m b_or_a_e, by(spyear charter)
	
	preserve
		u "${cleandata}LEAPiLEAP_assessment_LA", clear
		gen b_or_a_m = (math_a + math_m + ela_b)/100
		gen b_or_a_e = (ela_a + ela_m + ela_b)/100
		keep if inlist(leacode,"036","STATE")
		keep b_or_a_* spyear grade leacode
		collapse (mean) b_or_a_*, by(spyear leacode)
		ren leacode group
		replace group="OPSB" if group=="036"
		tempfile state
		sa `state'
	restore
	
	gen group = "RSD charter" if charter==1
	replace group = "RSD direct" if charter==0
	drop if charter==.
	drop charter
	
	append using `state'
	
	reshape long b_or_a_, i(spyear group) j(subject) string
	ren b_or_a_ pct
	foreach t in m e{
		preserve
			keep if subject=="`t'"
			twoway 	(line pct spyear if group=="RSD charter", lcolor(gs8) lwidth(medthick)) ///
						(line pct spyear if group=="RSD direct", lcolor(black) lwidth(medthick)) ///
						(line pct spyear if group=="OPSB", lcolor(gs8) lpattern(dash)) ///
						(line pct spyear if group=="STATE", lcolor(black) lpattern(dash)), ///
						xtitle("Spring year",size(small)) ///
						ytitle("Percent at basic or above",size(small)) ///
						scheme(s1mono) ylab(0.2(0.1)1.0, labsize(small) nogrid) xlab(2008(2)2014,labsize(small)) ///
						legend(	label(1 "RSD charter") label(2 "RSD direct-run") label(3 "OPSB") label(4 "Louisiana")   ///
								cols(2) size(small) )
												
			graph save "${graphs}RSD_LA_achievement_`t'", replace
			graph export "${graphs}RSD_LA_achievement_`t'.wmf", replace
		restore
	}
	
}
***********************************************
* Estimate charter fallback                   *
***********************************************
if `fallbackchange' == 1 {
	use "${cleandata}NOLA_rsdbl_gf2014", clear
	keep if matchedsample == 1

	expand 4
	bys sasid: gen copy = _n
	qui tab copy, gen(copy_)
	gen years_any = years_any1 if copy==1
	replace years_any = years_any2 if copy==2
	replace years_any = years_any3 if copy==3
	replace years_any = years_any4 if copy==4
	
	gen grade     = 5 if (samp4==1 & copy==1)
	replace grade = 6 if (samp4==1 & copy==2) | (samp5==1 & copy==1)
	replace grade = 7 if (samp4==1 & copy==3) | (samp5==1 & copy==2) | (samp6==1 & copy==1)
	replace grade = 8 if (samp4==1 & copy==4) | (samp5==1 & copy==3) | (samp6==1 & copy==2) | (samp7==1 & copy==1)
	
	qui tab grade, gen(outg_)
	
	gen yearout = .
	forval grade = 5/8 {
		replace yearout = spyear`grade' if grade==`grade'
	}
	qui tab yearout, gen(y_)
	
	forval out = 1/4 {
		replace incharter`out' = incharter`out'*(1-attend_any`out')
	}
	gen years_charter1 = incharter1
	gen temp = incharter1*(1+repeatout1)
	egen years_charter2 = rowtotal(temp incharter2)
	replace years_charter2 = . if incharter1==. & incharter2==.
	gen temp2 = incharter2*(1+repeatout2)
	egen years_charter3 = rowtotal(temp temp2 incharter3)
	replace years_charter3 = . if incharter1==. & incharter2==. & incharter3==.	
	gen temp3 = incharter3*(1+repeatout3)
	egen years_charter4 = rowtotal(temp temp2 temp3 incharter4)
	replace years_charter4 = . if incharter1==. & incharter2==. & incharter3==. & incharter4==.			
	drop temp temp2 temp3

	gen years_charter = years_charter1 if copy==1
	replace years_charter = years_charter2 if copy==2
	replace years_charter = years_charter3 if copy==3
	replace years_charter = years_charter4 if copy==4
		
	gen years_anycharter = years_any + years_charter
	
	gen baselineyr = .
	gen baselinespsbin = .
	gen baselinesps = .
	forval gr = 4/7 {
		local lgr = `gr' - 1
		replace baselineyr = spyear`lgr' if samp`gr'==1 
		replace baselinespsbin = int(sps`lgr'/5) if samp`gr'==1
		replace baselinesps = sps`lgr' if samp`gr'==1
		
	}

	egen temp = group(bl_sped baselinespsbin) 
	qui tab temp, gen(sd1_)
	drop temp
	egen temp = group(samp4 samp5 samp6 samp7)
	qui tab temp, gen(sd2_)
	drop temp
	qui tab bl_year, gen(sd3_)
	foreach var1 of varlist sd*_* {
		foreach var2 of varlist gf_any_elig {
			gen `var2'_`var1' = `var2'*`var1'
		}
	}
	gen temp = (bl_year > 2009)
	qui tab temp, gen(bd1_)
	foreach var1 of varlist bd*_* {
		foreach var2 of varlist gf_any_elig {
			gen `var2'_`var1' = `var2'*`var1'
		}
	}	

	local row = 1
	local col = 1
	foreach t in m e {

		gen `t'score  = `t'score1 - `t'gfscore0 if copy==1
		replace `t'score  = `t'score2 - `t'gfscore0 if copy==2
		replace `t'score  = `t'score3 - `t'gfscore0 if copy==3
		replace `t'score  = `t'score4 - `t'gfscore0 if copy==4
	
		ivreg2 `t'score (years_any = gf_*_elig_sd*_*) sd*_* d_* outg_* `covs' copy_* if `t'gfscore0!=., cluster(sasid) partial(d_* `covs' copy_* sd*_* outg_*) savefirst first
			matrix temp = e(first)
			matrix fallbackchange[`row',`col']    = _b[years_any]
			matrix fallbackchange[`row'+1,`col']  = _se[years_any]
			matrix fallbackchange[`row'+2,`col']  = temp[7,1]
			matrix fallbackchange[`row'+9,`col']  = e(exexog_ct)
			matrix fallbackchange[`row'+10,`col'] = e(jp)
			matrix fallbackchange[`row'+11,`col'] = e(N)
	
		local ++col
		ivreg2 `t'score (years_any years_charter = gf_*_elig_sd*_*) sd*_* d_* outg_* `covs' copy_* if `t'gfscore0!=., cluster(sasid) partial(d_* `covs' copy_* sd*_* outg_*) savefirst first
		matrix temp = e(first)
			matrix fallbackchange[`row',`col']    = _b[years_any]
			matrix fallbackchange[`row'+1,`col']  = _se[years_any]
			matrix fallbackchange[`row'+2,`col']  = temp[7,1]
			matrix fallbackchange[`row'+3,`col']  = _b[years_charter]
			matrix fallbackchange[`row'+4,`col']  = _se[years_charter]
			matrix fallbackchange[`row'+5,`col']  = temp[7,2]
			matrix fallbackchange[`row'+9,`col']  = e(exexog_ct)
			matrix fallbackchange[`row'+10,`col'] = e(jp)
			matrix fallbackchange[`row'+11,`col'] = e(N)

		local ++col
		ivreg2 `t'score (years_any years_charter = gf_*_elig_bd*_*) bd*_* d_* outg_* `covs' copy_* if `t'gfscore0!=. & copy==1, cluster(sasid) partial(d_* `covs' copy_* bd*_* outg_*) savefirst first
		matrix temp = e(first)
			matrix fallbackchange[`row',`col']    = _b[years_any]
			matrix fallbackchange[`row'+1,`col']  = _se[years_any]
			matrix fallbackchange[`row'+2,`col']  = temp[7,1]
			matrix fallbackchange[`row'+3,`col']  = _b[years_charter]
			matrix fallbackchange[`row'+4,`col']  = _se[years_charter]
			matrix fallbackchange[`row'+5,`col']  = temp[7,2]
			matrix fallbackchange[`row'+9,`col']  = e(exexog_ct)
			matrix fallbackchange[`row'+10,`col'] = e(jp)
			matrix fallbackchange[`row'+11,`col'] = e(N)
			
		local ++col
		ivreg2 `t'score (years_anycharter = gf_*_elig_sd*_*) sd*_* d_* `covs' copy_* outg_* if `t'gfscore0!=., cluster(sasid) partial(d_* `covs' copy_* outg_* sd*_*) savefirst first
			matrix temp = e(first)
			matrix fallbackchange[`row'+6,`col']  = _b[years_anycharter]
			matrix fallbackchange[`row'+7,`col']  = _se[years_anycharter]	
			matrix fallbackchange[`row'+8,`col']  = temp[7,1]
			matrix fallbackchange[`row'+9,`col']  = e(exexog_ct)
			matrix fallbackchange[`row'+10,`col'] = e(jp)
			matrix fallbackchange[`row'+11,`col'] = e(N)		
			
		local ++col
		ivreg2 `t'score (years_anycharter = gf_*_elig_bd*_*) bd*_* d_* `covs' copy_* outg_* if `t'gfscore0!=. & copy==1, cluster(sasid) partial(d_* `covs' copy_* outg_* bd*_*) savefirst first
			matrix temp = e(first)
			matrix fallbackchange[`row'+6,`col']  = _b[years_anycharter]
			matrix fallbackchange[`row'+7,`col']  = _se[years_anycharter]	
			matrix fallbackchange[`row'+8,`col']  = temp[7,1]
			matrix fallbackchange[`row'+9,`col']  = e(exexog_ct)
			matrix fallbackchange[`row'+10,`col'] = e(jp)
			matrix fallbackchange[`row'+11,`col'] = e(N)	
			
		local ++col
		ivreg2 `t'score (years_anycharter = gf_any_elig) d_* `covs' copy_* outg_* if `t'gfscore0!=., cluster(sasid) partial(d_* `covs' copy_* outg_*) savefirst first
			matrix temp = e(first)
			matrix fallbackchange[`row'+6,`col']  = _b[years_anycharter]
			matrix fallbackchange[`row'+7,`col']  = _se[years_anycharter]	
			matrix fallbackchange[`row'+8,`col']  = temp[7,1]
			matrix fallbackchange[`row'+9,`col']  = e(exexog_ct)
			matrix fallbackchange[`row'+10,`col'] = e(jp)
			matrix fallbackchange[`row'+11,`col'] = e(N)
			
		local row = `row' + 12
		local col = 1
	}
	
clear
svmat fallbackchange
br
}
***********************************************
* Peer quality                                *
***********************************************
if `switching' == 1 {
	use "${cleandata}NOLA_rsdbl_gf2014_wpeers", clear
	keep if matchedsample == 1
	
	gen grade     = 5 if samp4==1
	replace grade = 6 if samp5==1
	replace grade = 7 if samp6==1
	replace grade = 8 if samp7==1
	qui tab grade, gen(outg_)
	gen yearout = .
	forval grade = 5/8 {
		replace yearout = spyear`grade' if grade==`grade'
	}
	qui tab yearout, gen(y_)
	
	
	local row = 1
	local col = 1	
	
	foreach t in sped foodfred lep mscore escore {
		summ peer_`t'1 if gf_any_elig==0
			matrix switching[`row',`col']   = r(mean)
			matrix switching[`row'+2,`col'] = r(N)

		ivreg2 peer_`t'1 (attend_any1 = gf_any_elig) d_* `covs' , r partial(d_* `covs')
			matrix switching[`row',`col'+1]   = _b[attend_any1]
			matrix switching[`row'+1,`col'+1] = _se[attend_any1]
			matrix switching[`row'+2,`col'+1] = e(N)
		local row = `row' + 3
	
	}

	clear
	svmat switching
	br
}

***********************************************
* Estimate Score Persistence                  *
***********************************************
if `persistence' == 1 {
	use "${cleandata}NOLA_rsdbl_gf2014", clear
	keep if matchedsample == 1

	expand 4
	bys sasid: gen copy = _n
	qui tab copy, gen(copy_)
	gen years_any = years_any1 if copy==1
	replace years_any = years_any2 if copy==2
	replace years_any = years_any3 if copy==3
	replace years_any = years_any4 if copy==4
	
	gen grade     = 5 if (samp4==1 & copy==1)
	replace grade = 6 if (samp4==1 & copy==2) | (samp5==1 & copy==1)
	replace grade = 7 if (samp4==1 & copy==3) | (samp5==1 & copy==2) | (samp6==1 & copy==1)
	replace grade = 8 if (samp4==1 & copy==4) | (samp5==1 & copy==3) | (samp6==1 & copy==2) | (samp7==1 & copy==1)
	
	qui tab grade, gen(outg_)
	
	gen yearout = .
	forval grade = 5/8 {
		replace yearout = spyear`grade' if grade==`grade'
	}
	qui tab yearout, gen(y_)
	
	forval out = 1/4 {
		replace incharter`out' = incharter`out'*(1-attend_any`out')
	}
	gen years_charter1 = incharter1
	gen temp = incharter1*(1+repeatout1)
	egen years_charter2 = rowtotal(temp incharter2)
	replace years_charter2 = . if incharter1==. & incharter2==.
	gen temp2 = incharter2*(1+repeatout2)
	egen years_charter3 = rowtotal(temp temp2 incharter3)
	replace years_charter3 = . if incharter1==. & incharter2==. & incharter3==.	
	gen temp3 = incharter3*(1+repeatout3)
	egen years_charter4 = rowtotal(temp temp2 temp3 incharter4)
	replace years_charter4 = . if incharter1==. & incharter2==. & incharter3==. & incharter4==.			
	drop temp temp2 temp3

	gen years_charter = years_charter1 if copy==1
	replace years_charter = years_charter2 if copy==2
	replace years_charter = years_charter3 if copy==3
	replace years_charter = years_charter4 if copy==4
		
	gen years_anycharter = years_any + years_charter
	
	gen baselineyr = .
	gen baselinespsbin = .
	gen baselinesps = .
	forval gr = 4/7 {
		local lgr = `gr' - 1
		replace baselineyr = spyear`lgr' if samp`gr'==1 
		replace baselinespsbin = int(sps`lgr'/5) if samp`gr'==1
		replace baselinesps = sps`lgr' if samp`gr'==1
		
	}

	egen temp = group(bl_sped baselinespsbin) 
	qui tab temp, gen(sd1_)
	drop temp
	egen temp = group(samp4 samp5 samp6 samp7)
	qui tab temp, gen(sd2_)
	drop temp
	qui tab bl_year, gen(sd3_)
	foreach var1 of varlist sd*_* {
		foreach var2 of varlist gf_any_elig {
			gen `var2'_`var1' = `var2'*`var1'
		}
	}	
	
	local row = 1
	local col = 1
	foreach t in m e {

		gen `t'score  = `t'score1 if copy==1
		replace `t'score  = `t'score2 if copy==2
		replace `t'score  = `t'score3 if copy==3
		replace `t'score  = `t'score4 if copy==4
		gen `t'gain = `t'score - `t'gfscore0

		ivreg2 `t'gfscore0 gf_any_elig d_* outg_* `covs' copy_* if `t'score!=. & copy==1, cluster(sasid) partial(d_* `covs' copy_* outg_*)
			matrix persistence[`row',`col']    = _b[gf_any_elig]
			matrix persistence[`row'+1,`col']  = _se[gf_any_elig]
			matrix persistence[`row'+11,`col'] = e(N)
		
		local ++col
		ivreg2 `t'score (years_any `t'gfscore0 = gf_*_elig_sd*_*) sd*_* d_* outg_* `covs' copy_* if `t'gfscore0!=., cluster(sasid) partial(d_* `covs' copy_* sd*_* outg_*) savefirst first
			matrix temp = e(first)
			matrix persistence[`row'+2,`col']  = _b[years_any]
			matrix persistence[`row'+3,`col']  = _se[years_any]
			matrix persistence[`row'+4,`col']  = temp[7,1]
			matrix persistence[`row'+5,`col']  = _b[`t'gfscore0]
			matrix persistence[`row'+6,`col']  = _se[`t'gfscore0]
			matrix persistence[`row'+7,`col']  = temp[7,2]
			matrix persistence[`row'+10,`col'] = e(exexog_ct)
			matrix persistence[`row'+11,`col'] = e(N)
			
		local ++col		
		ivreg2 `t'gain (years_any = gf_*_elig_sd*_*) sd*_* d_* outg_* `covs' copy_* if `t'gfscore0!=., cluster(sasid) partial(d_* `covs' copy_* sd*_* outg_*) savefirst first
			matrix temp = e(first)
			matrix persistence[`row'+2,`col']  = _b[years_any]
			matrix persistence[`row'+3,`col']  = _se[years_any]
			matrix persistence[`row'+4,`col']  = temp[7,1]
			matrix persistence[`row'+10,`col'] = e(exexog_ct)
			matrix persistence[`row'+11,`col'] = e(N)		

		preserve
		keep if `t'gfscore0!=.
		gen one = 1
		expandcl 2, cluster(one) gen(est)
		foreach var of varlist years_any `t'gfscore0 gf_*_elig_sd*_* sd*_* d_* bl_sped bl_foodfred bl_lep bl_m bl_e bl_c bl_s y_* copy_* outg_* {
			forval est = 1/2 {
				gen est`est'_`var' = `var'*(est==`est')
			}
			drop `var'
		}
		gen Y = `t'score if est==1
		replace Y = `t'gain if est==2
		ivregress 2sls Y (est*_years_any est1_`t'gfscore0 = est*_gf_*_elig_sd*_*) est1_sd*_* est2_sd*_* est*_d_* est*_bl_foodfred est*_bl_lep est*_bl_m est*_bl_e est*_bl_c est*_bl_s est*_y_* est*_copy_* est*_outg_*
		lincom est1_years_any-est2_years_any
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
* Visual IV                                   *
***********************************************
if `fallbackgraph' == 1 {
	use "${cleandata}NOLA_rsdbl_gf2014", clear
	keep if matchedsample == 1

	expand 4
	bys sasid: gen copy = _n
	qui tab copy, gen(copy_)
	gen years_any = years_any1 if copy==1
	replace years_any = years_any2 if copy==2
	replace years_any = years_any3 if copy==3
	replace years_any = years_any4 if copy==4
	
	gen grade     = 5 if (samp4==1 & copy==1)
	replace grade = 6 if (samp4==1 & copy==2) | (samp5==1 & copy==1)
	replace grade = 7 if (samp4==1 & copy==3) | (samp5==1 & copy==2) | (samp6==1 & copy==1)
	replace grade = 8 if (samp4==1 & copy==4) | (samp5==1 & copy==3) | (samp6==1 & copy==2) | (samp7==1 & copy==1)
	
	qui tab grade, gen(outg_)
	
	gen yearout = .
	forval grade = 5/8 {
		replace yearout = spyear`grade' if grade==`grade'
	}
	qui tab yearout, gen(y_)
	
	forval out = 1/4 {
		replace incharter`out' = incharter`out'*(1-attend_any`out')
	}
	gen years_charter1 = incharter1
	gen temp = incharter1*(1+repeatout1)
	egen years_charter2 = rowtotal(temp incharter2)
	replace years_charter2 = . if incharter1==. & incharter2==.
	gen temp2 = incharter2*(1+repeatout2)
	egen years_charter3 = rowtotal(temp temp2 incharter3)
	replace years_charter3 = . if incharter1==. & incharter2==. & incharter3==.	
	gen temp3 = incharter3*(1+repeatout3)
	egen years_charter4 = rowtotal(temp temp2 temp3 incharter4)
	replace years_charter4 = . if incharter1==. & incharter2==. & incharter3==. & incharter4==.			
	drop temp temp2 temp3

	gen years_charter = years_charter1 if copy==1
	replace years_charter = years_charter2 if copy==2
	replace years_charter = years_charter3 if copy==3
	replace years_charter = years_charter4 if copy==4
		
	gen years_anycharter = years_any + years_charter
	
	gen baselineyr = .
	gen baselinespsbin = .
	gen baselinesps = .
	forval gr = 4/7 {
		local lgr = `gr' - 1
		replace baselineyr = spyear`lgr' if samp`gr'==1 
		replace baselinespsbin = int(sps`lgr'/5) if samp`gr'==1
		replace baselinesps = sps`lgr' if samp`gr'==1
		
	}

	egen temp = group(baselineyr bl_sped baselinespsbin)
	qui tab temp, gen(sd1_)
	drop sd1_13 /* omit one category */
	foreach var1 of varlist sd*_* {
		foreach var2 of varlist gf_any_elig {
			gen `var2'_`var1' = `var2'*`var1'
		}
	}	
	
	preserve
	local t = "m"
	gen `t'score  = `t'score1 - `t'gfscore0 if copy==1
	replace `t'score  = `t'score2 - `t'gfscore0 if copy==2
	replace `t'score  = `t'score3 - `t'gfscore0 if copy==3
	replace `t'score  = `t'score4 - `t'gfscore0 if copy==4

	matrix viv = J(31,5,.)
	ivreg2 `t'score (years_anycharter = gf_*_elig_sd*_*) sd*_* d_* outg_* `covs' copy_* if `t'gfscore0!=., cluster(sasid) partial(d_* `covs' copy_* sd*_* outg_*) savefirst saverf

	est restore _ivreg2_`t'score
	matrix temp = e(b)
	forval s = 1/31 {
		matrix viv[`s',1] = temp[1,`s']
	}
	est restore _ivreg2_years_anycharter
	matrix temp = e(b)
	matrix temp2 = e(V)
	forval s = 1/31 {
		matrix viv[`s',2] = temp[1,`s']
		matrix viv[`s',3] = temp2[`s',`s']^-0.5
	}
	
	clear
	svmat viv
	scatter viv1 viv2 [w=viv3], msymbol(circle_hollow) color(black) || function y = 0.371*x, color(gs8) lwidth(thick) ra(viv2) scheme(s1color) ///
		legend(off) xline(0, lcolor(gs8) lpattern(dash)) yline(0,lcolor(gs8) lpattern(dash)) ///
		ytitle("Reduced form") xtitle("Any charter first stage", color(bg)) ylab(-2(1)2) xlab(-1(0.5)1.6) ///
		text(-1.8 1.25 "2SLS slope: 0.371" "Overid. p-val.: 0.121", justification(right)) subtitle("Math")
	graph save "${graphs}RSD_GF_VIV_`t'", replace	

	restore
	local t = "e"
	gen `t'score  = `t'score1 - `t'gfscore0 if copy==1
	replace `t'score  = `t'score2 - `t'gfscore0 if copy==2
	replace `t'score  = `t'score3 - `t'gfscore0 if copy==3
	replace `t'score  = `t'score4 - `t'gfscore0 if copy==4

	matrix viv = J(31,5,.)
	ivreg2 `t'score (years_anycharter = gf_*_elig_sd*_*) sd*_* d_* outg_* `covs' copy_* if `t'gfscore0!=., cluster(sasid) partial(d_* `covs' copy_* sd*_* outg_*) savefirst saverf

	est restore _ivreg2_`t'score
	matrix temp = e(b)
	forval s = 1/31 {
		matrix viv[`s',1] = temp[1,`s']
	}
	est restore _ivreg2_years_anycharter
	matrix temp = e(b)
	matrix temp2 = e(V)
	forval s = 1/31 {
		matrix viv[`s',2] = temp[1,`s']
		matrix viv[`s',3] = temp2[`s',`s']^-0.5
	}
	
	clear
	svmat viv
	scatter viv1 viv2 [w=viv3], msymbol(circle_hollow) color(black) || function y = 0.214*x, color(gs8) lwidth(thick) ra(viv2) scheme(s1color) ///
		legend(off) xline(0, lcolor(gs8) lpattern(dash)) yline(0,lcolor(gs8) lpattern(dash)) ///
		ytitle("Reduced form") xtitle("Any charter first stage") ylab(-2(1)2) xlab(-1(0.5)1.6) ///
		text(-1.8 1.25 "2SLS slope: 0.214" "Overid. p-val.: 0.039", justification(right)) subtitle("ELA")
	graph save "${graphs}RSD_GF_VIV_`t'", replace

	graph combine "${graphs}RSD_GF_VIV_m" "${graphs}RSD_GF_VIV_e", scheme(s1color) imargin(vsmall) rows(2) xsize(6.5) ysize(8)
	graph export "${graphs}RSD_GF_VIV.wmf", replace
}
***********************************************
* Subgroup analysis                           *
***********************************************
if `subgroup' == 1 {
	use "${cleandata}NOLA_rsdbl_gf2014", clear
	keep if matchedsample == 1

	gen gf_nox_elig = .
	gen gf_x_elig = .
	
	forval g = 4/7 {	
	
		replace gf_nox_elig = (spyear`g' == 2010 & fallschoolcode`g' == "396005" & fallgrade`g' == `g') | ///
								(spyear`g' == 2010 & fallschoolcode`g' == "396013" & fallgrade`g' == `g') | ///
								(spyear`g' == 2010 & fallschoolcode`g' == "396014" & fallgrade`g' == `g') | ///
								(spyear`g' == 2011 & fallschoolcode`g' == "396021" & fallgrade`g' == `g') | ///
								(spyear`g' == 2012 & fallschoolcode`g' == "396001" & fallgrade`g' == `g') | ///
								(spyear`g' == 2013 & fallschoolcode`g' == "396034" & fallgrade`g' == `g') | ///
								(spyear`g' == 2013 & fallschoolcode`g' == "396203" & fallgrade`g' == `g') | ///
								(spyear`g' == 2013 & fallschoolcode`g' == "396009" & fallgrade`g' == `g') | ///
								(spyear`g' == 2013 & fallschoolcode`g' == "396010" & fallgrade`g' == `g') if samp`g'==1
		
		replace gf_x_elig = (spyear`g' == 2010 & fallschoolcode`g' == "396030" & fallgrade`g' == `g') |  ///
							(spyear`g' == 2010 & fallschoolcode`g' == "396040" & fallgrade`g' == `g') if samp`g'==1
		
		replace gf_nox_elig = . if spyear`g' == . & samp`g'==1	
		replace gf_x_elig = . if spyear`g' == . & samp`g'==1	
	}
	bys match_cell: egen nox_samp = max(gf_nox_elig)
	bys match_cell: egen x_samp = max(gf_x_elig)
	
	expand 4
	bys sasid: gen copy = _n
	qui tab copy, gen(copy_)
	gen years_any = years_any1 if copy==1
	replace years_any = years_any2 if copy==2
	replace years_any = years_any3 if copy==3
	replace years_any = years_any4 if copy==4
	
	gen grade     = 5 if (samp4==1 & copy==1)
	replace grade = 6 if (samp4==1 & copy==2) | (samp5==1 & copy==1)
	replace grade = 7 if (samp4==1 & copy==3) | (samp5==1 & copy==2) | (samp6==1 & copy==1)
	replace grade = 8 if (samp4==1 & copy==4) | (samp5==1 & copy==3) | (samp6==1 & copy==2) | (samp7==1 & copy==1)
	qui tab grade, gen(outg_)
	
	gen yearout = .
	forval grade = 5/8 {
		replace yearout = spyear`grade' if grade==`grade'
	}
	qui tab yearout, gen(y_)

	local row = 1
	local col = 1
	
	gen male = 1-female
	egen bl_me = rowmean(bl_m bl_e)
	xtile bl_me_terc = bl_me, nq(3)
	gen hi_bl = bl_me_terc == 3
	gen mid_bl = bl_me_terc == 2
	gen lo_bl = bl_me_terc == 1
	
	foreach t in m e {
	
		gen `t'score  = `t'score1 if copy==1
		replace `t'score = `t'score2 if copy==2
		replace `t'score = `t'score3 if copy==3
		replace `t'score = `t'score4 if copy==4
		gen `t'level = `t'score
		replace `t'score = `t'score1-`t'gfscore0 if copy==1		
		replace `t'score = `t'score2-`t'gfscore0 if copy==2
		replace `t'score = `t'score3-`t'gfscore0 if copy==3
		replace `t'score = `t'score4-`t'gfscore0 if copy==4	
		
		foreach sub of varlist female male hi_bl mid_bl lo_bl nox_samp x_samp {
			summ `t'level if gf_any_elig==0 & `t'gfscore0!=. & `sub'==1
					matrix gfing[`row',1]   = r(mean)
					matrix gfing[`row'+2,1] = r(N)

			/* OLS */
			reg `t'score years_any  d_* `covs' copy_* outg_* if `t'gfscore0!=. & `sub'==1, cluster(sasid)
				matrix gfing[`row',2]   = _b[years_any]
				matrix gfing[`row'+1,2] = _se[years_any]
				matrix gfing[`row'+2,2] = e(N)	
				
			/* RF */
			areg `t'score gf_any_elig `covs' copy_* outg_* if years_any!=. & `t'gfscore0!=. & `sub'==1, absorb(match_cell) cluster(sasid)
				matrix gfing[`row',3]   = _b[gf_any_elig]
				matrix gfing[`row'+1,3] = _se[gf_any_elig]
				matrix gfing[`row'+2,3] = e(N)
			
			/* FS */
			areg years_any gf_any_elig `covs' copy_* outg_* if `t'score!=. & `t'gfscore0!=. & `sub'==1, absorb(match_cell) cluster(sasid)
				matrix gfing[`row',4]   = _b[gf_any_elig]
				matrix gfing[`row'+1,4] = _se[gf_any_elig]
				matrix gfing[`row'+2,4] = e(N)
				
			/* 2SLS */
			ivregress 2sls `t'score (years_any = gf_any_elig) d_* `covs' copy_* outg_* if `t'gfscore0!=. & `sub'==1, cluster(sasid)
				matrix gfing[`row',5]   = _b[years_any]
				matrix gfing[`row'+1,5] = _se[years_any]
				matrix gfing[`row'+2,5] = e(N)
				
			local row = `row' + 3	
		}
		
	}
	clear
	svmat gfing
	br

}
***********************************************
* Fallback with alternative interactions      *
***********************************************
if `altfbchange' == 1 {
	use "${cleandata}NOLA_rsdbl_gf2014", clear
	keep if matchedsample == 1

	expand 4
	bys sasid: gen copy = _n
	qui tab copy, gen(copy_)
	gen years_any = years_any1 if copy==1
	replace years_any = years_any2 if copy==2
	replace years_any = years_any3 if copy==3
	replace years_any = years_any4 if copy==4
	
	gen grade     = 5 if (samp4==1 & copy==1)
	replace grade = 6 if (samp4==1 & copy==2) | (samp5==1 & copy==1)
	replace grade = 7 if (samp4==1 & copy==3) | (samp5==1 & copy==2) | (samp6==1 & copy==1)
	replace grade = 8 if (samp4==1 & copy==4) | (samp5==1 & copy==3) | (samp6==1 & copy==2) | (samp7==1 & copy==1)
	
	qui tab grade, gen(outg_)
	
	gen yearout = .
	forval grade = 5/8 {
		replace yearout = spyear`grade' if grade==`grade'
	}
	qui tab yearout, gen(y_)
	
	forval out = 1/4 {
		replace incharter`out' = incharter`out'*(1-attend_any`out')
	}
	gen years_charter1 = incharter1
	gen temp = incharter1*(1+repeatout1)
	egen years_charter2 = rowtotal(temp incharter2)
	replace years_charter2 = . if incharter1==. & incharter2==.
	gen temp2 = incharter2*(1+repeatout2)
	egen years_charter3 = rowtotal(temp temp2 incharter3)
	replace years_charter3 = . if incharter1==. & incharter2==. & incharter3==.	
	gen temp3 = incharter3*(1+repeatout3)
	egen years_charter4 = rowtotal(temp temp2 temp3 incharter4)
	replace years_charter4 = . if incharter1==. & incharter2==. & incharter3==. & incharter4==.			
	drop temp temp2 temp3

	gen years_charter = years_charter1 if copy==1
	replace years_charter = years_charter2 if copy==2
	replace years_charter = years_charter3 if copy==3
	replace years_charter = years_charter4 if copy==4
		
	gen years_anycharter = years_any + years_charter
	
	gen baselineyr = .
	gen baselinespsbin = .
	gen baselinesps = .
	forval gr = 4/7 {
		local lgr = `gr' - 1
		replace baselineyr = spyear`lgr' if samp`gr'==1 
		replace baselinespsbin = int(sps`lgr'/5) if samp`gr'==1
		replace baselinesps = sps`lgr' if samp`gr'==1
		
	}

	preserve
	
	egen temp = group(bl_sped baselinespsbin samp4 samp5 samp6 samp7 bl_year) 
	qui tab temp, gen(sd1_)

	local row = 1
	local col = 1
	
	foreach t in m e {
	
		gen `t'score  = `t'score1 - `t'gfscore0 if copy==1
		replace `t'score  = `t'score2 - `t'gfscore0 if copy==2
		replace `t'score  = `t'score3 - `t'gfscore0 if copy==3
		replace `t'score  = `t'score4 - `t'gfscore0 if copy==4

		ivreg2 `t'score (years_any = gf_any_elig) d_* `covs' copy_* outg_* if `t'gfscore0!=., cluster(sasid) partial(d_* `covs' copy_* outg_*) savefirst first
		gen samp = e(sample)
		bys temp: egen C_hat = mean(years_charter) if e(sample) == 1
		bys temp: egen N = count(years_charter) if e(sample) == 1
		replace C_hat = (C_hat*N - years_charter)/(N-1) if e(sample) == 1
		gen gf_elig_C_hat = gf_any_elig * C_hat	

	local ++col
	ivreg2 `t'score (years_any = gf_any_elig gf_elig_C_hat) sd*_* d_* outg_* `covs' copy_* if `t'gfscore0!=., cluster(sasid) partial(d_* `covs' copy_* sd*_* outg_*) savefirst first
		matrix temp = e(first)
		matrix fallbackchange[`row',`col']    = _b[years_any]
		matrix fallbackchange[`row'+1,`col']  = _se[years_any]
		matrix fallbackchange[`row'+2,`col']  = temp[7,1]
		matrix fallbackchange[`row'+9,`col']  = e(exexog_ct)
		matrix fallbackchange[`row'+10,`col'] = e(widstat)
		matrix fallbackchange[`row'+11,`col'] = e(N)
		
	local ++col
	ivreg2 `t'score (years_any years_charter = gf_any_elig gf_elig_C_hat) sd*_* d_* outg_* `covs' copy_* if `t'gfscore0!=., cluster(sasid) partial(d_* `covs' copy_* sd*_* outg_*) savefirst first
		matrix temp = e(first)
		matrix fallbackchange[`row',`col']    = _b[years_any]
		matrix fallbackchange[`row'+1,`col']  = _se[years_any]
		matrix fallbackchange[`row'+2,`col']  = temp[7,1]
		matrix fallbackchange[`row'+3,`col']  = _b[years_charter]
		matrix fallbackchange[`row'+4,`col']  = _se[years_charter]
		matrix fallbackchange[`row'+5,`col']  = temp[7,2]
		matrix fallbackchange[`row'+9,`col']  = e(exexog_ct)
		matrix fallbackchange[`row'+10,`col'] = e(widstat)
		matrix fallbackchange[`row'+11,`col'] = e(N)

	local ++col

	ivreg2 `t'score (years_anycharter = gf_any_elig gf_elig_C_hat) sd*_* d_* `covs' copy_* outg_* if `t'gfscore0!=., cluster(sasid) partial(d_* `covs' copy_* outg_* sd*_*) savefirst first
		matrix temp = e(first)
		matrix fallbackchange[`row'+6,`col']  = _b[years_anycharter]
		matrix fallbackchange[`row'+7,`col']  = _se[years_anycharter]	
		matrix fallbackchange[`row'+8,`col']  = temp[7,1]
		matrix fallbackchange[`row'+9,`col']  = e(exexog_ct)
		matrix fallbackchange[`row'+10,`col'] = e(widstat)
		matrix fallbackchange[`row'+11,`col'] = e(N)

	drop samp C_hat gf_elig_C_hat N
		
	local row = `row' + 12
	local col = 1
	}

	restore
	
	gen bl_schoolcode = ""
	gen bl_grade = .
	forval gr = 4/7 {
		local bl_gr = `gr' - 1
		replace bl_grade = `bl_gr' if samp`gr' == 1
		replace bl_schoolcode = schoolcode`bl_gr' if samp`gr' == 1
		
	}
	destring bl_schoolcode, replace

	merge m:1 bl_schoolcode bl_year using "${cleandata}RSD_closest_school_distance_short3"
	
	gen gf_elig_dist = gf_any_elig * m_chart_dist
	gen m_chart_dist2 = m_chart_dist^2
	gen gf_elig_dist2 = gf_any_elig * m_chart_dist^2
	
	egen temp = group(bl_sped baselinespsbin) 
	qui tab temp, gen(sd1_)
	drop temp
	egen temp = group(samp4 samp5 samp6 samp7)
	qui tab temp, gen(sd2_)
	drop temp
	qui tab bl_year, gen(sd3_)
	foreach var1 of varlist sd*_* {
		foreach var2 of varlist gf_any_elig {
			gen `var2'_`var1' = `var2'*`var1'
		}
	}
	
	
	local row = 1
	local col = 7

	foreach t in m e {
	
		gen `t'score  = `t'score1 - `t'gfscore0 if copy==1
		replace `t'score  = `t'score2 - `t'gfscore0 if copy==2
		replace `t'score  = `t'score3 - `t'gfscore0 if copy==3
		replace `t'score  = `t'score4 - `t'gfscore0 if copy==4

	ivreg2 `t'score (years_any = gf_any_elig gf_elig_dist*) m_chart_dist* d_* outg_* `covs' copy_* if `t'gfscore0!=., cluster(sasid) partial(d_* `covs' copy_* outg_*) savefirst first
		matrix temp = e(first)
		matrix fallbackchange[`row',`col']    = _b[years_any]
		matrix fallbackchange[`row'+1,`col']  = _se[years_any]
		matrix fallbackchange[`row'+2,`col']  = temp[7,1]
		matrix fallbackchange[`row'+9,`col']  = e(exexog_ct)
		matrix fallbackchange[`row'+10,`col'] = e(widstat)
		matrix fallbackchange[`row'+11,`col'] = e(N)
	
	local ++col
	ivreg2 `t'score (years_any years_charter = gf_any_elig gf_elig_dist*) m_chart_dist* d_* outg_* `covs' copy_* if `t'gfscore0!=., cluster(sasid) partial(d_* `covs' copy_* outg_*) savefirst first
		matrix temp = e(first)
		matrix fallbackchange[`row',`col']    = _b[years_any]
		matrix fallbackchange[`row'+1,`col']  = _se[years_any]
		matrix fallbackchange[`row'+2,`col']  = temp[7,1]
		matrix fallbackchange[`row'+3,`col']  = _b[years_charter]
		matrix fallbackchange[`row'+4,`col']  = _se[years_charter]
		matrix fallbackchange[`row'+5,`col']  = temp[7,2]
		matrix fallbackchange[`row'+9,`col']  = e(exexog_ct)
		matrix fallbackchange[`row'+10,`col'] = e(widstat)
		matrix fallbackchange[`row'+11,`col'] = e(N)

	local ++col
	ivreg2 `t'score (years_anycharter = gf_any_elig gf_elig_dist*) m_chart_dist* d_* `covs' copy_* outg_* if `t'gfscore0!=., cluster(sasid) partial(d_* `covs' copy_* outg_*) savefirst first
		matrix temp = e(first)
		matrix fallbackchange[`row'+6,`col']  = _b[years_anycharter]
		matrix fallbackchange[`row'+7,`col']  = _se[years_anycharter]	
		matrix fallbackchange[`row'+8,`col']  = temp[7,1]
		matrix fallbackchange[`row'+9,`col']  = e(exexog_ct)
		matrix fallbackchange[`row'+10,`col'] = e(widstat)
		matrix fallbackchange[`row'+11,`col'] = e(N)	
		
	local ++col
	ivreg2 `t'score (years_any = gf_elig_dist* gf_*_elig_sd*_*) m_chart_dist*  sd*_* d_* outg_* `covs' copy_* if `t'gfscore0!=., cluster(sasid) partial(d_* `covs' copy_* outg_*  sd*_*) savefirst first
		matrix temp = e(first)
		matrix fallbackchange[`row',`col']    = _b[years_any]
		matrix fallbackchange[`row'+1,`col']  = _se[years_any]
		matrix fallbackchange[`row'+2,`col']  = temp[7,1]
		matrix fallbackchange[`row'+9,`col']  = e(exexog_ct)
		matrix fallbackchange[`row'+10,`col'] = e(widstat)
		matrix fallbackchange[`row'+11,`col'] = e(N)
	
	local ++col
	ivreg2 `t'score (years_any years_charter = gf_elig_dist* gf_*_elig_sd*_*) m_chart_dist*  sd*_* d_* outg_* `covs' copy_* if `t'gfscore0!=., cluster(sasid) partial(d_* `covs' copy_* outg_*  sd*_*) savefirst first
		matrix temp = e(first)
		matrix fallbackchange[`row',`col']    = _b[years_any]
		matrix fallbackchange[`row'+1,`col']  = _se[years_any]
		matrix fallbackchange[`row'+2,`col']  = temp[7,1]
		matrix fallbackchange[`row'+3,`col']  = _b[years_charter]
		matrix fallbackchange[`row'+4,`col']  = _se[years_charter]
		matrix fallbackchange[`row'+5,`col']  = temp[7,2]
		matrix fallbackchange[`row'+9,`col']  = e(exexog_ct)
		matrix fallbackchange[`row'+10,`col'] = e(widstat)
		matrix fallbackchange[`row'+11,`col'] = e(N)

	local ++col
	ivreg2 `t'score (years_anycharter = gf_elig_dist* gf_*_elig_sd*_*) m_chart_dist*  sd*_* d_* `covs' copy_* outg_* if `t'gfscore0!=., cluster(sasid) partial(d_* `covs' copy_* outg_*  sd*_*) savefirst first
		matrix temp = e(first)
		matrix fallbackchange[`row'+6,`col']  = _b[years_anycharter]
		matrix fallbackchange[`row'+7,`col']  = _se[years_anycharter]	
		matrix fallbackchange[`row'+8,`col']  = temp[7,1]
		matrix fallbackchange[`row'+9,`col']  = e(exexog_ct)
		matrix fallbackchange[`row'+10,`col'] = e(widstat)
		matrix fallbackchange[`row'+11,`col'] = e(N)
		
	local row = `row' + 12
	local col = 7
	}
	
	
clear
svmat fallbackchange
br
}
