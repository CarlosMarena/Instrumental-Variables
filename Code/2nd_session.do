cls 
clear all 

local course_path "E:\OneDrive - UDEP\Cursos\Mixtape Sessions\Instrumental-Variables"
local data_path "`course_path'\Data"

use "`data_path'\stevenson", clear 

* Excercise 1 
keep if black == 1 
egen leniency = mean(jail3) , by(judge_pre_1-judge_pre_7)
summ leniency , d 
gen more_lenient = leniency > r(p50)

/* OLS Regression */
reg guilt jail3, r 

/* 2SLS Regression */
ivregress 2sls guilt (jail3=more_lenient), r 

