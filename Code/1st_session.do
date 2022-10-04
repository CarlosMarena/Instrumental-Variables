cls 
clear all 

local course_path "E:\OneDrive - UDEP\Cursos\Mixtape Sessions\Instrumental-Variables"
local data_path "`course_path'\Data"

use "`data_path'\angrist_krueger", clear 

* Excercise 1 
reg lwage educ, r 
scatter (lwage) (educ)

* Excercise 2 
ivregress 2sls lwage (educ=1.qob), r 

* Excercise 3 
gen firstqtr = qob == 1  

reg lwage firstqtr, r 
mean lwage , over(firstqtr)

reg educ firstqtr, r
mean educ , over(firstqtr)

* Excercise 4 
gen scondqtr = qob == 2
gen thirdqtr = qob == 3 

ivregress 2sls lwage (educ = firstqtr scondqtr thirdqtr), r first 
