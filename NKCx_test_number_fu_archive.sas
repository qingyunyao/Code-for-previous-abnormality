/*File name..: find the number of test taken ever since the first HPV test */
/*Study......: Phd project study 2 HPV-based screening among women with previous abnormality*/
/*Author.....: Qingyun Yao*/
/*Date.......: 2025/06/24*/
/*Updated....:*/
/*Purpose....: find the follow_up test after the first hpv date*/
/*Note.......: population defination: women's first valid HPV-based screening test from 2012 to the end of the date registered in the
registry (use hpv data). Women aged 50-70 at the time of their first valid HPV-based screening. Women without any screening test for the past 4 years (both ext and hpv, cell). 
Women had no hysterectomy before the first screening test (total hysterectomy or person).
test age was calculated by test_year-sample_year
*/
*------------------------------------------------------------------------;
/* Data used...:Socmob4.nkc_hpv_2022.dta
Socmob4.nkc_trans_cell_2022.dta
Socmob4.nkc_ext_hpv_2022.dta
Socmob4.nkc_pad_translated_2022.dta
Socmob4.nkc_person _2022.dta
Socmob4.nkc_deregister_2022.dta

libname NCSR  odbc complete="server=meb-sql02.meb.ki.se;driver=SQL Server Native Client 11.0;Trusted_Connection=Yes;database=NCSR" schema=NCSR ;
libname V_NCSR     odbc complete="server=meb-sql02.meb.ki.se;driver=SQL Server Native Client 11.0;Trusted_Connection=Yes;database=NCSR" schema=V_NCSR ;
libname PreAbn 'P:\ACCES\ACCES_Research\Qingyun\previous abnormality\Data';
libname SOCMOB4  odbc complete="dsn=kosmos;database=CERVIX_SOCMOB4" schema=CLEAN_2022 ;

There is one systematic mistake in ext_hpv in lab 411, x_sample_year=2016, x_sample_date were in 2002 delete it when analyzing

/* Data created.:  follow_up_20241015.dta*/
/*sas version.: SAS9.4*/
/*main program*/

libname NCSR odbc noprompt="dsn=kosmos; database=NCSR" schema=NCSR;
libname V_NCSR odbc noprompt="dsn=kosmos; database=NCSR" schema=V_NCSR ;
libname PreAbn 'P:\ACCES\ACCES_Research\Qingyun\previous abnormality\Data';


proc datasets library=work kill;quit;

data pop;
set preabn.NKCx_nohys_cohort_260202;
keep person_id first_hpv_date HPVrisk  noscr_his previous_pos previous_pos_cat previous_pos_date;
where noscr_his=0;
run;

/*all test*/

data hpv_24;
set V_ncsr.nkc_hpv_2024;
format test_date yymmdd10. regis_date yymmdd10.;
sample_year=input(substr(x_sample_date,1,4),8.);
sample_month=input(substr(x_sample_date,6,2),8.);
sample_day=input(substr(x_sample_date,9,2),8.);
test_date=MDY(sample_month,sample_day,sample_year);

reg_year=input(substr(x_reg_date,1,4),8.);
reg_month=input(substr(x_reg_date,6,2),8.);
reg_day=input(substr(x_reg_date,9,2),8.);
regis_date=MDY(reg_month,reg_day,reg_year);

if selftest=1 and x_reg_date^='' then test_date=regis_date;
keep person_id test_date hpvdiag hpv_type regis_date x_sample_yr;
rename regis_date=reg_date;
run;

data hpv_ext_24;
set V_ncsr.nkc_ext_hpv_2024;
format test_date yymmdd10. reg_date yymmdd10.;
sample_year=input(substr(x_sample_date,1,4),8.);
sample_month=input(substr(x_sample_date,6,2),8.);
sample_day=input(substr(x_sample_date,9,2),8.);
test_date=MDY(sample_month,sample_day,sample_year);

reg_year=input(substr(x_reg_date,1,4),8.);
reg_month=input(substr(x_reg_date,6,2),8.);
reg_day=input(substr(x_reg_date,9,2),8.);
reg_date=MDY(reg_month,reg_day,reg_year);

keep person_id test_date hpvdiag hpv_type reg_date x_sample_yr;
run;

data hpvall;
set hpv_24 hpv_ext_24;
run;

data hpvall;
set hpvall;
test_year=year(test_date);
reg_year=year(reg_date);
if abs(test_year-x_sample_yr)>1 then check=1;
if check=1 and reg_year=x_sample_yr then do; test_date=reg_date; change=1;end;
run;

proc freq data=hpvall;
tables change*check/missing;
run;


proc freq data=hpvall;
table hpv_type;
run;

data hpvall;
set hpvall;
HPV16=index(hpv_type,'16')>0;
HPV18=index(hpv_type,'18')>0;
if HPV_type in ('18 45' '1845')then HPV18_45=1;
HPV45=index(hpv_type,'45')>0;
HPV33=index(hpv_type,'33')>0;
HPV58=index(hpv_type,'58')>0;
HPV31=index(hpv_type,'31')>0;
HPV52=index(hpv_type,'52')>0;
HPV35=index(hpv_type,'35')>0;
HPV39=index(hpv_type,'39')>0;
HPV51=index(hpv_type,'51')>0;
HPV59=index(hpv_type,'59')>0;
HPV56=index(hpv_type,'56')>0;
HPV68=index(hpv_type,'68')>0;
HPV66=index(hpv_type,'66')>0;
if hpv_type=:'Dålig' then HPVmissing=1;
if hpv_type=:'Ej' then HPVmissing=1;
if hpv_type=:'ej' then HPVmissing=1;
if hpv_type=:'HPV DNA extracted' then HPVmissing=1;/*results with this hpv_type are all negative*/
if hpv_type=:'HPV Other' then do; HPV16=0;HPV18=0;HPVHr=1;end;
if hpv_type=:'Other HPV' then do; HPV16=0;HPV18=0;HPVHr=1;end;
if hpv_type=:'high risk type' then HPVHr=1;
if hpv_type=:'HPV high risk type' then HPVHr=1;
if hpv_type=:'HPV other' then HPVHr=1;
if hpv_type=:'HPV type missing' then HPVHr=1;
if hpv_type=:'HPV type unknown' then HPVHr=1;
if hpv_type=:'HPV unknown' then HPVHr=1;
if hpv_type='HPV Övrig' then HPVHr=1;
if hpv_type='HPVOVR' then HPVHr=1;
if hpv_type='Not available' then HPVmissing=1;
if hpv_type=:'Irrelevant' then HPVmissing=1;
if hpv_type=:'Se' then HPVmissing=1;
if hpv_type=:'human' then HPVmissing=1;
if hpv_type=:'brunn' then HPVmissing=1;
if HPV18_45=1 then do; HPV18=0; HPV45=0;end;
if HPV33=1 or HPV58=1 or HPV45=1 or HPV31=1 or HPV52=1 or HPV35=1 or HPV39=1 or HPV51=1 or HPV59=1 or HPV56=1 or HPV68=1 or HPV66=1 then
HPVHr=1;
run;

data hpvall;
set hpvall;
if HPV16=1 and HPVDIAG='POS' then do; HPVrisk=4;Type16=1;end;
else if HPV18=1 and HPVDIAG='POS' then do; HPVrisk=3;Type18=1;end;
else if HPV18_45=1 and HPVDIAG='POS' then do;HPVrisk=3; Type18_45=1;end;
else if HPVHr=1 and HPVDIAG='POS' then do; HPVrisk=2; TypeHr=1;end;
else if HPVmissing=1 and HPVDIAG='POS' then HPVrisk=.;
else if HPVDIAG='POS' then HPVrisk=0;
else if HPVDIAG^='NEG' then HPVrisk=.;
else if HPVDIAG='NEG' then HPVrisk=0;
drop HPV16 HPV18 HPV45 HPV33 HPV58 HPV31 HPV52 HPV35 HPV39 HPV51 HPV59 HPV56 HPV68 HPV66;
run;



data hpvtest;
set hpvall;
if hpvrisk>=2 then outcome=1;
else outcome=0;
if test_date=. and reg_date=. then delete;
if check=1 and change=. then delete;
if hpvrisk=. then delete;
keep person_id test_date outcome reg_date;
run; 

proc sort data=hpvtest;
by person_id test_date descending outcome;
run;
proc sort data=hpvtest nodupkey;
by person_id test_date;
run;
/*NOTE: There were 7031769 observations read from the data set WORK.HPVTEST.
NOTE: 1868502 observations with duplicate key values were deleted.
NOTE: The data set WORK.HPVTEST has 5163267 observations and 4 variables.
*/

data cytotest;
set V_ncsr.nkc_trans_cell_2024;
format test_date yymmdd10. regis_date yymmdd10.;
sample_year=input(substr(x_sample_date,1,4),8.);
sample_month=input(substr(x_sample_date,6,2),8.);
sample_day=input(substr(x_sample_date,9,2),8.);
test_date=MDY(sample_month,sample_day,sample_year);

reg_year=input(substr(x_reg_date,1,4),8.);
reg_month=input(substr(x_reg_date,6,2),8.);
reg_day=input(substr(x_reg_date,9,2),8.);
regis_date=MDY(reg_month,reg_day,reg_year);
if snomed_severity>5 then outcome=1;
else outcome=0;
/*should also delete the unavailable ones*/
keep person_id test_date outcome regis_date x_sample_yr;
rename regis_date=reg_date;
run;


data cytotest;
set cytotest;
test_year=year(test_date);
reg_year=year(reg_date);
if abs(test_year-x_sample_yr)>1 then check=1;
if check=1 and reg_year=x_sample_yr then do; test_date=reg_date; change=1;end;
run;
proc freq data=cytotest;
tables check*change/missing;
run;

data cytotest;
set cytotest;
if test_date=0 and reg_date=0 then delete;
if check=1 and change=. then delete;
keep person_id test_date outcome  reg_date;
run; 
/*NOTE: There were 25693461 observations read from the data set WORK.CYTOTEST.
NOTE: The data set WORK.CYTOTEST has 25690724 observations and 4 variables.
*/

proc sort data=cytotest;
by person_id test_date decending outcome;
run;
proc sort data=cytotest nodupkey;
by person_id test_date;
run;
/*NOTE: There were 25690724 observations read from the data set WORK.CYTOTEST.
NOTE: 109866 observations with duplicate key values were deleted.
NOTE: The data set WORK.CYTOTEST has 25580858 observations and 4 variables.
*/


data hpvtest;
set hpvtest;
testmethod=1;
run;
data cytotest;
set cytotest;
testmethod=0;
run;
/*can ignore cytology test*/
/*keep one test per one day*/
data Alltest;
set cytotest hpvtest;
sample_year=year(test_date);
run;
data alltest;
set alltest;
if sample_year<=2024;
run;

proc sort data=alltest;
by person_id test_date descending testmethod descending outcome;
run;
proc sort data=alltest nodupkey;
by person_id test_date testmethod;
run;


/*find the test number between previous abnormality and first hpv test*/
data test_between;
merge pop(in=a) alltest;
by person_id;
if previous_pos_date<test_date<first_hpv_date;
if a;
run;

proc sql;
create table between_test_num as select person_id, count(person_id) as num_test_bet from test_between group by person_id;
quit;

data pop_between_test;
merge pop (in=a) between_test_num;
by person_id;
if a;
run;

data pop_between_test;
set pop_between_test;
if num_test_bet=. then num_test_bet=0;
run;
proc freq data=pop_between_test;
tables previous_pos_cat/missing;
where num_test_bet=0;
run;

/*find all the test after first hpv test*/
data followup_all;
merge pop (in=a) alltest;
by person_id;
if test_date>first_hpv_date;
if a;
run;

proc sql;
create table follow_up_count as select person_id, count(person_id) as all_count, sum (testmethod) as num_hpv from followup_all group by person_id;
quit;
data follow_up_count;
set follow_up_count;
num_cyto=all_count-num_hpv;
run;


/*only use hpv test*/
proc sort data=hpvtest;
by person_id test_date descending outcome;
run;
proc sort data=hpvtest nodupkey;
by person_id test_date;
run;


/*find the test after first hpv test*/
data followup_hpv;
merge pop (in=a) hpvtest;
by person_id;
if test_date>first_hpv_date;
if a;
run;

proc sql;
create table hpv_count as select person_id, count(person_id) as hpv_count, sum(outcome) as followupresult from followup_hpv group by person_id;
quit;

proc sort data=followup_hpv;
by person_id descending test_date;
run;
proc sort data=followup_hpv nodupkey out=lasthpv;
by person_id;
run;

data hpv_numbers;
merge  lasthpv hpv_count;
by person_id;
rename test_date=last_hpv_date outcome=last_hpv_outcome;
drop reg_date;
run;

data hpv_numbers;
set hpv_numbers;
if HPVrisk>1 and  hpv_count-followupresult>0 then turn_negative=1;
else if HPVrisk>1 then turn_negative=0;
if HPVrisk=0 and followupresult>0 then turn_positive=1;
else if HPVrisk=0 then turn_positive=0;
followup=1;
run;


/*find first hpv test turn negative for positive women and find first positive test for negative women and keep the date*/

data positive;
set followup_hpv;
if HPVrisk>1 and outcome=0;
run;
proc sort data=positive;
by person_id test_date;
run;
proc sort data=positive nodupkey out=turnN;
by person_id ;
run;
data turnn;
set turnn;
keep person_id test_date;
rename test_date=turn_n_date;
run;


data negative;
set followup_hpv;
if HPVrisk=0 and outcome=1;
run;
proc sort data=negative;
by person_id test_date;
run;
proc sort data=negative nodupkey out=turnP;
by person_id ;
run;
data turnP;
set turnP;
keep person_id test_date;
rename test_date=turn_p_date;
run;








/*last test*/

proc sort data=followup_all;
by person_id descending test_date;
run;
proc sort data=followup_all nodupkey out=lastfollowup;
by person_id;
run;

data lastfollowup;
set lastfollowup;
rename test_date=last_test_date outcome=last_outcome testmethod=lastMethod;
drop reg_date;
run;

data followup_numbers;
merge  lastfollowup follow_up_count hpv_numbers;
by person_id;
followup=1;
run;


/*keep all followup test*/
proc transpose data=followup_all out=date prefix=followup_date_;
var test_date;
by person_id;
run;
proc transpose data=followup_all out=outcome prefix=followup_outcome_;
var outcome;
by person_id;
run;
proc transpose data=followup_all out=method prefix=followup_method_;
var testmethod;
by person_id;
run;


data followup_all;
merge pop (in=a) followup_numbers turnn turnP date outcome method;
by person_id;
drop _NAME_;
run;
data followup_all;
set followup_all;
if HPVrisk=0 then HPV=0;
else HPV=1;
follow_up_time=(last_test_date-first_hpv_date)/365.25;
run;

proc freq data=followup_all;
tables HPV*turn_negative HPV*turn_positive testmethod*HPV/missing;
run;


proc summary data=followup_all;
var all_count num_hpv num_cyto follow_up_time;
class turn_positive/missing;
output out=num median= q1=q3=/autoname;
where HPV=0 ;
run;


data pop;
set preabn.NKCx_nohys_cohort_260202;
keep person_id HPVrisk noscr_his previous_pos previous_pos_cat;
run;

data followup_all_merge;
merge pop followup_all;
by person_id;
run;
proc freq data=followup_all_merge;
table  previous_pos*noscr_his/missing;
run;
data followup_all_merge;
set followup_all_merge;
if previous_pos=0 and noscr_his=0 then ref=1;
else ref=0;
run;
proc freq data=followup_all_merge;
tables HPV*turn_negative HPV*turn_positive testmethod*HPV/missing;
where noscr_his=1;
run;

proc summary data=followup_all_merge;
var all_count follow_up_time;
*class turn_positive/missing;
output out=num median= q1=q3=/autoname;
where HPV=0 and ref=1 and testmethod=1;
run;

data preabn.NKcx_nohys_follow_up_test;
set followup_all_merge;
run;



/*table*/
data followup_all_merge;
set followup_all_merge;
if all_count=. then all_count=0;
run;

proc freq data=followup_all_merge;
table ref;
where HPVrisk=0;
run;

proc summary data=followup_all_merge;
var all_count;
class ref previous_pos previous_pos_cat noscr_his;
ways 1;
where HPV=1;
output out=fu_count_sum median= q1= q3=/autoname;
run;



/*supplementary analysis*/
data analysis;
set preabn.NKcx_nohys_follow_up_test;
where noscr_his=0;
run;

data analysis;
set analysis;
if followup=. then followup=0;
run;
data analysis;
set analysis;
format firstfollow_up yymmdd10.;
if followup=1 then firstfollow_up=min(of followup_date_:);
run;

data analysis;
set analysis;
if followup=1 then diff_first=(firstfollow_up-first_hpv_date)/365.25;
run;

proc summary data=analysis;
var diff_first;
class HPV previous_pos_cat previous_pos;
output out=fu_median median()= q1()= q3()=/autoname;
run;

proc summary data=analysis;
var followup;
class HPV previous_pos_cat;
*where HPV=0;
output out=fu_count_sum_neg_all N()= sum()= mean()= /autoname;
run;

data fu_count_sum_neg_all;
set fu_count_sum_neg_all;
nfu=cat(compress(put(followup_sum,8.0)),' (', compress(put(followup_mean*100,8.1)),'%)');
run;
%let mydir=P:\ACCES\ACCES_Research\Qingyun\previous abnormality\Doc\NKCx_260604_exnsh\;
title 'Sup table';
ods rtf file="&mydir.fu_test_num_&sysdate..rtf";
proc print data=fu_count_sum_neg_all noobs;run;
title '';
ods rtf close;
