/*File name..:  participation index*/
/*Study......: Phd project study 2 the participation before first HPV test age (50-70) in NKCx*/
/*Author.....: Qingyun Yao*/
/*Date.......: 2025.06-23*/
/*Updated....: */
/*Purpose....: /*add the screening compliance? index using registry data*/
/*from 1995 or age 23
because the NKCx has the all coverage of test since 1995 and women started to have screening from the age of 23*/
/*find all the test that women take in the time period of first test- to the first HPV date,*/
/*find earliest test date in NKCx system and mark the age of first screening test, and count from that test
screening time interval:

/*policy with cytology before 50 3 years interval after 50 5 years interval
hpv before 50 5 years interval after 50 7 years interval
if the outcome=1 then all within 1 year interval*/

/*


mark the screening start age compare to age 23
start counting screening participation from age 23
if first screening age>=23 screeningstartage-23=unscreened time
if first screening age<23 screeningstartage=23, start from the test that has test_age+interval>23 otherwise treat it as first category
*/

/*sort it by sample date*/
/*
libname NCSR  odbc complete="server=meb-sql02.meb.ki.se;driver=SQL Server Native Client 11.0;Trusted_Connection=Yes;database=NCSR" schema=NCSR ;
libname V_ncsr     odbc complete="server=meb-sql02.meb.ki.se;driver=SQL Server Native Client 11.0;Trusted_Connection=Yes;database=NCSR" schema=V_ncsr ;
libname PreAbn 'P:\ACCES\ACCES_Research\Qingyun\previous abnormality\Data';

				*/
/*Note.......: */
*------------------------------------------------------------------------;
/* Data used...: socmob4.nkc_hpv socmob4.nkc_trans_cell  self_samp_pop;

/* Data created.: pop_index_&date */
proc datasets library=work kill;quit;

libname NCSR odbc noprompt="dsn=kosmos; database=NCSR" schema=NCSR;
libname V_NCSR odbc noprompt="dsn=kosmos; database=NCSR" schema=V_NCSR ;
libname PreAbn 'P:\ACCES\ACCES_Research\Qingyun\previous abnormality\Data';




/*add the participation index
when calculating this, should including invalid test, based on behavior rather than the protection effect*/

/*pop*/
data pop;
set preabn.NKCx_nohys_cohort_260202;
keep person_id birth_date first_hpv_date;
run;
/*
NOTE: There were 742085 observations read from the data set PREABN.NKCX_NOHYS_COHORT_260202.
NOTE: The data set WORK.POP has 742085 observations and 3 variables.
*/


/*all test in the dataset*/

data hpv_24;
set v_ncsr.nkc_hpv_2024;
if selftest=1 and x_reg_date^='' then x_sample_date=x_reg_date;
keep person_id x_sample_date hpvdiag hpv_type x_reg_date x_sample_yr;
run;
data hpv_24;
set hpv_24;
format test_date yymmdd10. reg_date yymmdd10.;
sample_year=input(substr(x_sample_date,1,4),8.);
sample_month=input(substr(x_sample_date,6,2),8.);
sample_day=input(substr(x_sample_date,9,2),8.);
reg_year=input(substr(x_reg_date,1,4),8.);
reg_month=input(substr(x_reg_date,6,2),8.);
reg_day=input(substr(x_reg_date,9,2),8.);
test_date=MDY(sample_month,sample_day,sample_year);
reg_date=MDY(reg_month,reg_day,reg_year);
drop sample_year sample_month sample_day reg_year reg_month reg_day x_sample_date x_reg_date;
run;

data hpv_ext_24;
set  v_ncsr.nkc_ext_hpv_2024;
keep person_id x_sample_date hpvdiag hpv_type x_reg_date x_sample_yr;
run;

data hpv_ext_24;
set hpv_ext_24;
format test_date yymmdd10. reg_date yymmdd10.;
sample_year=input(substr(x_sample_date,1,4),8.);
sample_month=input(substr(x_sample_date,6,2),8.);
sample_day=input(substr(x_sample_date,9,2),8.);
reg_year=input(substr(x_reg_date,1,4),8.);
reg_month=input(substr(x_reg_date,6,2),8.);
reg_day=input(substr(x_reg_date,9,2),8.);
test_date=MDY(sample_month,sample_day,sample_year);
reg_date=MDY(reg_month,reg_day,reg_year);
drop sample_year sample_month sample_day reg_year reg_month reg_day x_sample_date x_reg_date;;
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
*if hpvrisk=. then delete;
keep person_id test_date outcome  reg_date;
run; 
/*
NOTE: There were 7074224 observations read from the data set WORK.HPVALL.
NOTE: The data set WORK.HPVTEST has 7070014 observations and 4 variables.

*/

data cytotest;
set V_ncsr.nkc_trans_cell_2024;
if snomed_severity>5 then outcome=1;
else outcome=0;
*if snomed_severity<=3 then delete;
keep person_id  outcome x_reg_date x_sample_yr x_sample_date;
run;


data cytotest;
set cytotest;
format test_date yymmdd10. reg_date yymmdd10.;
sample_year=input(substr(x_sample_date,1,4),8.);
sample_month=input(substr(x_sample_date,6,2),8.);
sample_day=input(substr(x_sample_date,9,2),8.);
reg_year=input(substr(x_reg_date,1,4),8.);
reg_month=input(substr(x_reg_date,6,2),8.);
reg_day=input(substr(x_reg_date,9,2),8.);
test_date=MDY(sample_month,sample_day,sample_year);
reg_date=MDY(reg_month,reg_day,reg_year);
drop sample_year sample_month sample_day reg_year reg_month reg_day x_sample_date x_reg_date;
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
/*delete test with invalid sample date*/
data cytotest;
set cytotest;
if test_date=. and reg_date=. then delete;
if check=1 and change=. then delete;
keep person_id test_date outcome  reg_date;
run; 
/*NOTE: There were 25693461 observations read from the data set WORK.CYTOTEST.
NOTE: The data set WORK.CYTOTEST has 25690723 observations and 4 variables.

*/


data hpvtest;
set hpvtest;
testmethod=1;
run;
data cytotest;
set cytotest;
testmethod=0;
run;

/*keep one test per one day*/
data Alltest;
set cytotest hpvtest;
run;

proc sort data=Alltest;
by person_id test_date descending outcome descending testmethod;
run;
proc sort data=ALLtest nodupkey;
by person_id test_date;
run;
/*NOTE: There were 32760737 observations read from the data set WORK.ALLTEST.
NOTE: 3901744 observations with duplicate key values were deleted.
NOTE: The data set WORK.ALLTEST has 28858993 observations and 5 variables.
*/

/*match pop and test*/
data pop_test;
merge pop (in=a) Alltest;
by person_id;
if test_date<=first_hpv_date;
if a;
run;


data pop_test;
set pop_test;
test_age=floor((test_date-birth_date)/365.25);
run;

proc sort data=pop_test;
by person_id test_age;run;
/*only keep test start from age 20*/
data pop_test_20;
set pop_test;
if test_age>=20;
run;
/*NOTE: There were 7082020 observations read from the data set WORK.POP_TEST.
NOTE: The data set WORK.POP_TEST_20 has 6993145 observations and 8 variables.
*/
proc sort data=pop_test_20 nodupkey out=screen_start;
by person_id;
run;

data screen_start;
set screen_start;
keep person_id test_age;
rename test_age=screen_start_age;
run;

data pop_test_20;
merge pop_test_20 screen_start;
by person_id;
run;

proc sort data=pop_test_20;
by person_id descending test_date;
run;

data pop_test_20;
set pop_test_20;
format next_sample_date yymmdd10.;
by person_id;
next_sample_date=lag(test_date);
if first_hpv_date=test_date then next_sample_date=test_date;
run;


/*policy with cytology before 50 3 years interval after 50 5 years interval
hpv before 50 5 years interval after 50 7 years interval
if the outcome=1 then all within 1 year interval*/
/*next sample date plan*/
/*UPDATE IN 20250121 as long as they have a test nomatter the result
under 50 3 years; over 50 5 years 
for HPV undr 50 3 years; over 50 7 years*/

data pop_test_plan;
set pop_test_20;
format plan_test_date yymmdd10.;
if test_age<50 and testmethod=0 then plan_test_date=test_date+365.25*3;
else if test_age>=50 and testmethod=0 then plan_test_date=test_date+365.25*5;
if test_age<50 and testmethod=1 then plan_test_date=test_date+365.25*3;
else if test_age>=50 and testmethod=1 then plan_test_date=test_date+365.25*7;
if test_age>=70 then do; next_sample_date=test_date; plan_test_date=test_date; end;
*if outcome=1 then plan_test_date=test_date+365.25;
run;




/*time covered with test, t_within
time without test,t_without*/
data time;
set pop_test_plan;
by person_id;
if next_sample_date<=plan_test_date then t_within=next_sample_date-test_date;
else  t_within=plan_test_date-test_date;
run;

/*add count*/

proc sql;
create table time_period  as
select  person_id,count(person_id)as count, sum(t_within) as within_all, max(test_date) as last_date, min(test_date) as first_date, first_hpv_date,
birth_date, screen_start_age from time
group by person_id;
quit;

data time_period;
set time_period;
format last_date yymmdd10. first_date yymmdd10.;
run;


proc sort data=time_period nodupkey;
by person_id;
run;


/*start time*/
data start_time;
set time_period;
format start_time yymmdd10.;
if screen_start_age>23 then theo_screening_time=(year(last_date)-year(birth_date)-23)*365.25;
else theo_screening_time=last_date-first_date;
run;
data proportion;
set start_time;
index=within_all/(theo_screening_time)*100;
run;

data proportion;
set proportion;
birth_yr=year(birth_date);
if birth_yr<=1948 then birthco=1;
else if 1948<birth_yr<=1953 then birthco=2; /*age75-71 in 2024*/
else if 1953<birth_yr<=1958 then birthco=3; /*age70-66 in 2024*/
else if 1958<birth_yr<=1963 then birthco=4; /*age65-61 in 2024*/
else if 1963<birth_yr<=1968 then birthco=5; /*age60-56 in 2024*/
else birthco=6; /* younger than 55 in 2024*/
run;
data proportion;
set proportion;
if index=0 then index=.;
run;
proc rank data=proportion out=index_rank groups=3;
var index;
ranks ranked_index;
run;

proc freq data=index_rank;
tables birthco*ranked_index/missing;
run;

proc summary data=index_rank;
var index;
class ranked_index;
output out=cutoff max=/autoname;
run;
/*55,73,100*/

data index_rank;
set index_rank;
num_pre_test=count-1;
if num_pre_test=0 then noscr_his=1;else noscr_his=0;
keep person_id num_pre_test birthco ranked_index screen_start_age index noscr_his;
run;

data HPV_cohort;
set preabn.NKCx_nohys_cohort_260202;
if previous_pos=. then previous_pos=0;
if previous_pos_cat=. then previous_pos_cat=0;
run;

data HPV_cohort;
merge HPV_cohort index_rank;
by person_id;
run;

proc freq data=HPV_cohort;
tables noscr_his*previous_pos;
run;

data preabn.NKCx_nohys_cohort_260202;
set HPV_cohort;
run;
