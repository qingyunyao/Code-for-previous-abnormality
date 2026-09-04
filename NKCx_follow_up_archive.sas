/*File name..: identify the longterm follow-up of the CIN2+ and ICC NKCx*/
/*Study......: Phd project study 2 the HPV positivity after one HPV test*/
/*Author.....: Qingyun Yao*/
/*Date.......: 2024/0623*/
/*Updated....: */
/*Purpose....: /*identify the end point of CIN2+ and ICC after first HPV */
/*The end of follow-up of CIN2+ is the last registered test date in NKCx or the detection of CIN2+
Then end of follow-up of ICC is the end of the cancer registry

/*
libname NCSR  odbc complete="server=meb-sql02.meb.ki.se;driver=SQL Server Native Client 11.0;Trusted_Connection=Yes;database=NCSR" schema=NCSR ;
libname V_ncsr     odbc complete="server=meb-sql02.meb.ki.se;driver=SQL Server Native Client 11.0;Trusted_Connection=Yes;database=NCSR" schema=V_ncsr ;
libname PreAbn 'P:\ACCES\ACCES_Research\Qingyun\previous abnormality\Data';
libname SOCMOB4  odbc complete="dsn=kosmos;database=CERVIX_SOCMOB4" schema=CLEAN_2022 ;


				*/
/*Note.......: */
*------------------------------------------------------------------------;
/* Data used...: socmob4.nkc_hpv socmob4.nkc_trans_cell  self_samp_pop;

/* Data created.: pop_index_&date */
proc datasets library=work kill;quit;

libname NCSR odbc noprompt="dsn=kosmos; database=NCSR" schema=NCSR;
libname V_NCSR odbc noprompt="dsn=kosmos; database=NCSR" schema=V_NCSR ;
libname PreAbn 'P:\ACCES\ACCES_Research\Qingyun\previous abnormality\Data';

/*get the population and the end point of CIN2+ and ICC*/
data pop;
set preabn.NKCx_nohys_cohort_260202;
keep person_id HPVrisk first_hpv_date previous_pos previous_pos_cat noscr_his hpv_age laboratory_id hpv_year birth_date selftest birthco;
run;


data pop;
set pop;
if laboratory_id in ('088' '999') then region='Stockholm-Gotland';
else if laboratory_id in ('211' '231' '237' '251') then region='Southeast Sweden';
else if laboratory_id in ('241' '271' '411' '417') then region='South Sweden';
else if laboratory_id in ('621' '631' '641' '651') then region='North Sweden';
else if laboratory_id in ('421' '427' '501' '507' '511' '517' '521' '527' '531' '537') then region='West Sweden';
else if laboratory_id in ('121' '127' '131' '541' '551' '561' '567' '571' '577' '611') then region='Middle Sweden';
if 50<=hpv_age<=59 then age_group=1;
else if 60<=hpv_age<=70 then age_group=2;
/*50-54, 55-59, 60-64, 65-70*/
if hpv_year<=2016 then sample_period=1;
else if hpv_year<=2019 then sample_period=2;
else sample_period=3;
if hpv_year<=2019 then calender_period=1;
else calender_period=2;
run;
/*find all the test after first analysis and keep the last_test_date*/
data cyto;
set v_ncsr.nkc_trans_cell_2024;
keep person_id x_sample_date snomed_severity  x_reg_date x_sample_yr;
run;
data hpv;
set v_ncsr.nkc_hpv_2024;
keep person_id x_sample_date HPVDIAG  x_reg_date x_sample_yr;
if selftest=1 and x_reg_date^='' then x_sample_date=x_reg_date;
run;
data hpv_ext;
set v_NCSR.nkc_ext_hpv_2024;
keep person_id x_sample_date HPVDIAG  x_reg_date x_sample_yr;
run;
data pad;
set v_ncsr.nkc_pad_translated_2024;
keep person_id pad_sev x_sample_date;
where TOPO3='T83';
run;

/*reorganize the sample date*/
data cyto;
set cyto;
format sample_date yymmdd10. reg_date yymmdd10.;
sample_year=input(substr(x_sample_date,1,4),8.);
sample_month=input(substr(x_sample_date,6,2),8.);
sample_day=input(substr(x_sample_date,9,2),8.);
sample_date=MDY(sample_month,sample_day,sample_year);
reg_year=input(substr(x_reg_date,1,4),8.);
reg_month=input(substr(x_reg_date,6,2),8.);
reg_day=input(substr(x_reg_date,9,2),8.);
reg_date=MDY(reg_month,reg_day,reg_year);
drop sample_year sample_month sample_day reg_year reg_month reg_day x_sample_date x_reg_date;
if snomed_severity>5 then outcome=1;
if snomed_severity<=3 then delete;
else outcome=0;
run;
/*check time */

data cyto;
set cyto;
test_year=year(sample_date);
reg_year=year(reg_date);
if abs(test_year-x_sample_yr)>1 then check=1;
if check=1 and reg_year=x_sample_yr then do; sample_date=reg_date; change=1;end;
run;
proc freq data=cyto;
tables check*change/missing;
run;
/*delete test with invalid sample date*/
data cyto;
set cyto;
if sample_date=. and reg_date=. then delete;
if check=1 and change=. then delete;
keep person_id sample_date outcome snomed_severity;
run; 
/*else*/
data hpv;
set hpv;
format sample_date yymmdd10. reg_date yymmdd10.;
sample_year=input(substr(x_sample_date,1,4),8.);
sample_month=input(substr(x_sample_date,6,2),8.);
sample_day=input(substr(x_sample_date,9,2),8.);
sample_date=MDY(sample_month,sample_day,sample_year);
reg_year=input(substr(x_reg_date,1,4),8.);
reg_month=input(substr(x_reg_date,6,2),8.);
reg_day=input(substr(x_reg_date,9,2),8.);
reg_date=MDY(reg_month,reg_day,reg_year);
drop sample_year sample_month sample_day reg_year reg_month reg_day x_sample_date x_reg_date;
if HPVDIAG='POS' then outcome=1;
else if HPVDIAG='NEG' then outcome=0;
run;
data hpv_ext;
set hpv_ext;
format sample_date yymmdd10. reg_date yymmdd10.;
sample_year=input(substr(x_sample_date,1,4),8.);
sample_month=input(substr(x_sample_date,6,2),8.);
sample_day=input(substr(x_sample_date,9,2),8.);
sample_date=MDY(sample_month,sample_day,sample_year);
reg_year=input(substr(x_reg_date,1,4),8.);
reg_month=input(substr(x_reg_date,6,2),8.);
reg_day=input(substr(x_reg_date,9,2),8.);
reg_date=MDY(reg_month,reg_day,reg_year);
drop sample_year sample_month sample_day reg_year reg_month reg_day x_sample_date x_reg_date;
if HPVDIAG='POS' then outcome=1;
else if HPVDIAG='NEG' then outcome=0;
run;

data pad;
set pad;
format sample_date yymmdd10.;
sample_year=input(substr(x_sample_date,1,4),8.);
sample_month=input(substr(x_sample_date,6,2),8.);
sample_day=input(substr(x_sample_date,9,2),8.);
sample_date=MDY(sample_month,sample_day,sample_year);
drop sample_year sample_month sample_day x_sample_date;
if pad_sev>=5 then pad_cin2=1;
run;
/*check date*/

data hpvall;
set hpv hpv_ext;
run;

data hpvall;
set hpvall;
test_year=year(sample_date);
reg_year=year(reg_date);
if abs(test_year-x_sample_yr)>1 then check=1;
if check=1 and reg_year=x_sample_yr then do; sample_date=reg_date; change=1;end;
run;

data hpvtest;
set hpvall;
if sample_date=. and reg_date=. then delete;
if check=1 and change=. then delete;
keep person_id sample_date outcome HPVDIAG;
run; 
/*get rid of the invalid test*/
/*cyto <4 delete, HPV outcome invalid delete*/
/*exclude invalid test*/
data cyto_valid;
set cyto;
if snomed_severity<=3 then delete;
run;
data hpv_valid;
set hpvtest;
if HPVDIAG not in ('POS' 'NEG') then delete;
run;

data pad_valid;
set pad;
if pad_sev=0 then delete;
run;
data test_date;
set cyto_valid hpv_valid  pad_valid;
keep person_id sample_date;
run;
proc sort data=test_date nodupkey;
by person_id descending sample_date;
run;
proc sort data=test_date nodupkey;
by person_id;
run;
data last_test_date;
set test_date;
rename sample_date=last_test_date;
run;
proc sort data=pop;
by person_id;
run;
proc sort data=last_test_date;
by person_id;
run;
data pop_last_date;
merge pop (in=a) last_test_date;
by person_id;
if a;
run;

/*NOTE: There were 742085 observations read from the data set WORK.POP.
NOTE: There were 4704878 observations read from the data set WORK.LAST_TEST_DATE.
NOTE: The data set WORK.POP_LAST_DATE has 742085 observations and 17 variables.
*/
/*find if they have a CIN2+ after first_hpv_date in pad*/
data pad_cin2;
set pad;
if pad_cin2>0;
drop x_sample_date;
rename sample_date=cin2_date;
run;
proc sort data=pad_cin2;
by person_id cin2_date;
run;
proc sort data=pop;
by person_id;
run;
data pop_cin2;
merge pop (in=a) pad_cin2;
by person_id;
if cin2_date>=first_hpv_date;
if a;
run;
/*keep the earliest data of CIN2*/
proc sort data=pop_cin2;
by person_id cin2_date;
run;
proc sort data=pop_cin2 nodupkey;
by person_id;
run;

data pop_cin2;
set pop_cin2;
fu_cin2=1;
if cin2_date<=first_hpv_date+365.25 then fu_cin2_12=1;
rename pad_sev=fu_cin2_sev cin2_date=fu_cin2_date;
keep person_id pad_sev cin2_date fu_cin2_12 fu_cin2;
run;

/*NOTE: There were 4917 observations read from the data set WORK.POP_CIN2.
NOTE: The data set WORK.POP_CIN2 has 4917 observations and 5 variables.
*/
/*find out if they have histopathology test after first HPV test */
proc sort data=pad_valid;
by person_id;
run;
data pop_pad;
merge pop (in=a) pad_valid;
by person_id;
if sample_date>=first_hpv_date;
if a;
run;
data pop_pad;
set pop_pad;
drop x_sample_date pad_cin2;
rename sample_date=pad_date;
run;
/*keep first pad*/
proc sort data=pop_pad;
by person_id pad_date;
run;
proc sort data=pop_pad nodupkey;
by person_id;
run;

data pop_pad;
set pop_pad;
if pad_date<=first_hpv_date+365.25 then pad_12=1;
run;
data pop_pad;
set pop_pad;
fu_pad=1;
rename pad_sev=fu_pad_sev pad_date=fu_pad_date pad_12=fu_pad_12;
keep person_id pad_sev pad_date pad_12 fu_pad;
run;

/*find cin2+ detection rate within 12 months*/
proc freq data=pop_pad;
tables fu_pad_12;
run;


/*find if they have a cancer after first_hpv_date in QGCR register*/

data cancer;
set v_ncsr.gyn_cancer_2024;
keep person_id cxca_date figo_stage histo_type;
run;


data cancer;
set cancer;
format diag_date yymmdd10.;
diag_year=input(substr(cxca_date,1,4),8.);
diag_month=input(substr(cxca_date,6,2),8.);
diag_day=input(substr(cxca_date,9,2),8.);
diag_date=MDY(diag_month, diag_day, diag_year);
run;
/*need to correct cancer date with pad date*/
data pad;
set v_ncsr.nkc_pad_translated_2024;
keep person_id PAD_sev x_sample_date;
if pad_sev>=5;
run;
data pad;
set pad;
format pad_date yymmdd10.;
sample_year=input(substr(x_sample_date,1,4),8.);
sample_month=input(substr(x_sample_date,6,2),8.);
sample_day=input(substr(x_sample_date,9,2),8.);
pad_date=MDY(sample_month,sample_day,sample_year);
drop x_sample_date;
run;
proc sort data=pad;
by person_id pad_date;
run;
proc sort data=pad nodupkey;
by person_id;
run;

proc sort data=cancer;
by person_id;
run;
data cancer_date;
merge cancer (in=a) pad;
by person_id;
if a;
run;

data cancer_date;
set cancer_date;
if diag_date<pad_date<diag_date+365.25 then cxca_date=pad_date;
run;

data cancer;
set cancer_date;
keep person_id diag_date cancer figo_stage histo_type;
cancer=1;
run;
proc sort data=cancer nodupkey;
by person_id diag_date;
run;

/*match the cancer found after first HPV test*/
data pop_cancer;
merge pop (in=a) cancer;
by person_id;
if diag_date>= first_hpv_date;
if a;
run;

proc sort data=pop_cancer;
by person_id diag_date;
run;
proc sort data=pop_cancer nodupkey;
by person_id ;
run;

data pop_cancer;
set pop_cancer;
keep person_id diag_date cancer figo_stage histo_type;
rename diag_date=cancer_date;
run;
/*merge everything*/
data pop_longterm_outcome;
merge pop pop_last_date pop_pad pop_cin2 pop_cancer ;
by person_id;
run;
/*
NOTE: There were 742085 observations read from the data set WORK.POP.
NOTE: There were 742085 observations read from the data set WORK.POP_LAST_DATE.
NOTE: There were 31943 observations read from the data set WORK.POP_PAD.
NOTE: There were 4917 observations read from the data set WORK.POP_CIN2.
NOTE: There were 408 observations read from the data set WORK.POP_CANCER.
NOTE: The data set WORK.POP_LONGTERM_OUTCOME has 742085 observations and 29 variables.

*/


/**/



/*deregistration*/
/*deregistration data*/

data dereg;
set v_ncsr.nkc_pop_2024;
keep person_id dereg_date dereg_reas;
run;

data dereg;
set dereg;
format dereg_date_final yymmdd10.;
dereg_year=input(substr(dereg_date,1,4),8.);
dereg_month=input(substr(dereg_date,5,2),8.);
dereg_day=input(substr(dereg_date,7,2),8.);
if dereg_month=0 then do;dereg_month=7;dereg_day=1;end;
else if dereg_day=0 then dereg_day=15;
if dereg_month=4 and dereg_day=31 then dereg_day=30;
dereg_date_final=MDY(dereg_month, dereg_day, dereg_year);
keep person_id  dereg_date_final dereg_reas;
where dereg_date^='';
run;
proc freq data=dereg;
tables dereg_reas;
run;

data Hys_nkc_2024;
set v_ncsr.nkc_deregister_2024;
keep person_id x_dereg_reason x_dereg_from_date;
run;


data hys_nkc;
set Hys_nkc_2024 ;
format dereg_date yymmdd10.;
dereg_year=input(substr(x_dereg_from_date,1,4),8.);
dereg_month=input(substr(x_dereg_from_date,6,2),8.);
dereg_day=input(substr(x_dereg_from_date,9,2),8.);
dereg_date=MDY(dereg_month, dereg_day, dereg_year);
drop x_dereg_from_date dereg_year dereg_month dereg_day;
rename x_dereg_reason=dereg_reas ;
run;

proc sort data=hys_nkc nodupkey;
by person_id dereg_date;
run;
proc freq data=hys_nkc;
tables dereg_reas;
run;


/*end of follow up for CIN2+ is the last test they have, 
end of cancer is cancer date, total hysterectomi date (both from hysterectomi and nkcx self report) and deregistration date(from pop registry move out or death)*/
proc sort data=dereg;
by person_id;
run;
data pop_dereg;
merge pop(in=a) dereg;
by person_id;
if dereg_date_final>=first_hpv_date;
if a;
run;
/*
NOTE: There were 742085 observations read from the data set WORK.POP.
NOTE: There were 2474403 observations read from the data set WORK.DEREG.
NOTE: The data set WORK.POP_DEREG has 13597 observations and 18 variables.
*/
proc sort data=pop_dereg;
by person_id dereg_date_final;
run;
proc sort data=pop_dereg nodupkey;
by person_id;
run;
/*
NOTE: There were 13597 observations read from the data set WORK.POP_DEREG.
NOTE: 0 observations with duplicate key values were deleted.
NOTE: The data set WORK.POP_DEREG has 13597 observations and 18 variables.
*/

proc sort data=hys_nkc;
by person_id;
run;
data pop_hyster_nkc;
merge pop(in=a) hys_nkc;
by person_id;
if dereg_date>=first_hpv_date;
if a;
if dereg_reas='Hysterektomi';
run;

proc sort data=pop_hyster_nkc;
by person_id dereg_date;
run;
proc sort data=pop_hyster_nkc nodupkey;
by person_id;
run;



/*organize deregistration information*/
data pop_dereg;
set pop_dereg;
keep person_id dereg_reas dereg_date_final pop_dereg;
rename dereg_reas=dereg_pop;
pop_dereg=1;
run;
data pop_hyster_nkc;
set pop_hyster_nkc;
keep person_id dereg_reas dereg_date nkc_hys;
rename dereg_reas=dereg_resa_nkc;
nkc_hys=1;
run;


data pop_dereg_all;
merge pop_dereg pop_hyster_nkc ;
by person_id;
run;
/*
NOTE: There were 13597 observations read from the data set WORK.POP_DEREG.
NOTE: There were 3557 observations read from the data set WORK.POP_HYSTER_NKC.
NOTE: The data set WORK.POP_DEREG_ALL has 17075 observations and 7 variables.
*/
/*organize deregistration data*/
data pop_dereg_final;
set pop_dereg_all;
format d_date yymmdd10.;
d_date=min(of dereg_date_final dereg_date );
run;

data pop_dereg_final;
set pop_dereg_final;
format d_reas $20.;
if d_date=dereg_date_final then d_reas=dereg_pop;
else d_reas='Hysterektomi';
run;

data pop_dereg_final;
set pop_dereg_final;
keep person_id d_date d_reas;
run;


/*keep all the information and compare later*/

data pop_longterm_outcome;
merge pop  pop_last_date pop_pad pop_cin2 pop_cancer pop_dereg_final;
by person_id;
run;
/*
NOTE: There were 742085 observations read from the data set WORK.POP.
NOTE: There were 742085 observations read from the data set WORK.POP_LAST_DATE.
NOTE: There were 31943 observations read from the data set WORK.POP_PAD.
NOTE: There were 4917 observations read from the data set WORK.POP_CIN2.
NOTE: There were 408 observations read from the data set WORK.POP_CANCER.
NOTE: There were 17075 observations read from the data set WORK.POP_DEREG_FINAL.
NOTE: The data set WORK.POP_LONGTERM_OUTCOME has 742085 observations and 31 variables.
*/

/*calculate per
son-times
1. CIN2+ ; if no cin2+ end of follow_up last test otherwise cin2_date
2. cancer: if no cancer end of follow_up 20221231/hyster_date or dereg_date otherwise cancer_date*/

data pop_sur_time;
set pop_longterm_outcome;
format eof_cin2 yymmdd10. eof_cancer yymmdd10.;
if fu_cin2=1 then eof_cin2=fu_cin2_date;
else eof_cin2=last_test_date;
if cancer=1 then eof_cancer=cancer_date;
else eof_cancer='31Dec2024'd;
if d_reas^=''  then eof_cancer= min(of eof_cancer d_date);
if first_hpv_date>'31Dec2024'd then eof_cancer=.;
sur_cin2=(eof_cin2-first_hpv_date)/365.25;
sur_cancer=(eof_cancer-first_hpv_date)/365.25;
keep person_id eof_cin2 eof_cancer sur_cin2 sur_cancer ;
run;
/*merge all*/
data follow_up_all;
merge pop_longterm_outcome pop_sur_time;
by person_id;
run;

data follow_up_all;
set follow_up_all;
if fu_pad=. then fu_pad=0;
if fu_cin2=. then fu_cin2=0;
if fu_cin2_12=. then fu_cin2_12=0;
if cancer=. then cancer=0;
run;


/*the information of follow up other test since the first HPV test using code follow_up_test for firstHPV*/

/*save the dataset*/
data preabn.NKCx_nohys_longterm_fu;
set follow_up_all;
run;

