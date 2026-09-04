/*File name..: building cohort for preliminary outcome (HPV positivity) use NKCx */
/*Study......: Phd project study 2 HPV-based screening among women with previous abnormality*/
/*Author.....: Qingyun Yao*/
/*Date.......: 2025/03/13*/
/*Updated....:update from the code written for Access dataset/lack cancer information
big update 20250627 new data from 2024*/
/*Purpose....: build up cohort for final analysis
20250707 make lsil and ascus one subgroup based on the clinical procedure*/
/*Note.......: population defination: women's first valid HPV-based screening test from 2012 to the end of the date registered in the
registry (use hpv data). Women aged 50-70 at the time of their first valid HPV-based screening. Women without any screening test for the past 4 years (both ext and hpv, cell). 
Women had no hysterectomy before the first screening test (total hysterectomy ).
test age was calculated by test_year-sample_year
*/
*------------------------------------------------------------------------;
/* Data used...:V_nscr.nkc_hpv_2024.dta
V_nscr.nkc_trans_cell_2024.dta
V_nscr.nkc_ext_hpv_2024.dta
V_nscr.nkc_pad_translated_2024.dta
V_nscr.nkc_person _2024.dta
V_nscr.nkc_deregister_2024.dta

libname NCSR  odbc complete="server=meb-sql02.meb.ki.se;driver=SQL Server Native Client 11.0;Trusted_Connection=Yes;database=NCSR" schema=NCSR ;
libname V_NCSR     odbc complete="server=meb-sql02.meb.ki.se;driver=SQL Server Native Client 11.0;Trusted_Connection=Yes;database=NCSR" schema=V_NCSR ;
libname PreAbn 'P:\ACCES\ACCES_Research\Qingyun\previous abnormality\Data';
libname S4 odbc noprompt="dsn=kosmos;database=Cervix_socmob4" schema=clean_2022;
libname Socmob4 odbc noprompt="dsn=kosmos;database=Cervix_socmob4" schema=clean_2023_2024;

There is one systematic mistake in ext_hpv in lab 411, x_sample_year=2016, x_sample_date were in 2002 delete it when analyzing

/* Data created.:  HPVpositivity_20241008.dta*/
/*sas version.: SAS9.4*/
/*main program*/
libname NCSR odbc noprompt="dsn=kosmos; database=NCSR" schema=NCSR;
libname V_NCSR odbc noprompt="dsn=kosmos; database=NCSR" schema=V_NCSR ;
libname PreAbn 'P:\ACCES\ACCES_Research\Qingyun\previous abnormality\Data';

/*indicator for primary HPV screening test*/
%let indicator='P:\ACCES\ACCES_Research\Qingyun\previous abnormality\Program\NKCx_remissserier_indikation_QY.xlsx'; 

proc datasets library=work kill;quit;
/*identify all women with primary HPV testing in each county*/

data HPV_24;
set V_ncsr.NKC_hpv_2024;
keep person_id age x_sample_date x_sample_yr x_reg_date county_id scr_type HPVDIAG laboratory_id referral_nr x_referral_type sample_type hpv_type selftest;
run;

/*clean the data*/
data HPV_24;
set HPV_24;
format sample_date yymmdd10. reg_date yymmdd10.;
sample_year=input(substr(x_sample_date,1,4),8.);
sample_month=input(substr(x_sample_date,6,2),8.);
sample_day=input(substr(x_sample_date,9,2),8.);
reg_year=input(substr(x_reg_date,1,4),8.);
reg_month=input(substr(x_reg_date,6,2),8.);
reg_day=input(substr(x_reg_date,9,2),8.);
sample_date=MDY(sample_month,sample_day,sample_year);
reg_date=MDY(reg_month,reg_day,reg_year);
drop sample_year sample_month sample_day reg_year reg_month reg_day;
if selftest=1  and x_reg_date^='' then sample_date=reg_date;
run;

/*check if there is systemmatic problem with sample_date*/

data hpv_check;
set HPV_24;
format test_date yymmdd10.;
test_year=year(sample_date);
reg_year=year(reg_date);
if abs(test_year-x_sample_yr)>1 then check=1;
if check=1 and reg_year=x_sample_yr then do; test_date=reg_date; change=1;end;
run;

proc freq data=hpv_check;
tables change change*check/missing;
run;
/*has 115 tests the sample year and register year are different; and 15 in sample_date;
the effect of wrong sample date will be small*/

data hpv_24;
set hpv_check;
if test_date^=. then sample_date=test_date;
test_year=year(sample_date);
drop  reg_year x_sample_yr check change test_date x_sample_date x_reg_date;
run;

data HPV_24;
set hpv_24;
if test_year<=2024;
run;
/*NOTE: There were 6189660 observations read from the data set WORK.HPV_24.
NOTE: The data set WORK.HPV_24 has 6179989 observations and 14 variables.
*/


/*find primay HPV test*/
/*check the referral type*/
proc freq data=hpv_24;
tables x_referral_type; 
run;
data hpv_24;
set hpv_24;
if x_referral_type='SY_G' then x_referral_type='SY-G';
if x_referral_type='ES_G' then x_referral_type='ES-G';
if x_referral_type='VG_G' then x_referral_type='VG-G';
if x_referral_type='VG_V' then x_referral_type='VG-V';
if laboratory_id in ('501' '507' '517' '527') then do;
referral_type_1=substr(x_referral_type,1,1);
referral_type_2=substr(x_referral_type,2,1);
referral_type_3=substr(x_referral_type,3,1);
if referral_type_2 IN ('1' '2' '3' '4' '5' '6' '7' '8' '9' '0') then do; x_referral_type=substr(x_referral_type,1,1);end;
else if referral_type_3 IN ('1' '2' '3' '4' '5' '6' '7' '8' '9' '0') then do;x_referral_type=substr(x_referral_type,1,2);end;
else x_referral_type=substr(x_referral_type,1,3);
drop referral_type_1 referral_type_2 referral_type_3;
end;
run;
/*link to population registry for birthday*/
data pop_valid;
set V_ncsr.nkc_person_2024;
keep person_id birth_date;
where valid_pnr=1;
run;
data dup;
set pop_valid;
by person_id;
if first.person_id and last.person_id then delete;
run;
/*noduplicate in person_id*/
data pop_valid;
set pop_valid;
by person_id;
if first.person_id and last.person_id ;
run;
/*
NOTE: There were 5252650 observations read from the data set WORK.POP_VALID.
NOTE: The data set WORK.POP_VALID has 5252650 observations and 2 variables.

*/

data pop_valid;
set pop_valid;
format b_date yymmdd10.;
b_year=input(substr(birth_date,1,4),8.);
b_month=input(substr(birth_date,6,2),8.);
b_day=input(substr(birth_date,9,2),8.);
b_date=MDY(b_month, b_day, b_year);
drop b_year b_month b_day;
run;

data pop_valid;
set pop_valid;
drop birth_date;
rename b_date=birth_date;
where b_date^=.;
run;

proc sort data=hpv_24;
by person_id;
run;

data hpv_24;
merge hpv_24 (in=a) pop_valid;
by person_id;
if a;
run;

/*calculate age*/
data hpv_24;
set hpv_24;
hpv_age=floor((sample_date-birth_date)/365.25);
run;

/*organize the HPV test with only most severe test per day*/
/*if a women have both 16 and HPV high risk types in one sample_date it means she is 16 positive
if a women has other HPV high risk not 16/18 positive and HPV high risk types it means she is other HrPositive*/

proc freq data=hpv_24;
tables hpv_type;
run;
data hpv_24_risk;
set hpv_24;
if hpv_type='16' then HPV16=1;
if hpv_type='18' then HPV18=1;
hpv16_18=index(hpv_type,'16 18')>0;
HPV18_45=index(hpv_type,'18 45')>0;
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
if hpv_type=:'Other' then HPVHr=1;
if hpv_type=:'HPV DNA' then HPVmissing=1;
if hpv_type=:'HPV high risk type' then HPVHr=1;
if hpv_type=:'HPV type unknown' then HPVmissing=1;
if hpv_type=:'HPV type missing' then HPVHr=1;
if hpv_type='Not available' then HPVmissing=1;
if hpv_type='Irrelevant code' then HPVmissing=1;
if HPV18_45=1 then HPV45=0;
if hpv16_18=1 then HPV16=1;
if HPV33=1 or HPV58=1 or HPV45=1 or HPV31=1 or HPV52=1 or HPV35=1 or HPV39=1 or HPV51=1 or HPV59=1 or HPV56=1 or HPV68=1 or HPV66=1 then
HPVHr=1;
run;
/*hpvtype missing but positive count as HrHPV positive, information from JW*/

proc freq data=hpv_24_risk;
tables HPV16 HPV18 HPV18_45 HPV45 HPVHr HPVmissing/missing;
run;


data hpv_24_risk;
set hpv_24_risk;
if HPV16=1 and HPVDIAG='POS' then do; HPVrisk=4;Type16=1;end;
else if HPV18=1 and HPVDIAG='POS' then do; HPVrisk=3;Type18=1;end;
else if HPV18_45=1 and HPVDIAG='POS' then do;HPVrisk=3; Type18_45=1;end;
else if HPVHr=1 and HPVDIAG='POS' then do; HPVrisk=2; TypeHr=1;end;
else if HPVmissing=1 and HPVDIAG='POS' then HPVrisk=.;
else if HPVDIAG='POS' then HPVrisk=0;/*HPV low risk count as HPV negative*/
else if HPVDIAG^='NEG' then HPVrisk=.;
else if HPVDIAG='NEG' then HPVrisk=0;
drop HPV16 HPV18  HPV33 HPV58 HPV31 HPV52 HPV35 HPV39 HPV51 HPV59 HPV56 HPV68 HPV66;
run;

proc freq data=hpv_24_risk;
tables HPVrisk HPVrisk*HPVDIAG/missing;
run;
/*one woman aged 36 with a HPV 34 positive  the risk is marked as invalid, not going to affect my population no correction is done here*/

proc freq data=hpv_24_risk;
tables hpv_type/missing;
where HPVrisk=0 and HPVDIAG='POS';
run;

proc sort data=hpv_24_risk;
by person_id sample_date descending HPVrisk;
run;

/*find out the multiple infections just in case we need this*/
/*add multi infection index*/
data multi_inf;
set hpv_24_risk;
keep person_id sample_date HPVrisk outcome hpv_type;
if HPVrisk>1 then outcome=1;
else outcome=0;
run;
/*multi-infection only exist in positive women*/
data multi_inf;
set multi_inf;
if outcome=1;
run;
/*delete same hpv_type per day*/
proc sort data=multi_inf nodupkey;
by person_id sample_date hpv_type descending HPVrisk;
run;
/*NOTE: There were 1465650 observations read from the data set WORK.MULTI_INF.
NOTE: 148148 observations with duplicate key values were deleted.
NOTE: The data set WORK.MULTI_INF has 1317502 observations and 5 variables.
*/
proc sort data=multi_inf;
by person_id sample_date HPVrisk hpv_type;
run;
/*try another way of coding*/
/*outcome with hpv high risk or hpv type missing represent early stage HPV screening ,later the type information is updated, so these
two don't count as a individual type outcome*/
data multi_inf_test;
set multi_inf;
if hpv_type in ('HPV high risk types' 'HPV type missing') then delete;
run;
data multi_inf_mark;
set multi_inf_test;
format b_date yymmdd10.;
retain b_type b_date id;
if person_id=id and sample_date=b_date and hpv_type^=b_type then multi_inf=1;
b_type=hpv_type;
b_date=sample_date;
id=person_id;
run;

proc sort data= multi_inf_mark;
by person_id sample_date descending multi_inf;
run;
proc sort data=multi_inf_mark nodupkey;
by person_id sample_date ;
run;
proc freq data=multi_inf_mark;
tables multi_inf;
run;
/*test*/
data multioutcomes;
set multi_inf_mark;
keep person_id sample_date multi_inf;
if multi_inf=. then multi_inf=0;
run;
/*merge with hpv risk information*/

proc sort data=hpv_24_risk;
by person_id sample_date;
proc sort data=multioutcomes;
by person_id sample_date;
run;


proc sort data=hpv_24_risk;
by person_id sample_date descending HPVrisk;
run;

proc sort data=hpv_24_risk nodupkey;
by person_id sample_date ;
run;
/*
NOTE: There were 6179989 observations read from the data set WORK.HPV_24_RISK.
NOTE: 1312692 observations with duplicate key values were deleted.
NOTE: The data set WORK.HPV_24_RISK has 4867297 observations and 26 variables.

*/

data hpv_24_risk;
merge hpv_24_risk multioutcomes;
by person_id sample_date;
run;
/*
NOTE: There were 4867297 observations read from the data set WORK.HPV_24_RISK.
NOTE: There were 751549 observations read from the data set WORK.MULTIOUTCOMES.
NOTE: The data set WORK.HPV_24_RISK has 4867297 observations and 27 variables.
*/



/*use the new indicator that has been updated*/

%let path='P:\ACCES\ACCES_Research\Qingyun\previous abnormality\Program\NKCx_remissserier_indikation_QY.xlsx'; /*'P:\GynCell\Gyncell_Docs\Rapporter\Rapport 2019\Nya rapporter\NKCx_remissserier_indikation_210827.xlsx';/*'P:\GynCell\Gyncell_Docs\Rapporter\Rapport 2019\Nya rapporter\NKCx_remissserier_indikation.xlsx'; */
proc import out = remiss 
			datafile= &path 
            dbms=xlsx replace;
     	    getnames=yes;
run;

data remiss1;
set remiss;
drop laboratory_id_2;
run;

data remiss2;
set remiss;
where laboratory_id_2 ne '';
drop laboratory_id;
run;

data remiss2;
set remiss2;
length laboratory_id $3.;
laboratory_id=put(laboratory_id_2, 3.);
drop laboratory_id_2;
run;
 
proc append base=remiss1 data=remiss2;
run;

/*not gonne include people with double testing in this analysis*/

data remiss1;
set remiss1;
code+1;
lab_id=input(laboratory_id,8.);
rename sample_type=r_sample_type x_referral_type=r_referral_type;
run;


proc sort data=remiss1 nodupkey;
by laboratory_id r_referral_type r_sample_type in_use in_use_from in_use_until below_30_CYT  scr_type_in_NKCx scr_type_matter ;
run;

data remiss_st remiss_nst remiss_nrt;
set remiss1;
if r_sample_type='' then output remiss_nst;
else if r_referral_type='' then output remiss_nrt;
else output remiss_st;
keep laboratory_id r_referral_type r_sample_type in_use_from in_use_until scr_type_in_NKCx below_30_cyt scr_type_matter screeening category primary_analysis lab_id code;
run;
/*only remiss_st and remiss_nst gonna put into use*/

/*match hpv and indicator*/
proc sql;
create table hpv_nst as select a.*, b.*  
from hpv_24_risk as a left join remiss_nst as b on a.laboratory_id=b.laboratory_id and a.x_referral_type=b.r_referral_type  ;
quit;

data hpv_nst;
set hpv_nst;
if in_use_from^=. and sample_date<in_use_from then delete;
else if in_use_until^=. and sample_date>in_use_until then delete;
rename code=nst_code;
if below_30_cyt=1 and hpv_age>=30 then delete;
if below_30_cyt=0 and hpv_age<30 then delete;
if scr_type_matter=1 and scr_type_in_NKCx^=scr_type then delete;
run;

proc sql;
create table hpv_st as select a.*, b.*  
from hpv_24_risk as a left join remiss_st as b on a.laboratory_id=b.laboratory_id and a.x_referral_type=b.r_referral_type and a.sample_type=b.r_sample_type;
quit;

data hpv_st;
set hpv_st;
if in_use_from^=. and sample_date<in_use_from then delete;
else if in_use_until^=. and sample_date>in_use_until then delete;
rename code=st_code;
if below_30_cyt=1 and hpv_age>=30 then delete;
if below_30_cyt=0 and hpv_age<30 then delete;
if scr_type_matter=1 and scr_type_in_NKCx^=scr_type then delete;
run;

proc sort data=hpv_nst;
by person_id sample_date;
proc sort data=hpv_st;
by person_id sample_date;
run;

data hpv_ref_all;
merge hpv_nst  hpv_st;
by person_id sample_date;
drop lab_id;
run;

proc sort data=hpv_ref_all;
by person_id sample_date;
run;

/*manage the code*/
data hpv_ref_all;
set hpv_ref_all;
if st_code^=. then code=st_code;
else code=nst_code;
run;
data hpv_ref_all;
set hpv_ref_all;
hpv_year=year(sample_date);
drop r_referral_type r_sample_type in_use_from in_use_until primary_analysis screeening category scr_type_in_NKCx  below_30_cyt scr_type_matter;
run;

/*merge the primary test method with the HPVref all data sets to have the screening hpv identified*/
proc sort data=hpv_ref_all;
by code;
run;
proc sort data=remiss1;
by code;
run;
data primary_method;
set remiss1;
keep code primary_analysis screeening category;
run;
data hpv_pm;
merge hpv_ref_all(in=a) primary_method;
by code;
if a ;
run;

data hpv_pm;
set hpv_pm;
if primary_analysis in ('HPV' 'Primary HPV') and screeening=1 then HPV_based=1;
else HPV_based=0;
run;


proc sort data=hpv_pm;
by person_id sample_date;
run;


/*find women with their first HPV-based screening test between age 50-70 years old*/
data hpv_based;
set hpv_pm;
if hpv_based=1 and HPVrisk^=.;
run;
/*NOTE: There were 4865763 observations read from the data set WORK.HPV_PM.
NOTE: The data set WORK.HPV_BASED has 3228868 observations and 35 variables.

*/
/*keep HPV test is after 2012 (>=2012)*/
data hpv_based;
set hpv_based;
hpv_year=year(sample_date);
where hpv_year>=2012;
drop test_year;
run;
/*NOTE: There were 3228726 observations read from the data set WORK.HPV_BASED.
      WHERE hpv_year>=2012;
NOTE: The data set WORK.HPV_BASED has 3228726 observations and 34 variables.

*/

proc sort data=hpv_based;
by person_id sample_date;
run;
proc sort data=hpv_based nodupkey out=first_hpv_based;
by person_id;
run;
/*NOTE: There were 3228726 observations read from the data set WORK.HPV_BASED.
NOTE: 884069 observations with duplicate key values were deleted.
NOTE: The data set WORK.FIRST_HPV_BASED has 2344657 observations and 34 variables.
*/



/*keep women who had their first HPV-based screening test at the age of 50-70*/
data first_hpv_5070;
set first_hpv_based;
if 50<=hpv_age<=70;
run;
/*NOTE: There were 2344657 observations read from the data set WORK.FIRST_HPV_BASED.
NOTE: The data set WORK.FIRST_HPV_5070 has 817349 observations and 34 variables.
*/




/*exclusion:
previous cancer
previous total hysterectomy
for women over 53 years old, any test for the past 4 years
for women under 53 years old (<=), any test for past 2.5 years
*/

/*build up a data set for find all the test
combine hpv ext_hpv and also cell_translated
need sample_date lopnr for hpv need risk system*/

data hpv;
set HPV_24;
hpv=1;
keep person_id hpvdiag sample_date reg_date hpv_type test_year hpv ;
if selftest=1 and reg_date^=. then sample_date=reg_date;
rename test_year=x_sample_yr;
run;
data ext_hpv;
set V_ncsr.nkc_ext_hpv_2024;
ext_hpv=1;
keep person_id hpvdiag x_sample_date x_reg_date x_sample_yr hpv_type ext_hpv;
where x_sample_yr<=2024;
run;
data ext_hpv;
set ext_hpv;
format sample_date yymmdd10. reg_date yymmdd10.;
sample_year=input(substr(x_sample_date,1,4),8.);
sample_month=input(substr(x_sample_date,6,2),8.);
sample_day=input(substr(x_sample_date,9,2),8.);
reg_year=input(substr(x_reg_date,1,4),8.);
reg_month=input(substr(x_reg_date,6,2),8.);
reg_day=input(substr(x_reg_date,9,2),8.);
sample_date=MDY(sample_month,sample_day,sample_year);
reg_date=MDY(reg_month,reg_day,reg_year);
drop sample_year sample_month sample_day reg_year reg_month reg_day;
drop x_sample_date x_reg_date;
run;

data all_hpv;/*used to find screening history and previous abnormality*/
set hpv ext_hpv;
run;
/*NOTE: There were 6179989 observations read from the data set WORK.HPV.
NOTE: There were 884564 observations read from the data set WORK.EXT_HPV.
NOTE: The data set WORK.ALL_HPV has 7064553 observations and 8 variables.

*/

data check;
set all_hpv;
test_year=year(sample_date);
reg_year=year(reg_date);
if abs(test_year-x_sample_yr)>1 then check=1;
if check=1 and reg_year=x_sample_yr then do; sample_date=reg_date; change=1;end;
run;

/*there is a systematic mistake of difference in sample_year and sample_date, 
cannot using this test, so here we delete them*/
data all_hpv;
set all_hpv;
if x_sample_yr=2016 and year(sample_date)=2002 then delete;
run;
/*
NOTE: There were 7064553 observations read from the data set WORK.ALL_HPV.
NOTE: The data set WORK.ALL_HPV has 7060385 observations and 8 variables.
*/


data all_hpv_risk;
set all_hpv;
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
if hpv_type='Not avai' then HPVmissing=1;
if hpv_type=:'Irrelevant' then HPVmissing=1;
if hpv_type=:'Se' then HPVmissing=1;
if hpv_type=:'human' then HPVmissing=1;
if hpv_type=:'brunn' then HPVmissing=1;
if HPV18_45=1 then do; HPV18=0; HPV45=0;end;
if HPV33=1 or HPV58=1 or HPV45=1 or HPV31=1 or HPV52=1 or HPV35=1 or HPV39=1 or HPV51=1 or HPV59=1 or HPV56=1 or HPV68=1 or HPV66=1 then
HPVHr=1;
if hpvdiag='Pos' then hpvdiag='POS';
run;
/*hpvtype missing but positive count as HrHPV positive*/

proc freq data=all_hpv_risk;
tables HPV16 HPV18 HPV45 HPV18_45 HPVHr HPVmissing/missing;
run;


data all_hpv_risk;
set all_hpv_risk;
if HPV16=1 and HPVDIAG='POS' then do; HPVrisk=4;Type16=1;end;
else if HPV18=1 and HPVDIAG='POS' then do; HPVrisk=3;Type18=1;end;
else if HPV18_45=1 and HPVDIAG='POS' then do;HPVrisk=3; Type18_45=1;end;
else if HPVHr=1 and HPVDIAG='POS' then do; HPVrisk=2; TypeHr=1;end;
else if HPVmissing=1 and HPVDIAG='POS' then HPVrisk=.;
else if HPVDIAG='POS' then HPVrisk=0;
else if HPVDIAG^='NEG' then HPVrisk=.;
else if HPVDIAG='NEG' then HPVrisk=0;
drop HPV16 HPV18  HPV33 HPV58 HPV31 HPV52 HPV35 HPV39 HPV51 HPV59 HPV56 HPV68 HPV66;
run;
proc freq data= all_hpv_risk;
table hpvrisk*HPVdiag/missing;
run;

/*delete all test without a valid result*/
data all_hpv_risk;
set all_hpv_risk;if HPVrisk=. then delete;
if sample_date=. then delete;
run;

/*
NOTE: There were 7060385 observations read from the data set WORK.ALL_HPV_RISK.
NOTE: The data set WORK.ALL_HPV_RISK has 7022239 observations and 17 variables.

*/
/*only keep one result per test per day per lopnr*/
proc sort data=all_hpv_risk;
by person_id sample_date descending HPVrisk;
run;
proc sort data=all_hpv_risk nodupkey;
by person_id sample_date;
run;
/*
NOTE: There were 7022239 observations read from the data set WORK.ALL_HPV_RISK.
NOTE: 1867220 observations with duplicate key values were deleted.
NOTE: The data set WORK.ALL_HPV_RISK has 5155019 observations and 17 variables.
*/

data hpv_his;
set all_hpv_risk;
keep person_id sample_date outcome hpv ext_hpv;
if HPVrisk=0 then outcome=0;else outcome=1;
run;

/*cytology until the end of 2024*/
data cyto;
set v_ncsr.nkc_trans_cell_2024;
keep person_id x_sample_date outcome snomed_severity;
if snomed_severity>=6 then outcome=1; else outcome=0;
if snomed_severity<=3 then delete;
where x_sample_yr<=2024;
run;
proc freq data=cyto;
table snomed_severity;
run;

data cyto;
set cyto;
format sample_date yymmdd10.;
sample_year=input(substr(x_sample_date,1,4),8.);
sample_month=input(substr(x_sample_date,6,2),8.);
sample_day=input(substr(x_sample_date,9,2),8.);
sample_date=MDY(sample_month,sample_day,sample_year);
drop sample_year sample_month sample_day ;
drop x_sample_date;
cyto=1;
run;

proc sort data=cyto;
by person_id sample_date descending outcome;run;
proc sort data=cyto nodupkey;
by person_id sample_date;
run;
/*NOTE: There were 25361658 observations read from the data set WORK.CYTO.
NOTE: 107685 observations with duplicate key values were deleted.
NOTE: The data set WORK.CYTO has 25253973 observations and 5 variables.
*/

data his_test;
set cyto hpv_his;
run;

proc sort data=his_test;
by person_id sample_date descending outcome cyto;
run;


/*exclude women with test within past 4 years*/
data pop;
set first_hpv_5070;
keep person_id sample_date hpv_age;
rename sample_date=first_hpv_date;
run;

data pop_exclude;
merge pop(in=a) his_test;
by person_id;
if hpv_age>53 and first_hpv_date-4*365.25<sample_date<first_hpv_date then exclude=1;
if hpv_age<=53 and first_hpv_date-2.5*365.25<sample_date<first_hpv_date then exclude=1;
if a;
run;
proc freq data=pop_exclude;
tables exclude;
run;
/*100713*/
proc sort data=pop_exclude;
by person_id descending exclude;
run;

proc sort data=pop_exclude nodupkey;
by person_id;
run;

data pop_exclude;
set pop_exclude;
keep person_id exclude;
run;

data population;
merge first_hpv_5070 pop_exclude;
by person_id;
run;

data population;
set population;
if exclude=1 then delete;
drop exclude;
run;
/*
NOTE: There were 817349 observations read from the data set WORK.POPULATION.
NOTE: The data set WORK.POPULATION has 742843 observations and 34 variables.

*/

/*cancer and total hysterectomy*/
data pop;
set population;
keep person_id sample_date birth_date;
rename sample_date=first_hpv_date;
run;

/*cancer*/
data cancer;
set v_ncsr.gyn_cancer_2024;
format diag_date yymmdd10.;
sample_year=input(substr(cxca_date,1,4),8.);
sample_month=input(substr(cxca_date,6,2),8.);
sample_day=input(substr(cxca_date,9,2),8.);
diag_date=MDY(sample_month,sample_day,sample_year);
drop sample_year sample_month sample_day ;
keep person_id diag_date figo_stage histo_type;
run;


data cancer;
set cancer;
cancer=1;
run;

/*keep first diagnosis*/
proc sort data=cancer nodupkey;
by person_id diag_date;
run;

proc sort data=cancer nodupkey;
by person_id ;
run;


/*using pad date adjust cancer date*/
data pad;
set v_ncsr.nkc_pad_translated_2024;
keep person_id PAD_sev x_sample_date;
where PAD_Sev>=5;
run;
data pad;
set pad;
format pad_date yymmdd10.;
sample_year=input(substr(x_sample_date,1,4),8.);
sample_month=input(substr(x_sample_date,6,2),8.);
sample_day=input(substr(x_sample_date,9,2),8.);
pad_date=MDY(sample_month,sample_day,sample_year);
drop sample_year sample_month sample_day ;
drop x_sample_date;
run;
proc sort data=pad;
by person_id pad_date;
run;
proc sort data=pad nodupkey;
by person_id ;
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
if diag_date<pad_date<diag_date+365.25 then diag_date=pad_date;
run;

data cancer;
set cancer_date;
drop pad_date pad_sev;
run;


/*the merge has no problem*/
data pre_cancer;
merge pop(in=a) cancer(in=b);
by person_id;
if diag_date<first_hpv_date;
if a and b;
keep person_id pre_cancer;
pre_cancer=1;
run;
/*NOTE: There were 742843 observations read from the data set WORK.POP.
NOTE: There were 6696 observations read from the data set WORK.CANCER.
NOTE: The data set WORK.PRE_CANCER has 117 observations and 2 variables.
*/



/*hysterectomy nkcx*/
data hyst_2024;
set v_ncsr.nkc_deregister_2024;
keep person_id x_dereg_from_date x_dereg_reason;
hyst=1;
if x_dereg_reason='Hysterektomi';
run;


data hyst;
set  hyst_2024;
format dereg_date yymmdd10.;
sample_year=input(substr(x_dereg_from_date,1,4),8.);
sample_month=input(substr(x_dereg_from_date,6,2),8.);
sample_day=input(substr(x_dereg_from_date,9,2),8.);
dereg_date=MDY(sample_month,sample_day,sample_year);
drop sample_year sample_month sample_day ;
drop x_dereg_from_date;
run;

proc sort data=hyst;
by person_id dereg_date;
run;
proc sort data=hyst nodupkey;
by person_id;
run;
proc sort data=pop;
by person_id;
run;

data hyst_pop;
merge hyst (in=a) pop (in=b);
by person_id;
if .<dereg_date<first_hpv_date;
if a and b;
run;

data pre_hyst_nkcx;
set hyst_pop;
keep person_id pre_hyst_nkc;
if x_dereg_reason='Hysterektomi';
pre_hyst_nkc=1;
run;


proc sort data=pre_hyst_nkcx nodupkey;
by person_id;
run;



data pre_hyst_all;
set pre_hyst_nkcx;
run;
/*
NOTE: There were 649 observations read from the data set WORK.PRE_HYST_NKCX.
NOTE: The data set WORK.PRE_HYST_ALL has 649 observations and 2 variables.
*/


/*exclusion data sets 
pre_cancer
pre_hyst_all

*/


/*previous abnormality*/
/*find women without any positive history for the past 10 years and*/

data test_10year;
merge pop(in=a) his_test;
by person_id;
if a;
if first_hpv_date-10*365.25<=sample_date<first_hpv_date;
run;


proc sort data=test_10year;
by person_id descending sample_date descending outcome;
run;

proc sort data=test_10year nodupkey;
by person_id descending sample_date;
run;
/*NOTE: There were 1230101 observations read from the data set WORK.TEST_10YEAR.
NOTE: 45832 observations with duplicate key values were deleted.
NOTE: The data set WORK.TEST_10YEAR has 1184269 observations and 9 variables.

*/

data test_10year;
set test_10year;
format next_test_date yymmdd10.;
by person_id;
next_test_date=lag(sample_date);
if first.lopnr then next_test_date=.;
run;

data test_10year;
set test_10year;
interval=(next_test_date-sample_date)/365.25;
run;



proc sql;
create table test_number  as
select  person_id,count(person_id)as test_count,sum(cyto) as num_cyto, sum(hpv) as num_hpv, sum(ext_hpv) as num_ext_hpv,  max(outcome) as positive, first_hpv_date, sum(interval) as total_interval,
birth_date from test_10year
group by person_id;
quit;

data test_number;
set test_number;
if positive=1 then delete;
run;
proc sort data=test_number nodupkey;
by person_id;
run;
data test_number;
set test_number;
if num_cyto<2  and total_interval<3 then delete;
run;
/*NOTE: There were 632981 observations read from the data set WORK.TEST_NUMBER.
NOTE: The data set WORK.TEST_NUMBER has 374386 observations and 9 variables.

*/

/*past 10 years have at least two negative test*/
data nopositive_10;
set test_number;
keep person_id nopositive_10;
nopositive_10=1;
run;

/*merge everything*/
proc sort data=pre_cancer;
by person_id;
proc sort data=pre_hyst_all;
by person_id;
run;
data population;
merge population (in=a) nopositive_10 pre_cancer  pre_hyst_all ;
by person_id;
if a;
run;

/*
NOTE: There were 742843 observations read from the data set WORK.POPULATION.
NOTE: There were 374386 observations read from the data set WORK.NOPOSITIVE_10.
NOTE: There were 117 observations read from the data set WORK.PRE_CANCER.
NOTE: There were 649 observations read from the data set WORK.PRE_HYST_ALL.
NOTE: The data set WORK.POPULATION has 742843 observations and 37 variables.


/*find previous positive HPV test (hpv his), keep the neareast test_date*/
Data previous_hpv;
merge hpv_his pop(in=a);
by person_id;
if sample_date<first_hpv_date;
if a;
run;

proc sort data=previous_hpv;
by person_id descending outcome descending sample_date;
run;
proc sort data=previous_hpv nodupkey out=previous_hpv_pos;
by person_id;
run;

data previous_hpv_pos;
set previous_hpv_pos;
if outcome=1;
previous_hpv_pos=1;
rename sample_date=previous_hpv_pos_date;
run;

data previous_hpv_pos;
set previous_hpv_pos;
keep person_id previous_hpv_pos previous_hpv_pos_date;
run;
/*
NOTE: There were 7410 observations read from the data set WORK.PREVIOUS_HPV_POS.
NOTE: The data set WORK.PREVIOUS_HPV_POS has 7410 observations and 3 variables.

*/

/*find previous positive cytology test (cyto), keep the neareast test_date*/

Data previous_cyto;
merge cyto pop(in=a);
by person_id;
if sample_date<first_hpv_date;
if a;
run;

proc sort data=previous_cyto;
by person_id descending outcome descending sample_date;
run;
proc sort data=previous_cyto nodupkey out=previous_cyto_pos;
by person_id;
run;

data previous_cyto_pos;
set previous_cyto_pos;
if outcome=1;
previous_cyto_pos=1;
rename sample_date=previous_cyto_pos_date;
run;
/*NOTE: There were 739864 observations read from the data set WORK.PREVIOUS_CYTO_POS.
NOTE: The data set WORK.PREVIOUS_CYTO_POS has 117573 observations and 8 variables.
*/

data previous_cyto_pos;
set previous_cyto_pos;
keep person_id previous_cyto_pos previous_cyto_pos_date;
run;
/*find previous LSIL*/
/*find previous CIN2*/
/*find previous CIN3; NKCx, cancer*/


/*find pad diagnosis before this hpv test*/ 
data pad_all;
set v_ncsr.nkc_pad_translated_2024;
keep pad_sev person_id pad_class x_sample_date TOPO3;
run;

data pad_all;
set pad_all;
format pad_date yymmdd10.;
sample_year=input(substr(x_sample_date,1,4),8.);
sample_month=input(substr(x_sample_date,6,2),8.);
sample_day=input(substr(x_sample_date,9,2),8.);
pad_date=MDY(sample_month,sample_day,sample_year);
drop sample_year sample_month sample_day ;
drop x_sample_date;
run;

proc sort data=pad_all;
by person_id descending TOPO3  pad_date;
run;


data LSIL HSIL CIN3 cancer;
set pad_all;
if pad_sev>=4 then output LSIL;
if pad_sev>=5 then output HSIL;
if pad_sev>=6 then output CIN3;
if pad_sev>=7 then output cancer;
run;

/*prioritize T83 first
when sorting by person_id descending TOPO3 descending pad_date*/


%macro padHis(type);
proc sort data=&type.;
by person_id;
run;
data pop_&type.;
merge pop (in=a) &type.;
by person_id;
if .<pad_date<first_hpv_date;
if a;
run;
proc sort data=pop_&type.;
by person_id descending TOPO3 descending pad_date;
proc sort data=pop_&type. nodupkey;
by person_id;
run;
data pop_&type.;
set pop_&type.;
rename pad_date=his_pad_&type._date pad_SEV=his_pad_&type._sev;
drop pad_class;
&type._diff=(first_hpv_date-pad_date)/365.25;
if .<&type._diff<=10 then his_&type._cat=1;
else if 10<&type._diff<=20 then his_&type._cat=2;
else if &type._diff>20 then his_&type._cat=3;
his_pad_&type.=1;
run;
%mend padhis;
%padhis(LSIL);
%padhis(HSIL);
%padhis(CIN3);
%padhis(cancer);


data pop_lsil;
set pop_lsil;
keep person_id his_pad_lsil_date his_pad_lsil;
run;
data pop_hsil;
set pop_hsil;
keep person_id his_pad_hsil_date his_pad_hsil;
run;
data pop_cin3;
set pop_cin3;
keep person_id his_pad_cin3_date his_pad_cin3;
run;
data pop_cancer;
set pop_cancer;
keep person_id his_pad_cancer_date his_pad_cancer;
run;

data previous_his_all;
merge pop (in=a) pop_LSIL pop_HSIL pop_CIN3 pop_cancer  previous_hpv_pos previous_cyto_pos ;
by person_id;
run;
/*
NOTE: There were 742843 observations read from the data set WORK.POP.
NOTE: There were 44676 observations read from the data set WORK.POP_LSIL.
NOTE: There were 30910 observations read from the data set WORK.POP_HSIL.
NOTE: There were 20401 observations read from the data set WORK.POP_CIN3.
NOTE: There were 2289 observations read from the data set WORK.POP_CANCER.
NOTE: There were 7410 observations read from the data set WORK.PREVIOUS_HPV_POS.
NOTE: There were 117573 observations read from the data set WORK.PREVIOUS_CYTO_POS.
NOTE: The data set WORK.PREVIOUS_HIS_ALL has 742843 observations and 15 variables.
*/

data previous_his_all;
set previous_his_all;
if his_pad_LSIL^=1 then his_pad_LSIL=0;
if his_pad_HSIL^=1 then his_pad_HSIL=0;
if his_pad_CIN3^=1 then his_pad_CIN3=0;
if his_pad_cancer^=1 then his_pad_cancer=0;
if previous_hpv_pos^=1 then previous_hpv_pos=0;
if previous_cyto_pos^=1 then previous_cyto_pos=0;
run;
/*find previous any abnormality*/

data previous_his_all;
set previous_his_all;
format previous_pos_date yymmdd10. ;
if previous_hpv_pos=1 or previous_cyto_pos=1 or his_pad_lsil=1   then do;
previous_pos=1;
previous_pos_date=max(of his_pad_LSIL_date his_pad_hSIL_date his_pad_cin3_date his_pad_cancer_date previous_hpv_pos_date previous_cyto_pos_date );
end;
diff_pos=(first_hpv_date-previous_pos_date)/365.25;
if .<diff_pos<=10 then previous_pos_cat=1;
else if 10<diff_pos<=20 then previous_pos_cat=2;
else if diff_pos>20 then previous_pos_cat=3;
drop diff_pos;

format pre_lsil_date yymmdd10.;
if (previous_hpv_pos=1 or previous_cyto_pos=1 or his_pad_lsil=1) and his_pad_hsil=0  then do;
pre_lsil=1;
pre_lsil_date=max(of his_pad_LSIL_date previous_hpv_pos_date previous_cyto_pos_date);
end;
diff_lsil=(first_hpv_date-pre_lsil_date)/365.25;
if .<diff_lsil<=10 then pre_lsil_cat=1;
else if 10<diff_lsil<=20 then pre_lsil_cat=2;
else if diff_lsil>20 then pre_lsil_cat=3;
drop diff_lsil;

format pre_hsil_date yymmdd10. ;
if his_pad_hsil=1 then do;
pre_hsil=1;
pre_hsil_date=max(of his_pad_hsil_date his_pad_cin3_date );
end;
diff_hsil=(first_hpv_date-pre_hsil_date)/365.25;
if .<diff_hsil<=10 then pre_hsil_cat=1;
else if 10<diff_hsil<=20 then pre_hsil_cat=2;
else if diff_hsil>20 then pre_hsil_cat=3;
drop diff_hsil;

format pre_cin3_date yymmdd10. ;
if his_pad_cin3=1 then do;
pre_cin3=1;
pre_cin3_date=max(of  his_pad_cin3_date );
end;
diff_cin3=(first_hpv_date-pre_cin3_date)/365.25;
if .<diff_cin3<=10 then pre_cin3_cat=1;
else if 10<diff_cin3<=20 then pre_cin3_cat=2;
else if diff_cin3>20 then pre_cin3_cat=3;
drop diff_cin3;
run;



/*here we already have the first round population and exposure 
population in population
exposure in previous_his_all*/


/*find second hpv test after the first one */
/*second hpv screening test, only do this among women with HPVrisk=0 
should be prmary HPV test (use hpv_based)*/
/*should be the first screening test after the first hpv test
defined it as 4 years apart and without any other test in between*/
/*match with pop*/
/*keep the first hpv test after first hpv*/
/*exclude women with other test in between*/

data secondhpv;
merge pop (in=a) hpv_based;
by person_id;
if sample_date>=first_hpv_date+4*365.25;
if a;
run;
/*
NOTE: Missing values were generated as a result of performing an operation on missing values.
      Each place is given by: (Number of times) at (Line):(Column).
      2352282 at 3025:31
NOTE: There were 742843 observations read from the data set WORK.POP.
NOTE: There were 3228726 observations read from the data set WORK.HPV_BASED.
NOTE: The data set WORK.SECONDHPV has 89634 observations and 35 variables.

*/

/*delete the invalid tests*/
data secondhpv;
set secondhpv;
if hpvrisk=. then delete;
run;
/*NOTE: There were 89634 observations read from the data set WORK.SECONDHPV.
NOTE: The data set WORK.SECONDHPV has 89634 observations and 35 variables.
*/

proc sort data=secondhpv;
by person_id sample_date;
run;
proc sort data=secondhpv nodupkey;
by person_id;
run;
/*NOTE: There were 89634 observations read from the data set WORK.SECONDHPV.
NOTE: 8563 observations with duplicate key values were deleted.
NOTE: The data set WORK.SECONDHPV has 81071 observations and 35 variables.
*/

data secondhpv;
set secondhpv;
rename sample_date=HPV_2nd_date;
run;

data secondround_exc;
merge secondhpv(in=a) his_test;
by person_id;
if first_hpv_date<sample_date<hpv_2nd_date;
if a;
run;
/*NOTE: There were 81071 observations read from the data set WORK.SECONDHPV.
NOTE: There were 30408992 observations read from the data set WORK.HIS_TEST.
NOTE: The data set WORK.SECONDROUND_EXC has 18318 observations and 41 variables.

*/

proc sort data=secondround_exc nodupkey;
by person_id;
run;
data secondround_exc;
set secondround_exc;
second_exclude=1;
keep person_id second_exclude;
run;
/*
NOTE: There were 7538 observations read from the data set WORK.SECONDROUND_EXC.
NOTE: The data set WORK.SECONDROUND_EXC has 7538 observations and 2 variables.

*/

data secondhpv;
merge secondhpv secondround_exc;
by person_id;
if second_exclude=1 then delete;
drop second_exclude;
run;
/*
NOTE: There were 81071 observations read from the data set WORK.SECONDHPV.
NOTE: There were 7538 observations read from the data set WORK.SECONDROUND_EXC.
NOTE: The data set WORK.SECONDHPV has 73533 observations and 35 variables.

*/

data secondhpv;
set secondhpv;
keep person_id hpv_2nd_date laboratory_id selftest  HPVrisk Type16 Type18 Type18_45 TypeHr multi_inf;
rename laboratory_id=lab_id_2nd selftest=selftest_2nd  HPVrisk=HPVrisk_2nd Type16=Type16_2nd Type18=Type18_2nd
Type18_45=Type18_45_2nd TypeHr=TypeHr_2nd multi_inf=multi_inf_2nd;
run;



/*here we have all the population informtion and exposure and exclusion information 
follow up use another code*/



data HPV_cohort;
merge population previous_his_all  secondhpv;
by person_id;
run;

data hpv_cohort;
set hpv_cohort;
drop st_code nst_code code test_year;
run;

data HPV_cohort;
set HPV_cohort;
format HPV_genotype $8. HPV_genotype_2nd $8.;
if Type16=1 then HPV_genotype='16';
else if Type18=1 then  HPV_genotype='18';
else if Type18_45=1 then  HPV_genotype='18/45';
else if TypeHr=1 then  HPV_genotype='Hr';
if Type16_2nd=1 then HPV_genotype_2nd='16';
else if Type18_2nd=1 then  HPV_genotype_2nd='18';
else if Type18_45_2nd=1 then  HPV_genotype_2nd='18/45';
else if TypeHr_2nd=1 then  HPV_genotype_2nd='Hr';
run;
data HPV_cohort;
set hpv_cohort;
if HPVrisk=0 then HPV=0;
else HPV=1;
if HPVrisk_2nd=0 then HPV_2nd=0;
else if HPVrisk_2nd^=. then HPV_2nd=1;
run;


/*check if more women with a previous abnormality had a total hysterectomy*/
proc freq data=hpv_cohort;
tables  pre_hyst_nkc*previous_pos/missing chisq;
where pre_cancer^=1;
run;
/*
women with previous_pos has higher chance with total hysterectomy chisq p<0.0001; 2%in women without, 3.8%in women with pre_abn, same with self reported one but that will be deleted anyway!
*/

data hpv_cohort;
set hpv_cohort;
where pre_cancer^=1 and pre_hyst_nkc^=1;
run;

proc freq data=hpv_cohort;
tables pre_cancer*pre_hyst_nkc/missing;
run;
/*
NOTE: There were 742085 observations read from the data set WORK.HPV_COHORT.
      WHERE (pre_cancer not = 1) and (pre_hyst_nkc not = 1);
NOTE: The data set WORK.HPV_COHORT has 742085 observations and 72 variables.
*/


data preabn.NKCx_nohys_cohort_260202;
set hpv_cohort;
run;

/*work on follow-up later next week*/

