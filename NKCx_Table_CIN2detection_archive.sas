/*File name..: Table CIN2+ detection*/
/*Study......: Phd project study 2 the HPV positivity following the first HPV-based screening test*/
/*Author.....: Qingyun Yao*/
/*Date.......: 2024/1104*/
/*Updated....: */
/*Purpose....: /*generate the table table for CIN2+ detection within 12 months*/
/*
/*20241105 detection rate of CIN2+ within 12 months*/

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
%let mydir=P:\ACCES\ACCES_Research\Qingyun\previous abnormality\Doc\NKCx_260604_exnsh\;





data analysis;
set preabn.NKCx_nohys_cohort_260202;
where noscr_his=0;
run;

data pad_outcome;
set preabn.NKCx_nohys_longterm_fu;
keep person_id fu_pad_12 fu_cin2_12;
if fu_pad_12=. then fu_pad_12=0;
where noscr_his=0;
run;
data analysis;
merge analysis pad_outcome;
by person_id;
run;
/*sensitivity analysis*/
data analysis;
set analysis;
if hpv_age<60;
run;
/*skip above data code for main analysis*/

data analysis;
set analysis;
if HPVrisk>1 then HPVpositive=1;
else HPVpositive=0;
if HPVrisk_2nd>1 then HPVpositive_2nd=1;
else HPVpositive_2nd=0;
if laboratory_id in ('088' '999') then region='Stockholm-Gotland';
else if laboratory_id in ('211' '231' '237' '251') then region='Southeast Sweden';
else if laboratory_id in ('241' '271' '411' '417') then region='South Sweden';
else if laboratory_id in ('621' '631' '641' '651') then region='North Sweden';
else if laboratory_id in ('421' '427' '501' '507' '511' '517' '521' '527' '531' '537') then region='West Sweden';
else if laboratory_id in ('121' '127' '131' '541' '551' '561' '567' '571' '577' '611') then region='Middle Sweden';
test_age_2nd=year(HPV_2nd_date)-year(birth_date);
if 50<=hpv_age<=59 then age_group=1;
else if 60<=hpv_age<=70 then age_group=2;
/*50-54, 55-59, 60-64, 65-70*/
if 50<=hpv_age<=54 then age_group_4=1;
else if 55<=hpv_age<=59 then age_group_4=2;
else if 60<=hpv_age<=64 then age_group_4=3;
else if 65<=hpv_age<=70 then age_group_4=4;
if hpv_year<=2016 then sample_period=1;
else if hpv_year<=2019 then sample_period=2;
else sample_period=3;
if hpv_year<=2019 then calender_period=1;
else calender_period=2;
run;

data analysis;
set analysis;
if HPVrisk=0 and HPVrisk_2nd^=. then secondround=1;
else secondround=0;
run;

data analysis;
set analysis;
if selftest=. then selftest=0;
if pre_ascus=. then pre_ascus=0;
if pre_lsil=. then pre_lsil=0;
if pre_hsil=. then pre_hsil=0;
if pre_cin3=. then pre_cin3=0;
if noscr_his=0 and previous_pos=0 then ref=1;
format self_hpv best12.;
self_hpv=selftest;
drop selftest;
run;
proc freq data=analysis;
table ref;
run;




/*for the outcome of HSIL detection rate within 12 months*/
data analysis;
set analysis;
if first_hpv_date<='01Jan2023'd;
rename self_hpv=selftest;
run;
/*
NOTE: There were 722195 observations read from the data set WORK.ANALYSIS.
NOTE: The data set WORK.ANALYSIS has 592448 observations and 93 variables.
*/

/*need to add the number of women with histopathology test*/
/*all population*/
proc summary data=analysis;
var selftest fu_pad_12 fu_cin2_12 ;
output out=cin2_all sum= mean= /autoname;
run;

data cin2_all;
set cin2_all;
format char $24.;
char='All';
run;
%macro cin2(group);
proc summary data=analysis;
var selftest fu_pad_12 fu_cin2_12 ;
class &group.;
output out=cin2_&group. sum= mean=/autoname;
run;
data cin2_&group;
set cin2_&group;
char=cat("&group.",&group.);
if _type_=0 then delete;
drop _type_ &group.;
run;
%mend;
%cin2(ref);
%cin2(previous_pos);
%cin2(pre_ascus);
%cin2(pre_lsil);
%cin2(pre_hsil);
%cin2(pre_cin3);
%cin2(noscr_his);

%cin2(previous_pos_cat);
%cin2(pre_ascus_cat);
%cin2(pre_lsil_cat);
%cin2(pre_hsil_cat);
%cin2(pre_cin3_cat);


/*nopositive_10*/
%macro cin2_10(group);
proc summary data=analysis;
var self_HPV hispa_12 pad_cin2_12 ;
class &group. nopositive_10;
ways 2;

output out=cin2_&group._10 sum= mean=/autoname;
run;

data cin2_&group._10;
set cin2_&group._10;
char=cat("&group.",&group.,'-',nopositive_10);
if _type_=0 then delete;
drop _type_ &group. nopositive_10;
run;


%mend;
%cin2_10(previous_pos)
%cin2_10(pre_ascus);
%cin2_10(pre_lsil);
%cin2_10(pre_hsil);
%cin2_10(pre_cin3);
data pad_cin2_12;
format char $30.;
set cin2_:;
run;
data pad_cin2_12;
set pad_cin2_12;
if char='All' then type=1;
else if char='ref1' then type=2;
else type=3;
run;

proc sort data=pad_cin2_12;
by type char;
run;
data pad_cin2_12;
set pad_cin2_12;
if index(char,'1-1') then delete;
if index(char,'0') then delete;
run;

data table_4_a;
set pad_cin2_12;
keep char _freq_ self_hpv hispa_12 pad_cin2_12 ;
*det_test=pad_cin2_6_sum/hispa_6_sum*100;
rename _freq_=N;
self_hpv= cat(compress(put (selftest_sum,8.)),' (',compress(put(selftest_mean*100,8.1)),'%)');
hispa_12= cat(compress(put (fu_pad_12_sum,8.)),' (',compress(put(fu_pad_12_mean*100,8.1)),'%)');
pad_cin2_12= cat(compress(put (fu_cin2_12_sum,8.)),' (',compress(put(fu_cin2_12_mean*100,8.2)),'%)');
run;

data analysis;
set analysis;
rename hpv_year= sample_year;
run;



/*run the logistic regression*/
/*adjusted for sample year, region, birth-cohort*/

%macro cin2rr(exposure, outcome, num);
data regression;
set analysis;
keep person_id &exposure._cat &exposure. ref region sample_year birthco &outcome. nopositive_10 selftest;
if &exposure._cat=. and ref=1 then &exposure._cat=0;
if &exposure.=0 and ref=1 then &exposure.=0;
else if &exposure.=0 then &exposure.=.;
run;

proc genmod data=regression;
class &exposure.(ref='0') ;
model &outcome.(event='1')= &exposure./dist=binomial link=log;
ods output ParameterEstimates=estimate_org_3;
run;

proc genmod data=regression;
class &exposure.(ref='0') birthco  sample_year ;
model &outcome.(event='1')= &exposure. birthco  sample_year /dist=binomial link=log;
ods output ParameterEstimates=estimate_adj_3;
run;

data estimate_adj_3;
set estimate_adj_3;
if parameter="&exposure.";
run;

proc genmod data=regression;
class &exposure._cat(ref='0') ;
model &outcome.(event='1')= &exposure._cat/dist=binomial link=log;
ods output ParameterEstimates=estimate_org_1;
run;
proc genmod data=regression;
class &exposure._cat(ref='0') birthco  sample_year ;
model &outcome.(event='1')= &exposure._cat birthco  sample_year /dist=binomial link=log;
ods output ParameterEstimates=estimate_adj_1;
run;
data estimate_adj_1;
set estimate_adj_1;
if parameter="&exposure._cat";
run;
proc genmod data=regression;
class &exposure._cat(ref='0') ;
model &outcome.(event='1')= &exposure._cat/dist=binomial link=log;
ods output ParameterEstimates=estimate_org_2;
where nopositive_10=1 ;
run;
data estimate_org_2;
set estimate_org_2;
if parameter="&exposure._cat";
nopositive_10=1;
if level1=1 then delete;
run;
proc genmod data=regression;
class &exposure._cat(ref='0') birthco  sample_year ;
model &outcome.(event='1')= &exposure._cat birthco  sample_year /dist=binomial link=log;
ods output ParameterEstimates=estimate_adj_2;
where nopositive_10=1 ;
run;
data estimate_adj_2;
set estimate_adj_2;
nopositive_10=1;
if parameter="&exposure._cat";
if level1=1 then delete;
run;

data outcome_org;
set estimate_org:;
format char $35. crude_rr $35.;
char=cat(compress(parameter),compress(level1));
if nopositive_10=1 then char=cat(compress(parameter),compress(level1),'-',compress(nopositive_10));
crude_rr=cat(compress(put(exp(estimate),8.2)),' (',compress(put(exp(lowerwaldcl),8.2)),', ',compress(put(exp(upperwaldcl),8.2)),')');
keep char crude_rr;
run;
data outcome_adj;
set estimate_adj:;
format char $35. adj_rr $35.;
char=cat(compress(parameter),compress(level1));
if nopositive_10=1 then char=cat(compress(parameter),compress(level1),'-',compress(nopositive_10));
adj_rr=cat(compress(put(exp(estimate),8.2)),' (',compress(put(exp(lowerwaldcl),8.2)),', ',compress(put(exp(upperwaldcl),8.2)),')');
rename probchisq=adj_p;
keep char adj_rr probchisq;
run;
proc sort data=outcome_org;
by char;
proc sort data=outcome_adj;
by char;
run;
data outcome_&num.;
merge outcome_org outcome_adj;
by char;
if index(char,'0') then delete;
if adj_rr='' then delete;
run;
%mend;


%cin2rr(previous_pos, fu_cin2_12, 1);
%cin2rr(pre_ascus, fu_cin2_12, 2);
%cin2rr(pre_lsil, fu_cin2_12, 3);
%cin2rr(pre_hsil, fu_cin2_12, 4);
%cin2rr(pre_cin3, fu_cin2_12, 5);



data regression;
set analysis;
keep person_id noscr_his ref region sample_year birthco fu_cin2_12 nopositive_10 selftest;
if noscr_his=0 and ref=1 then noscr_his=0;
else if noscr_his=0 then noscr_his=.;
run;

proc genmod data=regression;
class noscr_his(ref='0') ;
model fu_cin2_12(event='1')= noscr_his/dist=binomial link=log;
ods output ParameterEstimates=estimate_org_3;
run;

proc genmod data=regression;
class noscr_his(ref='0') birthco  sample_year ;
model fu_cin2_12(event='1')= noscr_his birthco  sample_year /dist=binomial link=log;
ods output ParameterEstimates=estimate_adj_3;
run;

data estimate_adj_3;
set estimate_adj_3;
if parameter='noscr_his';
run;

data outcome_org;
set estimate_org_3;
char=cat(compress(parameter),compress(level1));
crude_rr=cat(compress(put(exp(estimate),8.2)),' (',compress(put(exp(lowerwaldcl),8.2)),', ',compress(put(exp(upperwaldcl),8.2)),')');
keep char crude_rr;
run;

data outcome_adj;
set estimate_adj_3:;
format char $35. adj_rr $35.;
char=cat(compress(parameter),compress(level1));
adj_rr=cat(compress(put(exp(estimate),8.2)),' (',compress(put(exp(lowerwaldcl),8.2)),', ',compress(put(exp(upperwaldcl),8.2)),')');
rename probchisq=adj_p;
keep char adj_rr probchisq;
run;
proc sort data=outcome_org;
by char;
proc sort data=outcome_adj;
by char;
run;
data outcome_6;
merge outcome_org outcome_adj;
by char;
if index(char,'0') then delete;
if adj_rr='' then delete;
run;


data outcome_all;
set outcome_1 - outcome_6;
run;


proc sql;
create table table_4 as select a.*, b.* from  table_4_a as a left join outcome_all as b
on  upcase(a.char)=upcase(b.char);
quit;



title 'table4 HSIL detection within 12 month';
ods rtf file="&mydir.pre_abn_Table4_CIN2 detection rate&sysdate..rtf";
proc print data=table_4 noobs;run;
title '';
ods rtf close;
/*export results*/
