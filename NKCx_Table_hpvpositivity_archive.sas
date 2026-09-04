/*File name..: Table 2 HPVpositivity*/
/*Study......: Phd project study 2 the HPV positivity and genotype after one HPV/two HPV test*/
/*Author.....: Qingyun Yao*/
/*Date.......: 2024/1031*/
/*Updated....: */
/*Purpose....: /*generate the HPV positivity in first round HPV-based screening*/
/*from 1995 or age 23

/*
libname NCSR  odbc complete="server=meb-sql02.meb.ki.se;driver=SQL Server Native Client 11.0;Trusted_Connection=Yes;database=NCSR" schema=NCSR ;
libname V_ncsr     odbc complete="server=meb-sql02.meb.ki.se;driver=SQL Server Native Client 11.0;Trusted_Connection=Yes;database=NCSR" schema=V_ncsr ;
libname PreAbn 'P:\ACCES\ACCES_Research\Qingyun\previous abnormality\Data';
libname Socmob4 odbc noprompt="dsn=kosmos;database=Cervix_socmob4" schema=clean_2023_2024;
 


				*/
/*Note.......: */
*------------------------------------------------------------------------;
/* Data used...: socmob4.nkc_hpv socmob4.nkc_trans_cell  self_samp_pop;

/* Data created.: pop_index_&date */
proc datasets library=work kill;quit;

libname NCSR odbc noprompt="dsn=kosmos; database=NCSR" schema=NCSR;
libname V_NCSR odbc noprompt="dsn=kosmos; database=NCSR" schema=V_NCSR ;
libname PreAbn 'P:\ACCES\ACCES_Research\Qingyun\previous abnormality\Data';




data analysis;
set preabn.NKCx_nohys_cohort_260202;
where noscr_his=0;
run;

proc freq data= preabn.NKCx_nohys_cohort_260202;
tables noscr_his;
run;

/*only used when do sensitivity analysis*/
data analysis;
set analysis;
if hpv_age<60;/*change this criteria based on which age group you are focusing on*/
run;

/*organize analysis table*/

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
run;



proc freq data=analysis;
table ref;
run;
/*%HPVp(ref);
%HPVp(previous_pos_cat);
%HPVp(previous_pos_not_LSIL_cat);
%HPVp(only_lsil_cat);
%HPVp(his_hsil_cat);
%HPVp(his_cin3_cat);
%HPVp(noscr_his);
rename the group*/




/*all population*/
proc summary data=analysis;
var HPVpositive;
output out=hpvp_all sum= mean=/autoname;
run;
data hpvp_all;
set hpvp_all;
char='All';
run;
%macro HPVp(group);
proc summary data=analysis;
var HPVpositive;
class &group.;
output out=hpvp_&group. sum= mean=/autoname;
run;
data hpvp_&group;
set hpvp_&group;
char=cat("&group.",&group.);
if _type_=0 then delete;
drop _type_ &group.;
run;
%mend;
%HPVp(ref);
%HPVp(previous_pos);
%HPVp(pre_ascus);
%HPVp(pre_lsil);
%HPVp(pre_hsil);
%HPVp(pre_cin3);
%HPVp(noscr_his);


%HPVp(previous_pos_cat);
%HPVp(pre_ascus_cat);
%HPVp(pre_lsil_cat);
%HPVp(pre_hsil_cat);
%HPVp(pre_cin3_cat);




/*nopositive_10*/
%macro hpvp10(group);
proc summary data=analysis;
var HPVpositive;
class &group. nopositive_10;
ways 2;
output out=hpvp_&group._10 sum= mean=/autoname;
run;


data hpvp_&group._10;
set hpvp_&group._10;
char=cat("&group.",&group.,'-',nopositive_10);
if _type_=0 then delete;
drop _type_ &group. nopositive_10;
run;

%mend;
%hpvp10(previous_pos_cat)
%hpvp10(pre_ascus_cat);
%hpvp10(pre_lsil_cat);
%hpvp10(pre_hsil_cat);
%hpvp10(pre_cin3_cat);



data hpvpositivity;
format char $30.;
set hpvp_:;
run;
data hpvpositivity;
set hpvpositivity;
if char='All' then type=1;
else if char='ref1' then type=2;
else type=3;
run;

proc sort data=hpvpositivity;
by type char;
run;
data hpvpositivity;
set hpvpositivity;
if index(char,'1-1') then delete;
if index(char,'0') then delete;
run;

data table_2_a;
set HPVpositivity;
keep char _freq_ HPVpositive;
rename _freq_=N;
HPVpositive= cat(compress(put (hpvpositive_sum,8.)),' (',compress(put(HPVpositive_mean*100,8.1)),'%)');
run;


/*run the logistic regression*/
/*adjusted for sample year, region, birth-cohort*/

%macro HPVrr(exposure, outcome, num);
data regression;
set analysis;
keep person_id &exposure._cat &exposure. ref region hpv_year birthco &outcome. nopositive_10 selftest;
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
class &exposure.(ref='0') birthco  hpv_year selftest;
model &outcome.(event='1')= &exposure. birthco  hpv_year /dist=binomial link=log;
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
class &exposure._cat(ref='0') birthco  hpv_year selftest;
model &outcome.(event='1')= &exposure._cat birthco  hpv_year /dist=binomial link=log;
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
class &exposure._cat(ref='0') birthco  hpv_year selftest;
model &outcome.(event='1')= &exposure._cat birthco  hpv_year /dist=binomial link=log;
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


%HPVrr(previous_pos, HPVpositive, 1);
%HPVrr(pre_ascus, HPVpositive, 2);
%HPVrr(pre_lsil, HPVpositive, 3);
%HPVrr(pre_hsil, HPVpositive, 4);
%HPVrr(pre_cin3, HPVpositive, 5);





data regression;
set analysis;
keep person_id noscr_his ref region hpv_year birthco HPVpositive nopositive_10;
if noscr_his=0 and ref=1 then noscr_his=0;
else if noscr_his=0 then noscr_his=.;
run;

proc genmod data=regression;
class noscr_his(ref='0') ;
model HPVpositive(event='1')= noscr_his/dist=binomial link=log;
ods output ParameterEstimates=estimate_org_3;
run;

proc genmod data=regression;
class noscr_his(ref='0') birthco region hpv_year ;
model HPVpositive(event='1')= noscr_his birthco  hpv_year /dist=binomial link=log;
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
create table table_2 as select a.*, b.* from  table_2_a as a left join outcome_all as b
on  upcase(a.char)=upcase(b.char);
quit;


/*remember to change the output data file name*/
%let mydir=P:\ACCES\ACCES_Research\Qingyun\previous abnormality\Doc\NKCx_260604_exnsh\;
title 'Table2 first round HPVpositivity ';
ods rtf file="&mydir.pre_abn_chr_table2_newpop_&sysdate..rtf";
proc print data=table_2 noobs;run;
title '';
ods rtf close;
/*export results*/
