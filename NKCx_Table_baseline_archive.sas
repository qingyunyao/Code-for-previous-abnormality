/*File name..: Table 1 baseline characteristics*/
/*Study......: Phd project study 2 the HPV positivity and genotype after one HPV/two HPV test*/
/*Author.....: Qingyun Yao*/
/*Date.......: 2024/1031*/
/*Updated....: */
/*Purpose....: /*generate the table 1 for the second study*/
/*from 1995 or age 23

/*
libname NCSR  odbc complete="server=meb-sql02.meb.ki.se;driver=SQL Server Native Client 11.0;Trusted_Connection=Yes;database=NCSR" schema=NCSR ;
libname V_ncsr     odbc complete="server=meb-sql02.meb.ki.se;driver=SQL Server Native Client 11.0;Trusted_Connection=Yes;database=NCSR" schema=V_ncsr ;
libname PreAbn 'P:\ACCES\ACCES_Research\Qingyun\previous abnormality\Data';
libname Socmob4 odbc noprompt="dsn=kosmos;database=Cervix_socmob4" schema=clean_2023_2024;
libname S4 odbc noprompt="dsn=kosmos;database=Cervix_socmob4" schema=clean_2022;


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




/*number of women in each group*/
/*using the cancer registry data*/
/*the end of HPV test including in */

/*check birth country information*/

data analysis;
set preabn.NKCx_nohys_cohort_260202;
where noscr_his=0;
run;
/*add birth country and education information*/
data pop;
set analysis;
keep person_id;
run;
proc sql;
create table home_country as select a.*, b.birth_other_country from pop as a inner join v_ncsr.nkc_pop_2024 as b on a.person_id=b.person_id;
quit;



data home_country;
set home_country;
if birth_other_country in ('FINLAND' 'DANMARK' 'NORGE' 'ISLAND') then birth_country='NORDIC';
else if birth_other_country not in ('SVERIGE' ' ') then birth_country='OTHER';
else birth_country="SWEDEN";
run;
proc freq data=home_country;
table birth_other_country birth_country;
run;


data home_country;
set home_country;
keep person_id birth_country;
run;


proc sort data=home_country;
by person_id;
run;


data analysis;
merge analysis home_country ;
by person_id;
run;
proc freq data=analysis;
table birth_country;
run;

data analysis;
set analysis;
if birth_country='' then birth_country='OTHER';
run;

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
else ref=0;
run;

proc freq data=analysis;
table ref*noscr_his*previous_pos/missing;
run;



/*women with previous abnormality/ women without previous abnormality/women with previous abn but not LSIL*/
/*women with only LSIL*/ /*women with cin2+*/ /*women with CIN3+*/
proc summary data=analysis;
var ref previous_pos pre_ascus pre_lsil pre_hsil pre_cin3 noscr_his ;
class age_group  birth_country region selftest sample_period secondround birthco; ways 1;
output out=baseline  mean=/autoname;
run;
proc summary data=analysis;
var ref previous_pos pre_ascus pre_lsil pre_hsil pre_cin3 noscr_his ;
output out=all  mean=/autoname;
run;
proc summary data=analysis;
var hpv_age test_age_2nd;
class ref previous_pos pre_ascus pre_lsil pre_hsil pre_cin3 noscr_his ;
ways 1;
output out=age  mean= std=/autoname;
run;
proc summary data=analysis;
var hpv_age test_age_2nd;
output out=age_all  mean= std=/autoname;
run;
data age;
set age_all age;
run;
proc sort data=age;
by _type_ descending previous_pos descending pre_ascus descending pre_lsil descending pre_hsil descending pre_cin3
descending noscr_his descending ref;
run;
proc sort data=age nodupkey;
by _type_;
run; 

data age_table;
set age;
format var $24.;
keep var age1 age2;
if _type_=0 then var='all';
if _type_=1 then var='noscr_his';
if _type_=2 then var='pre_cin3';
if _type_=4 then var='pre_hsil';
if _type_=8 then var='pre_lsil';
if _type_=16 then var='pre_ascus';
if _type_=32 then var='previous_pos';
if _type_=64 then var='ref';
age1=cat(compress(put(hpv_age_mean,8.1)),' (',compress(put(hpv_age_stddev,8.1)),')');
age2=cat(compress(put(test_age_2nd_mean,8.1)),' (',compress(put(test_age_2nd_stddev,8.1)),')');
run;

proc transpose data=age_table out=age1;
var age1;
id var;
run;
proc transpose data=age_table out=age2;
var age2;
id var;
run;

data age_trans;
set age1 age2;
rename _NAME_=character;
run;


data baseline;
set baseline;
number_preivous_pos=_freq_*previous_pos_mean;
number_pre_ascus=_freq_*pre_ascus_mean;
number_pre_LSIL=_freq_*pre_lsil_mean;
number_pre_hsil=_freq_*pre_hsil_mean;
number_pre_cin3=_freq_*pre_cin3_mean;
number_noscr_his=_freq_*noscr_his_mean;
number_ref=_freq_*ref_mean;
run;
data all;
set all;
number_preivous_pos=_freq_*previous_pos_mean;
number_pre_ascus=_freq_*pre_ascus_mean;
number_pre_LSIL=_freq_*pre_lsil_mean;
number_pre_hsil=_freq_*pre_hsil_mean;
number_pre_cin3=_freq_*pre_cin3_mean;
number_noscr_his=_freq_*noscr_his_mean;
number_ref=_freq_*ref_mean;
run;

data baseline;
set all baseline;
run;

data baseline;
set baseline;
format character $20.;
if _type_=0 then character='All';
if age_group=1 then character='age_50-59';
if age_group=2 then character='age_60-70';
if birth_country^='' then character=cat('birth-',substr(birth_country,1,8));
if region^='' then character=region;
if selftest=0 then character='clinical sampling';
if selftest=1 then character='self sampling';
if sample_period=1 then character='2012-2016';
if sample_period=2 then character='2017-2019';
if sample_period=3 then character='2020-2024';
if secondround=0 then character='----';
if secondround=1 then character='secondround';
if birthco=1 then character='birthco1';
if birthco=2 then character='birthco2';
if birthco=3 then character='birthco3';
if birthco=4 then character='birthco4';
if birthco=5 then character='birthco5';
if birthco=6 then character='birthco6';
run;

data table1;
set baseline;
keep character All previous_pos pre_ascus pre_lsil pre_hsil pre_cin3 noscr_his ref;
all=cat(compress(put(_freq_,8.)), ' (', compress(put(_freq_/742085*100,8.1)),'%)');
previous_pos=cat(compress(put(number_preivous_pos,8.)),' (',compress(put(previous_pos_mean*100,8.1)),'%)');
pre_ascus=cat(compress(put(number_pre_ascus,8.)),' (',compress(put(pre_ascus_mean*100,8.1)),'%)');
pre_lsil=cat(compress(put(number_pre_lsil,8.)),' (',compress(put(pre_lsil_mean*100,8.1)),'%)');
pre_hsil=cat(compress(put(number_pre_hsil,8.)),' (',compress(put(pre_hsil_mean*100,8.1)),'%)');
pre_cin3=cat(compress(put(number_pre_cin3,8.)),' (',compress(put(pre_cin3_mean*100,8.1)),'%)');
noscr_his=cat(compress(put(number_noscr_his,8.)),' (',compress(put(noscr_his_mean*100,8.1)),'%)');
ref=cat(compress(put(number_ref,8.)),' (',compress(put(ref_mean*100,8.1)),'%)');
run;
/*this table is not col% it presented row%, in the final table we present col%, was calculated using excel*/

data table1 ;
set table1 age_trans;
run;

/*export results*/
title 'Table1 baseline characteristics';
ods rtf file="&mydir.pre_ abn_chr_table1_&sysdate..rtf";
proc print data=table1 noobs;run;
title '';
ods rtf close;
