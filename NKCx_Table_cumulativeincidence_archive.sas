/*File name..: identify the longterm follow-up of the CIN2+ and ICC*/
/*Study......: Phd project study 2 the HPV positivity following first HPV-based screening test*/
/*Author.....: Qingyun Yao*/
/*Date.......: 2024/1121*/
/*Updated....: */
/*Purpose....: /*identify the end point of CIN2+ and ICC after first HPV */
/*The end of follow-up of CIN2+ is the last registered test date in NKCx or the detection of CIN2+
Then end of follow-up of ICC is the end of the cancer registry

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



data  pop_sur_time ;
set preabn.NKCx_nohys_longterm_fu;
where noscr_his=0;
run;


/*limit population when do the analysis*/

/*sensitivity analysis*/


/*data presentation of person_time*/
data pop_sur_time;
set pop_sur_time;
if previous_pos=. then previous_pos=0;
if previous_pos_cat=. then previous_pos_cat=0;
if previous_pos=0 and noscr_his=0 then ref=1;
else ref=0;

if sur_cin2=0 then sur_cin2=.;
if hpvrisk=0 then hpv=0;
else hpv=1;
run;
data pop_sur_time;
set pop_sur_time;
if Figo in ('IA1' 'IA2') then earlystage=1;
else if figo in ('IB1' 'IB2' 'IB3') then earlystage=1;
else if figo^='' then earlystage=3;
run;
data pop_sur_time;
set pop_sur_time;
if Figo in ('IB1' 'IB2' 'IB3') then IB=1;
else if figo^='' then IB=0;
if hpvrisk=0 then hpv=0;
else hpv=1;
run;
proc freq data=pop_sur_time;
tables IB IB*previous_pos IB*previous_pos_cat IB*noscr_his IB*ref/missing;
where cancer=1;
run;
/*all women*/
proc summary data=pop_sur_time;
var sur_cin2 sur_cancer fu_cin2 cancer;
class HPV;
output out=person_time_whopop sum= median(sur_cin2 sur_cancer)= q1(sur_cin2 sur_cancer)= q3(sur_cin2 sur_cancer)= /autoname;
run;
data person_time_whopop;
set person_time_whopop;
CIN2_IR=fu_cin2_sum/sur_cin2_sum*100000;
cancer_IR=cancer_sum/sur_cancer_sum*100000;
run;
/*negative women*/
proc summary data=pop_sur_time;
var sur_cin2 sur_cancer fu_cin2 cancer;
class previous_pos noscr_his ref;
ways 1;
output out=person_time_negative sum= median(sur_cin2 sur_cancer)= q1(sur_cin2 sur_cancer)= q3(sur_cin2 sur_cancer)= /autoname;
where HPVrisk=0;
run;
proc sort data=person_time_negative;
by descending ref descending noscr_his descending previous_pos;
run;
proc sort data=person_time_negative nodupkey;
by _TYPE_;
run;

data person_time_negative;
set person_time_negative;
CIN2_IR=fu_cin2_sum/sur_cin2_sum*100000;
cancer_IR=cancer_sum/sur_cancer_sum*100000;
run;

/*positive*/

proc summary data=pop_sur_time;
var sur_cin2 sur_cancer fu_cin2 cancer;
class previous_pos noscr_his ref;
ways 1;
output out=person_time_positive sum= median(sur_cin2 sur_cancer)= q1(sur_cin2 sur_cancer)= q3(sur_cin2 sur_cancer)= /autoname;
where HPVrisk>0;
run;
proc sort data=person_time_positive;
by descending ref descending noscr_his descending previous_pos;
run;
proc sort data=person_time_positive nodupkey;
by _TYPE_;
run;

data person_time_positive;
set person_time_positive;
CIN2_IR=fu_cin2_sum/sur_cin2_sum*100000;
cancer_IR=cancer_sum/sur_cancer_sum*100000;
run;

/*all*/

proc summary data=pop_sur_time;
var sur_cin2 sur_cancer fu_cin2 cancer;
class previous_pos noscr_his ref;
ways 1;
output out=person_time_all sum= median(sur_cin2 sur_cancer)= q1(sur_cin2 sur_cancer)= q3(sur_cin2 sur_cancer)= max(sur_cin2 sur_cancer)= min(sur_cin2 sur_cancer)=/autoname;
run;
proc sort data=person_time_all;
by descending ref descending noscr_his descending previous_pos;
run;
proc sort data=person_time_all nodupkey;
by _TYPE_;
run;

data person_time_all;
set person_time_all;
CIN2_IR=fu_cin2_sum/sur_cin2_sum*100000;
cancer_IR=cancer_sum/sur_cancer_sum*100000;
run;



/*summary based on the time elapse*/

proc summary data=pop_sur_time;
var sur_cin2 sur_cancer fu_cin2 cancer;
class previous_pos_cat;
ways 1;
output out=person_time_cat_all sum= median(sur_cin2 sur_cancer)= q1(sur_cin2 sur_cancer)= q3(sur_cin2 sur_cancer)= /autoname;
*where noscr_his=0;
run;

data person_time_cat_all;
set person_time_cat_all;
CIN2_IR=fu_cin2_sum/sur_cin2_sum*100000;
cancer_IR=cancer_sum/sur_cancer_sum*100000;
run;

/*among women tested negative*/
proc summary data=pop_sur_time;
var sur_cin2 sur_cancer fu_cin2 cancer;
class previous_pos_cat;
ways 1;
output out=person_time_cat_negative sum= median(sur_cin2 sur_cancer)= q1(sur_cin2 sur_cancer)= q3(sur_cin2 sur_cancer)= /autoname;
where  HPVrisk=0;
run;

data person_time_cat_negative;
set person_time_cat_negative;
CIN2_IR=fu_cin2_sum/sur_cin2_sum*100000;
cancer_IR=cancer_sum/sur_cancer_sum*100000;
run;
/*among women tested positive*/

proc summary data=pop_sur_time;
var sur_cin2 sur_cancer fu_cin2 cancer;
class previous_pos_cat;
ways 1;
output out=person_time_cat_positive sum= median(sur_cin2 sur_cancer)= q1(sur_cin2 sur_cancer)= q3(sur_cin2 sur_cancer)= /autoname;
where  HPVrisk>0;
run;

data person_time_cat_positive;
set person_time_cat_positive;
CIN2_IR=fu_cin2_sum/sur_cin2_sum*100000;
cancer_IR=cancer_sum/sur_cancer_sum*100000;
run;


/*merge all of those/name the group*/
data person_time_all;
set person_time_all;
format char $30.;
if _TYPE_=1 then char='ref';
if _TYPE_=2 then char='no screening history';
if _TYPE_=4 then char='previous_pos';
drop previous_pos noscr_his ref;
run;
data person_time_negative;
set person_time_negative;
format char $30. HPVoutcome $30.;
if _TYPE_=1 then char='ref';
if _TYPE_=2 then char='no screening history';
if _TYPE_=4 then char='previous_pos';
drop previous_pos noscr_his ref;
HPVoutcome='Negative';
run;
data person_time_positive;
set person_time_positive;
format char $30. HPVoutcome $30.;
if _TYPE_=1 then char='ref';
if _TYPE_=2 then char='no screening history';
if _TYPE_=4 then char='previous_pos';
drop previous_pos noscr_his ref;
HPVoutcome='Positive';
run;

data person_time_cat_all;
set person_time_cat_all;
format char $30.;
char=cat('previous_pos','-', previous_pos_cat);
run;

data person_time_cat_negative;
set person_time_cat_negative;
format char $30. HPVoutcome $30.;
char=cat('previous_pos','-', previous_pos_cat);
HPVoutcome='Negative';
run;
data person_time_cat_positive;
set person_time_cat_positive;
format char $30. HPVoutcome $30.;
char=cat('previous_pos','-', previous_pos_cat);
HPVoutcome='Positive';
run;
data person_time_whopop;
set person_time_whopop;
format char $30. HPVoutcome $30.;
char='All';
if hpv=1 then HPVoutcome='Positive';
if hpv=0 then HPVoutcome='Negative';
run;

/*calculate 95%CI to ref group*/
data person_time_all;
set person_time_all;
ref_cin2_ir=lag(cin2_ir); ref_cin2_sum=lag(fu_cin2_sum); ref_cancer_ir=lag(cancer_ir); ref_cancer_sum=lag(cancer_sum);
if char='ref' then do; ref_cin2_ir=cin2_ir; ref_cin2_sum=fu_cin2_sum; ref_cancer_ir=cancer_ir; ref_cancer_sum=cancer_sum;end;
ref_cin2_ir=lag(ref_cin2_ir); ref_cin2_sum=lag(ref_cin2_sum); ref_cancer_ir=lag(ref_cancer_ir); ref_cancer_sum=lag(ref_cancer_sum);
if char='ref' then do; ref_cin2_ir=cin2_ir; ref_cin2_sum=fu_cin2_sum; ref_cancer_ir=cancer_ir; ref_cancer_sum=cancer_sum;end;
run;
data person_time_all;
set person_time_all;
CIN2_IRR=CIN2_IR/(ref_CIN2_IR);
low_IRR_CIN2=exp(log(CIN2_IRR)-1.96*sqrt(1/fu_cin2_sum+1/ref_cin2_sum));
Up_IRR_CIN2=exp(log(CIN2_IRR)+1.96*sqrt(1/fu_cin2_sum+1/ref_cin2_sum));
cancer_IRR=cancer_IR/(ref_cancer_IR);
low_IRR_cancer=exp(log(cancer_IRR)-1.96*sqrt(1/cancer_sum+1/ref_cancer_sum));
Up_IRR_cancer=exp(log(cancer_IRR)+1.96*sqrt(1/cancer_sum+1/ref_cancer_sum));
run;
/*too all groups*/

data person_time_negative;
set person_time_negative;
ref_cin2_ir=lag(cin2_ir); ref_cin2_sum=lag(fu_cin2_sum); ref_cancer_ir=lag(cancer_ir); ref_cancer_sum=lag(cancer_sum);
if char='ref' then do; ref_cin2_ir=cin2_ir; ref_cin2_sum=fu_cin2_sum; ref_cancer_ir=cancer_ir; ref_cancer_sum=cancer_sum;end;
ref_cin2_ir=lag(ref_cin2_ir); ref_cin2_sum=lag(ref_cin2_sum); ref_cancer_ir=lag(ref_cancer_ir); ref_cancer_sum=lag(ref_cancer_sum);
if char='ref' then do; ref_cin2_ir=cin2_ir; ref_cin2_sum=fu_cin2_sum; ref_cancer_ir=cancer_ir; ref_cancer_sum=cancer_sum;end;
run;
data person_time_negative;
set person_time_negative;
CIN2_IRR=CIN2_IR/(ref_CIN2_IR);
low_IRR_CIN2=exp(log(CIN2_IRR)-1.96*sqrt(1/fu_cin2_sum+1/ref_cin2_sum));
Up_IRR_CIN2=exp(log(CIN2_IRR)+1.96*sqrt(1/fu_cin2_sum+1/ref_cin2_sum));
cancer_IRR=cancer_IR/(ref_cancer_IR);
low_IRR_cancer=exp(log(cancer_IRR)-1.96*sqrt(1/cancer_sum+1/ref_cancer_sum));
Up_IRR_cancer=exp(log(cancer_IRR)+1.96*sqrt(1/cancer_sum+1/ref_cancer_sum));
run;

data person_time_positive;
set person_time_positive;
ref_cin2_ir=lag(cin2_ir); ref_cin2_sum=lag(fu_cin2_sum); ref_cancer_ir=lag(cancer_ir); ref_cancer_sum=lag(cancer_sum);
if char='ref' then do; ref_cin2_ir=cin2_ir; ref_cin2_sum=fu_cin2_sum; ref_cancer_ir=cancer_ir; ref_cancer_sum=cancer_sum;end;
ref_cin2_ir=lag(ref_cin2_ir); ref_cin2_sum=lag(ref_cin2_sum); ref_cancer_ir=lag(ref_cancer_ir); ref_cancer_sum=lag(ref_cancer_sum);
if char='ref' then do; ref_cin2_ir=cin2_ir; ref_cin2_sum=fu_cin2_sum; ref_cancer_ir=cancer_ir; ref_cancer_sum=cancer_sum;end;
run;
data person_time_positive;
set person_time_positive;
CIN2_IRR=CIN2_IR/(ref_CIN2_IR);
low_IRR_CIN2=exp(log(CIN2_IRR)-1.96*sqrt(1/fu_cin2_sum+1/ref_cin2_sum));
Up_IRR_CIN2=exp(log(CIN2_IRR)+1.96*sqrt(1/fu_cin2_sum+1/ref_cin2_sum));
cancer_IRR=cancer_IR/(ref_cancer_IR);
low_IRR_cancer=exp(log(cancer_IRR)-1.96*sqrt(1/cancer_sum+1/ref_cancer_sum));
Up_IRR_cancer=exp(log(cancer_IRR)+1.96*sqrt(1/cancer_sum+1/ref_cancer_sum));
run;
/*category outcome*/


data person_time_cat_all;
set person_time_cat_all;
ref_cin2_ir=lag(cin2_ir); ref_cin2_sum=lag(fu_cin2_sum); ref_cancer_ir=lag(cancer_ir); ref_cancer_sum=lag(cancer_sum);
if previous_pos_cat=0 then do; ref_cin2_ir=cin2_ir; ref_cin2_sum=fu_cin2_sum; ref_cancer_ir=cancer_ir; ref_cancer_sum=cancer_sum;end;
ref_cin2_ir=lag(ref_cin2_ir); ref_cin2_sum=lag(ref_cin2_sum); ref_cancer_ir=lag(ref_cancer_ir); ref_cancer_sum=lag(ref_cancer_sum);
if previous_pos_cat=0 then do; ref_cin2_ir=cin2_ir; ref_cin2_sum=fu_cin2_sum; ref_cancer_ir=cancer_ir; ref_cancer_sum=cancer_sum;end;
ref_cin2_ir=lag(ref_cin2_ir); ref_cin2_sum=lag(ref_cin2_sum); ref_cancer_ir=lag(ref_cancer_ir); ref_cancer_sum=lag(ref_cancer_sum);
if previous_pos_cat=0 then do; ref_cin2_ir=cin2_ir; ref_cin2_sum=fu_cin2_sum; ref_cancer_ir=cancer_ir; ref_cancer_sum=cancer_sum;end;
run;
data person_time_cat_all;
set person_time_cat_all;
CIN2_IRR=CIN2_IR/(ref_CIN2_IR);
low_IRR_CIN2=exp(log(CIN2_IRR)-1.96*sqrt(1/fu_cin2_sum+1/ref_cin2_sum));
Up_IRR_CIN2=exp(log(CIN2_IRR)+1.96*sqrt(1/fu_cin2_sum+1/ref_cin2_sum));
cancer_IRR=cancer_IR/(ref_cancer_IR);
low_IRR_cancer=exp(log(cancer_IRR)-1.96*sqrt(1/cancer_sum+1/ref_cancer_sum));
Up_IRR_cancer=exp(log(cancer_IRR)+1.96*sqrt(1/cancer_sum+1/ref_cancer_sum));
run;
/*negative or positive*/

data person_time_cat_negative;
set person_time_cat_negative;
ref_cin2_ir=lag(cin2_ir); ref_cin2_sum=lag(fu_cin2_sum); ref_cancer_ir=lag(cancer_ir); ref_cancer_sum=lag(cancer_sum);
if previous_pos_cat=0 then do; ref_cin2_ir=cin2_ir; ref_cin2_sum=fu_cin2_sum; ref_cancer_ir=cancer_ir; ref_cancer_sum=cancer_sum;end;
ref_cin2_ir=lag(ref_cin2_ir); ref_cin2_sum=lag(ref_cin2_sum); ref_cancer_ir=lag(ref_cancer_ir); ref_cancer_sum=lag(ref_cancer_sum);
if previous_pos_cat=0 then do; ref_cin2_ir=cin2_ir; ref_cin2_sum=fu_cin2_sum; ref_cancer_ir=cancer_ir; ref_cancer_sum=cancer_sum;end;
ref_cin2_ir=lag(ref_cin2_ir); ref_cin2_sum=lag(ref_cin2_sum); ref_cancer_ir=lag(ref_cancer_ir); ref_cancer_sum=lag(ref_cancer_sum);
if previous_pos_cat=0 then do; ref_cin2_ir=cin2_ir; ref_cin2_sum=fu_cin2_sum; ref_cancer_ir=cancer_ir; ref_cancer_sum=cancer_sum;end;
run;
data person_time_cat_negative;
set person_time_cat_negative;
CIN2_IRR=CIN2_IR/(ref_CIN2_IR);
low_IRR_CIN2=exp(log(CIN2_IRR)-1.96*sqrt(1/fu_cin2_sum+1/ref_cin2_sum));
Up_IRR_CIN2=exp(log(CIN2_IRR)+1.96*sqrt(1/fu_cin2_sum+1/ref_cin2_sum));
cancer_IRR=cancer_IR/(ref_cancer_IR);
low_IRR_cancer=exp(log(cancer_IRR)-1.96*sqrt(1/cancer_sum+1/ref_cancer_sum));
Up_IRR_cancer=exp(log(cancer_IRR)+1.96*sqrt(1/cancer_sum+1/ref_cancer_sum));
run;


data person_time_cat_positive;
set person_time_cat_positive;
ref_cin2_ir=lag(cin2_ir); ref_cin2_sum=lag(fu_cin2_sum); ref_cancer_ir=lag(cancer_ir); ref_cancer_sum=lag(cancer_sum);
if previous_pos_cat=0 then do; ref_cin2_ir=cin2_ir; ref_cin2_sum=fu_cin2_sum; ref_cancer_ir=cancer_ir; ref_cancer_sum=cancer_sum;end;
ref_cin2_ir=lag(ref_cin2_ir); ref_cin2_sum=lag(ref_cin2_sum); ref_cancer_ir=lag(ref_cancer_ir); ref_cancer_sum=lag(ref_cancer_sum);
if previous_pos_cat=0 then do; ref_cin2_ir=cin2_ir; ref_cin2_sum=fu_cin2_sum; ref_cancer_ir=cancer_ir; ref_cancer_sum=cancer_sum;end;
ref_cin2_ir=lag(ref_cin2_ir); ref_cin2_sum=lag(ref_cin2_sum); ref_cancer_ir=lag(ref_cancer_ir); ref_cancer_sum=lag(ref_cancer_sum);
if previous_pos_cat=0 then do; ref_cin2_ir=cin2_ir; ref_cin2_sum=fu_cin2_sum; ref_cancer_ir=cancer_ir; ref_cancer_sum=cancer_sum;end;
run;
data person_time_cat_positive;
set person_time_cat_positive;
CIN2_IRR=CIN2_IR/(ref_CIN2_IR);
low_IRR_CIN2=exp(log(CIN2_IRR)-1.96*sqrt(1/fu_cin2_sum+1/ref_cin2_sum));
Up_IRR_CIN2=exp(log(CIN2_IRR)+1.96*sqrt(1/fu_cin2_sum+1/ref_cin2_sum));
cancer_IRR=cancer_IR/(ref_cancer_IR);
low_IRR_cancer=exp(log(cancer_IRR)-1.96*sqrt(1/cancer_sum+1/ref_cancer_sum));
Up_IRR_cancer=exp(log(cancer_IRR)+1.96*sqrt(1/cancer_sum+1/ref_cancer_sum));
run;



/*outcome for long-term follow up of CIN2+ and cancer*/
data longoutcome;
set person_time_:;
run;

/*incidence rate*/
data longoutcome;
set longoutcome;
Low_cin2_ir=exp(log(CIN2_ir/100000)-1.96*(1/sqrt(fu_cin2_sum)))*100000;
High_cin2_ir=exp(log(CIN2_ir/100000)+1.96*(1/sqrt(fu_cin2_sum)))*100000;
Low_cancer_ir=exp(log(cancer_ir/100000)-1.96*(1/sqrt(cancer_sum)))*100000;
High_cancer_ir=exp(log(cancer_ir/100000)+1.96*(1/sqrt(cancer_sum)))*100000;
run;



/*organize the outcomes*/
data table_6;
set longoutcome;
format N_cin2 12. sur_time_cin2 $30. IR_CIN2 $30. IRR_CIN2 $30.  N_cancer 12. sur_time_cancer $30. IR_cancer $30. IRR_cancer $30.;
keep char HPVoutcome _FREQ_ N_cin2  sur_time_cin2  IR_CIN2  IRR_CIN2   N_cancer  sur_time_cancer  IR_cancer  IRR_cancer ;
N_cin2=fu_cin2_sum;
N_cancer=cancer_sum;
sur_time_cin2=cat(compress(put(sur_cin2_median,8.1)),' (',compress(put(sur_cin2_q1,8.1)),', ',compress(put(sur_cin2_q3,8.1)),')');
sur_time_cancer=cat(compress(put(sur_cancer_median,8.1)),' (',compress(put(sur_cancer_q1,8.1)),', ',compress(put(sur_cancer_q3,8.1)),')');
IR_cin2=cat(compress(put(cin2_ir,8.1)),' (',compress(put(low_cin2_ir,8.1)),', ',compress(put(high_cin2_ir,8.1)),')');
IR_cancer=cat(compress(put(cancer_ir,8.2)),' (',compress(put(low_cancer_ir,8.2)),', ',compress(put(high_cancer_ir,8.2)),')');
IRR_cin2=cat(compress(put(cin2_irr,8.1)),' (',compress(put(low_IRR_cin2,8.1)),', ',compress(put(up_IRR_CIN2,8.1)),')');
IRR_cancer=cat(compress(put(cancer_irr,8.1)),' (',compress(put(low_IRR_cancer,8.1)),', ',compress(put(up_IRR_cancer,8.1)),')');
run;

proc sort data=table_6;
by HPVoutcome;
run;

data table_6;
set table_6;
if char='ref' or index(char,'0')>0 then do; IRR_cin2=.;IRR_cancer=.;end;
run;
proc sort data=table_6;
by HPVoutcome char;
run;


title 'Table 6. long-term follow up of CIN2+ and cancer after first HPV screening test.';
ods rtf file="&mydir.Table5_IR_CIN2_ICC_&sysdate..rtf";
proc print data=table_6 noobs;run;
ods rtf close;
title;



/*calculate cumulative incidence of CIN2+*/

data pop_sur_time;
set pop_sur_time;
if cancer=. then cancer=0;
if fu_cin2=. then fu_cin2=0;
run;
/*all*/

ods exclude ProductLimitEstimates;
proc lifetest data=pop_sur_time nelson
plot=(survival) outsurv= cin2_sur_all ;
time sur_cin2*fu_cin2(0);
*strata even/test=all;
*ods output ProductLimitEstimates=can_sur_all;
run;



data cin2_sur_all;
set cin2_sur_all;
format char $30.;
if survival^=.;
char='all';
run;

proc sort data=cin2_sur_all nodupkey;
by  survival;
run;
proc sort data=cin2_sur_all nodupkey;
by  char;
run;
/*negative*/
ods exclude ProductLimitEstimates;
proc lifetest data=pop_sur_time nelson
plot=(survival) outsurv= cin2_sur_neg ;
time sur_cin2*fu_cin2(0);
where HPVrisk=0;
*strata even/test=all;
*ods output ProductLimitEstimates=can_sur_all;
run;

data cin2_sur_neg;
set cin2_sur_neg;
format char $30.;
if survival^=.;
char='neg';
run;

proc sort data=cin2_sur_neg nodupkey;
by  survival;
run;
proc sort data=cin2_sur_neg nodupkey;
by  char;
run;
/*positive*/
ods exclude ProductLimitEstimates;
proc lifetest data=pop_sur_time nelson
plot=(survival) outsurv= cin2_sur_pos ;
time sur_cin2*fu_cin2(0);
where HPVrisk>0;
*strata even/test=all;
*ods output ProductLimitEstimates=can_sur_all;
run;

data cin2_sur_pos;
set cin2_sur_pos;
format char $30.;
if survival^=.;
char='pos';
run;

proc sort data=cin2_sur_pos nodupkey;
by  survival;
run;
proc sort data=cin2_sur_pos nodupkey;
by  char;
run;
/*ref*/
/*all*/
ods exclude ProductLimitEstimates;
proc lifetest data=pop_sur_time nelson
plot=(survival) outsurv= cin2_sur_all_ref ;
time sur_cin2*fu_cin2(0);
where noscr_his=0 and previous_pos=0;
*strata even/test=all;
*ods output ProductLimitEstimates=can_sur_all;
run;

data cin2_sur_all_ref;
set cin2_sur_all_ref;
format char $30.;
if survival^=.;
char='all_ref';
run;

proc sort data=cin2_sur_all_ref nodupkey;
by  survival;
run;
proc sort data=cin2_sur_all_ref nodupkey;
by  char;
run;
/*negative*/
ods exclude ProductLimitEstimates;
proc lifetest data=pop_sur_time nelson
plot=(survival) outsurv= cin2_sur_neg_ref ;
time sur_cin2*fu_cin2(0);
where HPVrisk=0 and noscr_his=0 and previous_pos=0;
*strata even/test=all;
*ods output ProductLimitEstimates=can_sur_all;
run;

data cin2_sur_neg_ref;
set cin2_sur_neg_ref;
format char $30.;
if survival^=.;
char='neg_ref';
run;

proc sort data=cin2_sur_neg_ref nodupkey;
by  survival;
run;
proc sort data=cin2_sur_neg_ref nodupkey;
by  char;
run;
/*positive*/
ods exclude ProductLimitEstimates;
proc lifetest data=pop_sur_time nelson
plot=(survival) outsurv= cin2_sur_pos_ref ;
time sur_cin2*fu_cin2(0);
where HPVrisk>0 and noscr_his=0 and previous_pos=0;
*strata even/test=all;
*ods output ProductLimitEstimates=can_sur_all;
run;

data cin2_sur_pos_ref;
set cin2_sur_pos_ref;
format char $30.;
if survival^=.;
char='pos_ref';
run;

proc sort data=cin2_sur_pos_ref nodupkey;
by  survival;
run;
proc sort data=cin2_sur_pos_ref nodupkey;
by  char;
run;
/*noscr_his*/
/*all*/
ods exclude ProductLimitEstimates;
proc lifetest data=pop_sur_time nelson
plot=(survival) outsurv= cin2_sur_all_nohis ;
time sur_cin2*fu_cin2(0);
where noscr_his=1;
*strata even/test=all;
*ods output ProductLimitEstimates=can_sur_all;
run;

data cin2_sur_all_nohis;
set cin2_sur_all_nohis;
format char $30.;
if survival^=.;
char='all_nohis';
run;

proc sort data=cin2_sur_all_nohis nodupkey;
by  survival;
run;
proc sort data=cin2_sur_all_nohis nodupkey;
by  char;
run;
/*negative*/
ods exclude ProductLimitEstimates;
proc lifetest data=pop_sur_time nelson
plot=(survival) outsurv= cin2_sur_neg_nohis ;
time sur_cin2*fu_cin2(0);
where HPVrisk=0 and noscr_his=1;
*strata even/test=all;
*ods output ProductLimitEstimates=can_sur_all;
run;


data cin2_sur_neg_nohis;
set cin2_sur_neg_nohis;
format char $30.;
if survival^=.;
char='neg_nohis';
run;

proc sort data=cin2_sur_neg_nohis nodupkey;
by  survival;
run;
proc sort data=cin2_sur_neg_nohis nodupkey;
by  char;
run;
/*positive*/
ods exclude ProductLimitEstimates;
proc lifetest data=pop_sur_time nelson
plot=(survival) outsurv= cin2_sur_pos_nohis ;
time sur_cin2*fu_cin2(0);
where HPVrisk>0 and noscr_his=1;
*strata even/test=all;
*ods output ProductLimitEstimates=can_sur_all;
run;

data cin2_sur_pos_nohis;
set cin2_sur_pos_nohis;
format char $30.;
if survival^=.;
char='pos_nohis';
run;

proc sort data=cin2_sur_pos_nohis nodupkey;
by  survival;
run;
proc sort data=cin2_sur_pos_nohis nodupkey;
by  char;
run;
/*previous pos*/
/*all*/
ods exclude ProductLimitEstimates;
proc lifetest data=pop_sur_time nelson
plot=(survival) outsurv= cin2_sur_all_prepos ;
time sur_cin2*fu_cin2(0);
where previous_pos=1;
*strata previous_pos_cat/test=all;
*ods output ProductLimitEstimates=can_sur_all;
run;

data cin2_sur_all_prepos;
set cin2_sur_all_prepos;
format char $30.;
if survival^=.;
char='all_prepos';
run;

proc sort data=cin2_sur_all_prepos nodupkey;
by  survival;
run;
proc sort data=cin2_sur_all_prepos nodupkey;
by  char;
run;
/*negative*/
ods exclude ProductLimitEstimates;
proc lifetest data=pop_sur_time nelson
plot=(survival) outsurv= cin2_sur_neg_prepos ;
time sur_cin2*fu_cin2(0);
where HPVrisk=0 and previous_pos=1;
*strata even/test=all;
*ods output ProductLimitEstimates=can_sur_all;
run;

data cin2_sur_neg_prepos;
set cin2_sur_neg_prepos;
format char $30.;
if survival^=.;
char='neg_prepos';
run;

proc sort data=cin2_sur_neg_prepos nodupkey;
by  survival;
run;
proc sort data=cin2_sur_neg_prepos nodupkey;
by  char;
run;
/*positive*/
ods exclude ProductLimitEstimates;
proc lifetest data=pop_sur_time nelson
plot=(survival) outsurv= cin2_sur_pos_prepos ;
time sur_cin2*fu_cin2(0);
where HPVrisk>0 and previous_pos=1;
*strata even/test=all;
*ods output ProductLimitEstimates=can_sur_all;
run;

data cin2_sur_pos_prepos;
set cin2_sur_pos_prepos;
format char $30.;
if survival^=.;
char='pos_prepos';
run;

proc sort data=cin2_sur_pos_prepos nodupkey;
by  survival;
run;
proc sort data=cin2_sur_pos_prepos nodupkey;
by  char;
run;
/*cat*/
/*all*/
ods exclude ProductLimitEstimates;
proc lifetest data=pop_sur_time nelson
plot=(survival) outsurv= cin2_sur_all_cat ;
time sur_cin2*fu_cin2(0);
where previous_pos=1;
strata previous_pos_cat/test=all;
*ods output ProductLimitEstimates=can_sur_all;
run;

data cin2_sur_all_cat;
set cin2_sur_all_cat;
format char $30.;
if survival^=.;
char=cat('all_prepos-',previous_pos_cat);
run;

proc sort data=cin2_sur_all_cat nodupkey;
by  char survival;
run;
proc sort data=cin2_sur_all_cat nodupkey;
by  char;
run;
/*negative*/
ods exclude ProductLimitEstimates;
proc lifetest data=pop_sur_time nelson
plot=(survival) outsurv= cin2_sur_neg_cat ;
time sur_cin2*fu_cin2(0);
where HPVrisk=0 and previous_pos=1;
strata previous_pos_cat/test=all;
*ods output ProductLimitEstimates=can_sur_all;
run;

data cin2_sur_neg_cat;
set cin2_sur_neg_cat;
format char $30.;
if survival^=.;
char=cat('neg_prepos-',previous_pos_cat);
run;

proc sort data=cin2_sur_neg_cat nodupkey;
by  char survival;
run;
proc sort data=cin2_sur_neg_cat nodupkey;
by  char;
run;
/*positive*/
ods exclude ProductLimitEstimates;
proc lifetest data=pop_sur_time nelson
plot=(survival) outsurv= cin2_sur_pos_cat ;
time sur_cin2*fu_cin2(0);
where HPVrisk>0 and previous_pos=1;
strata previous_pos_cat/test=all;
*ods output ProductLimitEstimates=can_sur_all;
run;

data cin2_sur_pos_cat;
set cin2_sur_pos_cat;
format char $30.;
if survival^=.;
char=cat('pos_prepos-',previous_pos_cat);
run;

proc sort data=cin2_sur_pos_cat nodupkey;
by  char survival;
run;
proc sort data=cin2_sur_pos_cat nodupkey;
by  char;
run;



/*cancer*/
/*calculate cumulative incidence of cancer*/


/*all*/
ods exclude ProductLimitEstimates;
proc lifetest data=pop_sur_time nelson
plot=(survival) outsurv= cancer_sur_all ;
time sur_cancer*cancer(0);
*strata even/test=all;
*ods output ProductLimitEstimates=can_sur_all;
run;

data cancer_sur_all;
set cancer_sur_all;
format char $30.;
if survival^=.;
char='all';
run;

proc sort data=cancer_sur_all nodupkey;
by  survival;
run;
proc sort data=cancer_sur_all nodupkey;
by  char;
run;
/*negative*/
ods exclude ProductLimitEstimates;
proc lifetest data=pop_sur_time nelson
plot=(survival) outsurv= cancer_sur_neg ;
time sur_cancer*cancer(0);
where HPVrisk=0;
*strata even/test=all;
*ods output ProductLimitEstimates=can_sur_all;
run;

data cancer_sur_neg;
set cancer_sur_neg;
format char $30.;
if survival^=.;
char='neg';
run;

proc sort data=cancer_sur_neg nodupkey;
by  survival;
run;
proc sort data=cancer_sur_neg nodupkey;
by  char;
run;
/*positive*/
ods exclude ProductLimitEstimates;
proc lifetest data=pop_sur_time nelson
plot=(survival) outsurv= cancer_sur_pos ;
time sur_cancer*cancer(0);
where HPVrisk>0;
*strata even/test=all;
*ods output ProductLimitEstimates=can_sur_all;
run;

data cancer_sur_pos;
set cancer_sur_pos;
format char $30.;
if survival^=.;
char='pos';
run;

proc sort data=cancer_sur_pos nodupkey;
by  survival;
run;
proc sort data=cancer_sur_pos nodupkey;
by  char;
run;
/*ref*/
/*all*/
ods exclude ProductLimitEstimates;
proc lifetest data=pop_sur_time nelson
plot=(survival) outsurv= cancer_sur_all_ref ;
time sur_cancer*cancer(0);
where noscr_his=0 and previous_pos=0;
*strata even/test=all;
*ods output ProductLimitEstimates=can_sur_all;
run;

data cancer_sur_all_ref;
set cancer_sur_all_ref;
format char $30.;
if survival^=.;
char='all_ref';
run;

proc sort data=cancer_sur_all_ref nodupkey;
by  survival;
run;
proc sort data=cancer_sur_all_ref nodupkey;
by  char;
run;
/*negative*/
ods exclude ProductLimitEstimates;
proc lifetest data=pop_sur_time nelson
plot=(survival) outsurv= cancer_sur_neg_ref ;
time sur_cancer*cancer(0);
where HPVrisk=0 and noscr_his=0 and previous_pos=0;
*strata even/test=all;
*ods output ProductLimitEstimates=can_sur_all;
run;

data cancer_sur_neg_ref;
set cancer_sur_neg_ref;
format char $30.;
if survival^=.;
char='neg_ref';
run;

proc sort data=cancer_sur_neg_ref nodupkey;
by  survival;
run;
proc sort data=cancer_sur_neg_ref nodupkey;
by  char;
run;
/*positive*/
ods exclude ProductLimitEstimates;
proc lifetest data=pop_sur_time nelson
plot=(survival) outsurv= cancer_sur_pos_ref ;
time sur_cancer*cancer(0);
where HPVrisk>0 and noscr_his=0 and previous_pos=0;
*strata even/test=all;
*ods output ProductLimitEstimates=can_sur_all;
run;

data cancer_sur_pos_ref;
set cancer_sur_pos_ref;
format char $30.;
if survival^=.;
char='pos_ref';
run;

proc sort data=cancer_sur_pos_ref nodupkey;
by  survival;
run;
proc sort data=cancer_sur_pos_ref nodupkey;
by  char;
run;
/*noscr_his*/
/*all*/
ods exclude ProductLimitEstimates;
proc lifetest data=pop_sur_time nelson
plot=(survival) outsurv= cancer_sur_all_nohis ;
time sur_cancer*cancer(0);
where noscr_his=1;
*strata even/test=all;
*ods output ProductLimitEstimates=can_sur_all;
run;

data cancer_sur_all_nohis;
set cancer_sur_all_nohis;
format char $30.;
if survival^=.;
char='all_nohis';
run;

proc sort data=cancer_sur_all_nohis nodupkey;
by  survival;
run;
proc sort data=cancer_sur_all_nohis nodupkey;
by  char;
run;
/*negative*/
ods exclude ProductLimitEstimates;
proc lifetest data=pop_sur_time nelson
plot=(survival) outsurv= cancer_sur_neg_nohis ;
time sur_cancer*cancer(0);
where HPVrisk=0 and noscr_his=1;
*strata even/test=all;
*ods output ProductLimitEstimates=can_sur_all;
run;

data cancer_sur_neg_nohis;
set cancer_sur_neg_nohis;
format char $30.;
if survival^=.;
char='neg_nohis';
run;

proc sort data=cancer_sur_neg_nohis nodupkey;
by  survival;
run;
proc sort data=cancer_sur_neg_nohis nodupkey;
by  char;
run;
/*positive*/
ods exclude ProductLimitEstimates;
proc lifetest data=pop_sur_time nelson
plot=(survival) outsurv= cancer_sur_pos_nohis ;
time sur_cancer*cancer(0);
where HPVrisk>0 and noscr_his=1;
*strata even/test=all;
*ods output ProductLimitEstimates=can_sur_all;
run;

data cancer_sur_pos_nohis;
set cancer_sur_pos_nohis;
format char $30.;
if survival^=.;
char='pos_nohis';
run;

proc sort data=cancer_sur_pos_nohis nodupkey;
by  survival;
run;
proc sort data=cancer_sur_pos_nohis nodupkey;
by  char;
run;
/*previous pos*/
/*all*/
ods exclude ProductLimitEstimates;
proc lifetest data=pop_sur_time nelson
plot=(survival) outsurv= cancer_sur_all_prepos ;
time sur_cancer*cancer(0);
where previous_pos=1;
*strata previous_pos_cat/test=all;
*ods output ProductLimitEstimates=can_sur_all;
run;

data cancer_sur_all_prepos;
set cancer_sur_all_prepos;
format char $30.;
if survival^=.;
char='all_prepos';
run;

proc sort data=cancer_sur_all_prepos nodupkey;
by  survival;
run;
proc sort data=cancer_sur_all_prepos nodupkey;
by  char;
run;
/*negative*/
ods exclude ProductLimitEstimates;
proc lifetest data=pop_sur_time nelson
plot=(survival) outsurv= cancer_sur_neg_prepos ;
time sur_cancer*cancer(0);
where HPVrisk=0 and previous_pos=1;
*strata even/test=all;
*ods output ProductLimitEstimates=can_sur_all;
run;

data cancer_sur_neg_prepos;
set cancer_sur_neg_prepos;
format char $30.;
if survival^=.;
char='neg_prepos';
run;

proc sort data=cancer_sur_neg_prepos nodupkey;
by  survival;
run;
proc sort data=cancer_sur_neg_prepos nodupkey;
by  char;
run;
/*positive*/
ods exclude ProductLimitEstimates;
proc lifetest data=pop_sur_time nelson
plot=(survival) outsurv= cancer_sur_pos_prepos ;
time sur_cancer*cancer(0);
where HPVrisk>0 and previous_pos=1;
*strata even/test=all;
*ods output ProductLimitEstimates=can_sur_all;
run;

data cancer_sur_pos_prepos;
set cancer_sur_pos_prepos;
format char $30.;
if survival^=.;
char='pos_prepos';
run;

proc sort data=cancer_sur_pos_prepos nodupkey;
by  survival;
run;
proc sort data=cancer_sur_pos_prepos nodupkey;
by  char;
run;
/*cat*/
/*all*/
ods exclude ProductLimitEstimates;
proc lifetest data=pop_sur_time nelson
plot=(survival) outsurv= cancer_sur_all_cat ;
time sur_cancer*cancer(0);
where previous_pos=1;
strata previous_pos_cat/test=all;
*ods output ProductLimitEstimates=can_sur_all;
run;

data cancer_sur_all_cat;
set cancer_sur_all_cat;
format char $30.;
if survival^=.;
char=cat('all_prepos-',previous_pos_cat);
run;

proc sort data=cancer_sur_all_cat nodupkey;
by  char survival;
run;
proc sort data=cancer_sur_all_cat nodupkey;
by  char;
run;
/*negative*/
ods exclude ProductLimitEstimates;
proc lifetest data=pop_sur_time nelson
plot=(survival) outsurv= cancer_sur_neg_cat ;
time sur_cancer*cancer(0);
where HPVrisk=0 and previous_pos=1;
strata previous_pos_cat/test=all;
*ods output ProductLimitEstimates=can_sur_all;
run;

data cancer_sur_neg_cat;
set cancer_sur_neg_cat;
format char $30.;
if survival^=.;
char=cat('neg_prepos-',previous_pos_cat);
run;

proc sort data=cancer_sur_neg_cat nodupkey;
by  char survival;
run;
proc sort data=cancer_sur_neg_cat nodupkey;
by  char;
run;
/*positive*/
ods exclude ProductLimitEstimates;
proc lifetest data=pop_sur_time nelson
plot=(survival) outsurv= cancer_sur_pos_cat ;
time sur_cancer*cancer(0);
where HPVrisk>0 and previous_pos=1;
strata previous_pos_cat/test=all;
*ods output ProductLimitEstimates=can_sur_all;
run;

data cancer_sur_pos_cat;
set cancer_sur_pos_cat;
format char $30.;
if survival^=.;
char=cat('pos_prepos-',previous_pos_cat);
run;

proc sort data=cancer_sur_pos_cat nodupkey;
by  char survival;
run;
proc sort data=cancer_sur_pos_cat nodupkey;
by  char;
run;

/*arrange*/
data cin2_sur_table;
set cin2_sur_:;
run;

data sur_out_cin2;
set cin2_sur_table;
cum_inc_CIN2=1-survival;
high_cin2=1-sdf_lcl;
low_cin2=1-sdf_ucl;
keep char cum_inc_cin2 high_cin2 low_cin2;
run;

data cancer_sur_table;
set cancer_sur_:;
run;

data sur_out_cancer;
set cancer_sur_table;
cum_inc_cancer=1-survival;
high_cancer=1-sdf_lcl;
low_cancer=1-sdf_ucl;
keep char cum_inc_cancer high_cancer low_cancer;
run;

/*arrange*/
data table_7a ;
set sur_out_cin2;
format cum_inc_95_cin2 $30.;
cum_inc_95_cin2=cat(compress(put(cum_inc_cin2*100,8.1)),'% (',compress(put(low_cin2*100,8.1)),'%, ',compress(put(high_cin2*100,8.1)),'%)');
keep char cum_inc_95_cin2;
run;
data table_7b ;
set sur_out_cancer;
format cum_inc_95_cancer $30.;
cum_inc_95_cancer=cat(compress(put(cum_inc_cancer*100,8.3)),'% (',compress(put(low_cancer*100,8.3)),'%, ',compress(put(high_cancer*100,8.3)),'%)');
keep char cum_inc_95_cancer;
run;

proc sort data=table_7a;
by char;
proc sort data=table_7b;
by char;
run;
data table_7;
merge table_7a table_7b;
by char;
run;



title 'Table 6. Cumulative incidence of CIN2+ and ICC ';
ods rtf file="&mydir.table6_CI_CIN_ICC_&sysdate..rtf";
proc print data=table_7 noobs;run;
ods rtf close;
title;





/*export dataset for stata*/
data pop_km;
set pop_sur_time;
keep person_id HPVrisk previous_pos previous_pos_cat noscr_his fu_cin2 sur_cin2 cancer sur_cancer hpv_age;
run;

data pop_km;
set pop_km;
format level $6.;
if previous_pos=0 and noscr_his=0 then level='A';
else if previous_pos_cat=1 then level='B1';
else if previous_pos_cat=2 then level='B2';
else if previous_pos_cat=3 then level='B3';
run;

proc freq data=pop_km;
tables level;
run;


data pop_km;
set pop_km;
if sur_cin2=0 then sur_cin2=0.000000000001;
if sur_cancer=0 then sur_can=0.000000000001;
run;

proc export data=pop_km outfile="&mydir.pop_km_260604_notruncate.dta" replace dbms=stata replace; run; 


/*log rank*/

ods exclude ProductLimitEstimates;
proc lifetest data=pop_km nelson
plot=(survival) outsurv= logrank ;
time sur_cin2*fu_cin2(0);
where  level in ('A' 'B2') and hpvrisk=0;
strata level/test=all;
*ods output ProductLimitEstimates=can_sur_all;
run;

















/*cox compare*/
/*no scr to ref*/
/*need to differentiate positive and negative, and the*/

data pop_sur_time_cox;
set pop_sur_time;
rename pad_cin2=cin2 previous_pos=pos previous_pos_cat=pos_cat;
run;
options mprint;
options mlogic;

%macro cox(HPVDIAG,outcome,exposure);
proc phreg data=pop_sur_time_cox simple;
class &exposure. (param=ref ref="0") ;
model  sur_&outcome.*&outcome.(0)=&exposure. /ties=exact rl=pl type3(all);
where
	%if &exposure.=noscr_his %then %do; pos=0 %end;
	%else %do; noscr_his=0 %end; 
	%if &HPVDIAG.=POS %then %do; and HPVrisk>0 %end;
	%else %if &HPVDIAG.=NEG %then %do; and HPVrisk=0 %end;
;
ods output ParameterEstimates=hazard_org_&outcome._&exposure._&HPVDIAG. ;
run;

data hazard_org_&outcome._&exposure._&HPVDIAG.;
set hazard_org_&outcome._&exposure._&HPVDIAG.;
format char $30.;
char=cat(compress(label),'-',"&HPVDIAG.");
HRLowerPLCL=exp(estimate-1.96*stderr);
HRUpperPLCL=exp(estimate+1.96*stderr);
run;


proc phreg data=pop_sur_time_cox simple;
class &exposure. (param=ref ref="0") 
	%if &outcome.=cin2 %then %do; sample_year %end;
	%else %do; sample_period %end;
birthco ;
model  sur_&outcome.*&outcome.(0)=&exposure. 
	%if &outcome.=cin2 %then %do; sample_year %end;
	%else %do; sample_period %end;
  birthco/ties=exact rl=pl type3(all);
where 
	%if &exposure.=noscr_his %then %do; pos=0 %end;
	%else %do; noscr_his=0 %end;
	%if &HPVDIAG.=POS %then %do; and HPVrisk>0 %end;
	%else %if &HPVDIAG.=NEG %then %do; and HPVrisk=0 %end;
;
ods output ParameterEstimates=hazard_adj_&outcome._&exposure._&HPVDIAG.;
run;
data hazard_adj_&outcome._&exposure._&HPVDIAG.;
set hazard_adj_&outcome._&exposure._&HPVDIAG.;
format char $30.;
char=cat(compress(label),'-',"&HPVDIAG.",);
if parameter="&exposure.";
HRLowerPLCL=exp(estimate-1.96*stderr);
HRUpperPLCL=exp(estimate+1.96*stderr);
run;
%mend;
%cox(NEG,cin2,noscr_his);
%cox(POS,cin2,noscr_his);
%cox(ALL,cin2,noscr_his);
%cox(NEG,cancer,noscr_his);
%cox(POS,cancer,noscr_his);
%cox(ALL,cancer,noscr_his);

%cox(NEG,cin2,pos);
%cox(POS,cin2,pos);
%cox(ALL,cin2,pos);
%cox(NEG,cancer,pos);
%cox(POS,cancer,pos);
%cox(ALL,cancer,pos);


%cox(NEG,cin2,pos_cat);
%cox(POS,cin2,pos_cat);
%cox(ALL,cin2,pos_cat);
%cox(NEG,cancer,pos_cat);
%cox(POS,cancer,pos_cat);
%cox(ALL,cancer,pos_cat);


%macro coxcat(HPVDIAG,outcome,level);
proc phreg data=pop_sur_time_cox simple;
class pos_cat (param=ref ref="0") ;
model  sur_&outcome.*&outcome.(0)=pos_cat /ties=exact rl=pl type3(all);
where noscr_his=0 
	%if &HPVDIAG.=POS %then %do; and HPVrisk>0 %end;
	%else %if &HPVDIAG.=NEG %then %do; and HPVrisk=0 %end;
	and pos_cat in (0 &level.)
;
ods output ParameterEstimates=hazard_org_&outcome._&level._&HPVDIAG. ;
run;

data hazard_org_&outcome._&level._&HPVDIAG.;
set hazard_org_&outcome._&level._&HPVDIAG.;
format char $30.;
char=cat(compress(label),'-',"&HPVDIAG.");
HRLowerPLCL=exp(estimate-1.96*stderr);
HRUpperPLCL=exp(estimate+1.96*stderr);
run;


proc phreg data=pop_sur_time_cox simple;
class pos_cat (param=ref ref="0") 
	%if &outcome.=cin2 %then %do; sample_year %end;
	%else %do; sample_period %end;
birthco ;
model  sur_&outcome.*&outcome.(0)=pos_cat 
	%if &outcome.=cin2 %then %do; sample_year %end;
	%else %do; sample_period %end;
  birthco/ties=exact rl=pl type3(all);
where noscr_his=0 
	%if &HPVDIAG.=POS %then %do; and HPVrisk>0 %end;
	%else %if &HPVDIAG.=NEG %then %do; and HPVrisk=0 %end;
	and pos_cat in (0 &level.)
;
ods output ParameterEstimates=hazard_adj_&outcome._&level._&HPVDIAG.;
run;
data hazard_adj_&outcome._&level._&HPVDIAG.;
set hazard_adj_&outcome._&level._&HPVDIAG.;
format char $30.;
char=cat(compress(label),'-',"&HPVDIAG.");
if parameter='pos_cat';
HRLowerPLCL=exp(estimate-1.96*stderr);
HRUpperPLCL=exp(estimate+1.96*stderr);
run;
%mend;
%coxcat(NEG,cin2,1);
%coxcat(NEG,cin2,2);
%coxcat(NEG,cin2,3);
%coxcat(POS,cin2,1);
%coxcat(POS,cin2,2);
%coxcat(POS,cin2,3);
%coxcat(ALL,cin2,1);
%coxcat(ALL,cin2,2);
%coxcat(ALL,cin2,3);
%coxcat(NEG,cancer,1);
%coxcat(NEG,cancer,2);
%coxcat(NEG,cancer,3);
%coxcat(POS,cancer,1);
%coxcat(POS,cancer,2);
%coxcat(POS,cancer,3);
%coxcat(ALL,cancer,1);
%coxcat(ALL,cancer,2);
%coxcat(ALL,cancer,3);




proc freq data=	pop_sur_time;
table previous_pos*cancer;
where HPVrisk=0 and noscr_his=0;
run;
 
/*keep pop_sur_time as reference*/
data preabn.pop_sur_time_newpop;
set pop_sur_time;
run;


/*previous_pos to ref*/
proc phreg data=pop_sur_time simple;
class previous_pos (param=ref ref="0") ;
model  sur_cin2*pad_cin2(0)=previous_pos /ties=exact rl=pl type3(all);
where noscr_his=0 and HPVrisk>0;
ods output ParameterEstimates=hazard_cin2_noscr_adj ;
run;

proc phreg data=pop_sur_time simple;
class previous_pos (param=ref ref="0") sample_year  birthco ;
model  sur_cin2*pad_cin2(0)=previous_pos sample_year  birthco/ties=exact rl=pl type3(all);
*assess ph/crpanel resample;
where noscr_his=0 ;
ods output ParameterEstimates=hazard_cin2_noscr_adj ;
run;

/*previous_pos_cat to ref*/
proc phreg data=pop_sur_time simple;
class previous_pos_cat (param=ref ref="0") ;
model  sur_cin2*pad_cin2(0)=previous_pos_cat /ties=exact rl=pl type3(all);
assess ph/crpanel resample;
where noscr_his=0 and previous_pos_cat in (0 3) and HPVrisk>0;
ods output ParameterEstimates=hazard_cin2_noscr_adj ProportionalHazardsSupTest=ppt_cin2_noscr_adj;
run;

proc phreg data=pop_sur_time simple;
class previous_pos_cat (param=ref ref="0") sample_year  birthco ;
model  sur_cin2*pad_cin2(0)=previous_pos_cat sample_year  birthco/ties=exact rl=pl type3(all);
*assess ph/crpanel resample;
where noscr_his=0 and previous_pos_cat in (0 3) ;
ods output ParameterEstimates=hazard_cin2_noscr_adj ProportionalHazardsSupTest=ppt_cin2_noscr_adj;
run;
/*cancer*/

/*correction only compare the cumulative incidence use proc life*/


/*export data for survival analysis*/

data analysis_km;
set preabn.longterm_250123;
if previous_pos=. and noscr_his=0 then ref=1;
keep person_id HPVrisk noscr_his previous_pos previous_pos_cat age_group ref sur_cin2 pad_cin2 sur_cancer cancer birth_yr sample_year;
run;

data analysis_km;
set analysis_km;
if noscr_his=1 then noscr=1;
else if ref=1 then noscr=0;
if previous_pos=1 then pos_cat=previous_pos_cat;
else if ref=1 then pos_cat=0;
rename pad_cin2=cin2;
if HPVrisk=0 then HPV=0;
else HPV=1;
if cancer=. then cancer=0;
if pad_cin2=. then pad_cin2=0;
if previous_pos=. then previous_pos=0;
run;

proc freq data=analysis_km;
table noscr*pos_cat ref HPVrisk*HPV/missing;
run;

data analysis_km;
set analysis_km;
drop HPVrisk noscr_his previous_pos_cat;
run;


/*export data try to analysis with STATA*/

data analysis_stata;
set preabn.longterm_newpop53_250304;
keep person_id first_hpv_date pad_cin2 cancer eof_cin2 eof_cancer birth_yr sample_year year_cin2 year_ICC previous_pos previous_pos_cat ref HPV age_group sur_cancer sur_cin2;
if previous_pos=. and noscr_his=0 then ref=1;
if pad_cin2=. then pad_cin2=0;
if cancer=. then cancer=0;
year_cin2=year(eof_cin2);
year_ICC=year(eof_cancer);
if previous_pos=. and ref=1 then previous_pos=0;
if previous_pos_cat=. and ref=1 then previous_pos_cat=0;
if HPVrisk=0 then HPV=0;
if HPVrisk>0 then HPV=1;
run;

/*extract birth day for the population*/

data birth_date;
set socmob4.nkc_person_2022;
format bday yymmdd10.;
birth_year=input(substr(birth_date,1,4),8.);
birth_month=input(substr(birth_date,6,2),8.);
birth_day=input(substr(birth_date,9,2),8.);
bday=mdy(birth_month,birth_day, birth_year);
if valid_pnr=1;
run;
data birth_date;
set birth_date;
keep person_id bday;
run;
proc sort data=birth_date nodupkey out=nodup_bday;
by person_id;
run;
data multi;
set birth_date;
by person_id;
if first.person_id and last.person_id then delete;
run;
data analysis_stata;
merge analysis_stata(in=a) nodup_bday;
by person_id;
if a;
run;


proc export data=analysis_stata outfile="&mydir.stata_250304.xlsx" dbms=xlsx replace; sheet='all'; run;



data pop_sur_time;
set pop_sur_time;
keep person_id first_hpv_date pad_cin2 cancer eof_cin2 eof_cancer birth_yr sample_year year_cin2 year_ICC previous_pos previous_pos_cat ref HPV age_group sur_cancer sur_cin2;
if previous_pos=. and noscr_his=0 then ref=1;
if pad_cin2=. then pad_cin2=0;
if cancer=. then cancer=0;
year_cin2=year(eof_cin2);
year_ICC=year(eof_cancer);
if previous_pos=. and ref=1 then previous_pos=0;
if previous_pos_cat=. and ref=1 then previous_pos_cat=0;
if HPVrisk=0 then HPV=0;
if HPVrisk>0 then HPV=1;
run;

ods exclude ProductLimitEstimates;
proc lifetest data=pop_sur_time nelson plot=(survival)  ;
time sur_cin2*pad_cin2(0);
strata previous_pos_cat/test=all;
where HPV=0  and previous_pos_cat in (0 1) ;
*ods output ProductLimitEstimates=can_sur_all;
run;

/*newpop53 cancer p=0.0009/0.1076 all   not category 0.0150/0.1628
cin2= 0.1051/0.1754 all not cat 0.2591/0.5299
0-1
cin2=0.0150/0.0494
can=0.0197/0.1301
0-2
cin2=0.8998/0.3847
can=0.0005/0.0425
0-3
cin2=0.7957/0.6888
can=0.7397/0.7906
*/


proc freq data=analysis_stata;
tables cancer;
run;

data analysis_stata;
set analysis_stata;
hpv_age=year (first_hpv_date)- birth_yr;
run;
proc freq data=analysis_stata;
tables hpv_age*previous_pos hpv_age*previous_pos_cat;
run;
