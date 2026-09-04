/* Data created.: pop_index_&date */
proc datasets library=work kill;quit;

libname NCSR odbc noprompt="dsn=kosmos; database=NCSR" schema=NCSR;
libname V_NCSR odbc noprompt="dsn=kosmos; database=NCSR" schema=V_NCSR ;
libname PreAbn 'P:\ACCES\ACCES_Research\Qingyun\previous abnormality\Data';
%let mydir=P:\ACCES\ACCES_Research\Qingyun\previous abnormality\Doc\NKCx_260604_exnsh\;


/*time split in the survival data*/

data analysis;
set preabn.NKCx_nohys_longterm_fu;
keep person_id HPVrisk noscr_his previous_pos previous_pos_cat cancer sur_cancer;
where noscr_his=0;
run;

/*0-<3,3-<6,>=6*/
data time_split;
set analysis;
if .<sur_cancer<3 then do;
sur_cancer_3=sur_cancer;
sur_cancer_6=.;
sur_cancer_over=.;
end;
else if 3<=sur_cancer<6 then do;
sur_cancer_3=3;
sur_cancer_6=sur_cancer-3;
sur_cancer_over=.;
end;
else if sur_cancer>=6 then do;
sur_cancer_3=3;
sur_cancer_6=3;
sur_cancer_over=sur_cancer-6;
end;
if .<sur_cancer<3 and cancer=1 then do;
cancer_3=1;
cancer_6=0;
cancer_over=0;
end;
else if 3<=sur_cancer<6 and cancer=1 then do;
cancer_3=0;
cancer_6=1;
cancer_over=0;
end;
else if sur_cancer>=6 and cancer=1 then do;
cancer_3=0;
cancer_6=0;
cancer_over=1;
end;
if cancer_3=. then cancer_3=0;
if cancer_6=.  then cancer_6=0;
if cancer_over=. then cancer_over=0;
run;

data time_split;
set time_split;
if HPVrisk=0 then HPV=0;
else HPV=1;
if previous_pos=0 and noscr_his=0 then ref=1;
else ref=0;
where sur_cancer^=.;
run;

data time_split;
set time_split;
if sur_cancer_3^=. then sur_3=1;else sur_3=0;
if sur_cancer_6^=. then sur_6=1;else sur_6=0;
if sur_cancer_over^=. then sur_over=1;else sur_over=0;
run;


proc summary data=time_split;
var sur_cancer sur_cancer_3 sur_cancer_6 sur_cancer_over cancer_3 cancer_6 cancer_over sur_3 sur_6 sur_over;
class ref previous_pos previous_pos_cat noscr_his HPV;
ways 1;
output out=All_split sum()=/autoname;
run;
proc summary data=time_split;
var sur_cancer sur_cancer_3 sur_cancer_6 sur_cancer_over cancer_3 cancer_6 cancer_over sur_3 sur_6 sur_over;
output out=All_pop_split sum()=/autoname;
run;
data all_split; 
set All_pop_split all_split;
run;

data All_IR;
set All_split;
IR_3=(cancer_3_sum/sur_cancer_3_sum)*100000;
Low_ir_3=exp(log(IR_3/100000)-1.96*(1/sqrt(cancer_3_sum)))*100000;
High_ir_3=exp(log(IR_3/100000)+1.96*(1/sqrt(cancer_3_sum)))*100000;

IR_6=(cancer_6_sum/sur_cancer_6_sum)*100000;
Low_ir_6=exp(log(IR_6/100000)-1.96*(1/sqrt(cancer_6_sum)))*100000;
High_ir_6=exp(log(IR_6/100000)+1.96*(1/sqrt(cancer_6_sum)))*100000;

IR_over=(cancer_over_sum/sur_cancer_over_sum)*100000;
Low_ir_over=exp(log(IR_over/100000)-1.96*(1/sqrt(cancer_over_sum)))*100000;
High_ir_over=exp(log(IR_over/100000)+1.96*(1/sqrt(cancer_over_sum)))*100000;
run;

proc summary data=time_split;
var sur_cancer sur_cancer_3 sur_cancer_6 sur_cancer_over cancer_3 cancer_6 cancer_over sur_3 sur_6 sur_over;
class ref previous_pos previous_pos_cat noscr_his;
ways 1;
where HPV=0;
output out=Neg_split sum()=/autoname;
run;

data Neg_IR;
set Neg_split;
IR_3=(cancer_3_sum/sur_cancer_3_sum)*100000;
Low_ir_3=exp(log(IR_3/100000)-1.96*(1/sqrt(cancer_3_sum)))*100000;
High_ir_3=exp(log(IR_3/100000)+1.96*(1/sqrt(cancer_3_sum)))*100000;

IR_6=(cancer_6_sum/sur_cancer_6_sum)*100000;
Low_ir_6=exp(log(IR_6/100000)-1.96*(1/sqrt(cancer_6_sum)))*100000;
High_ir_6=exp(log(IR_6/100000)+1.96*(1/sqrt(cancer_6_sum)))*100000;

IR_over=(cancer_over_sum/sur_cancer_over_sum)*100000;
Low_ir_over=exp(log(IR_over/100000)-1.96*(1/sqrt(cancer_over_sum)))*100000;
High_ir_over=exp(log(IR_over/100000)+1.96*(1/sqrt(cancer_over_sum)))*100000;
run;

proc summary data=time_split;
var sur_cancer sur_cancer_3 sur_cancer_6 sur_cancer_over cancer_3 cancer_6 cancer_over sur_3 sur_6 sur_over;
class ref previous_pos previous_pos_cat noscr_his;
ways 1;
where HPV=1;
output out=Pos_split sum()=/autoname;
run;

data Pos_IR;
set Pos_split;
IR_3=(cancer_3_sum/sur_cancer_3_sum)*100000;
Low_ir_3=exp(log(IR_3/100000)-1.96*(1/sqrt(cancer_3_sum)))*100000;
High_ir_3=exp(log(IR_3/100000)+1.96*(1/sqrt(cancer_3_sum)))*100000;

IR_6=(cancer_6_sum/sur_cancer_6_sum)*100000;
Low_ir_6=exp(log(IR_6/100000)-1.96*(1/sqrt(cancer_6_sum)))*100000;
High_ir_6=exp(log(IR_6/100000)+1.96*(1/sqrt(cancer_6_sum)))*100000;

IR_over=(cancer_over_sum/sur_cancer_over_sum)*100000;
Low_ir_over=exp(log(IR_over/100000)-1.96*(1/sqrt(cancer_over_sum)))*100000;
High_ir_over=exp(log(IR_over/100000)+1.96*(1/sqrt(cancer_over_sum)))*100000;
run;


/*organize table*/
data table_a;
format  HPVoutcome $30. char $30.;
set all_ir;
HPVoutcome="All";
if _TYPE_=0 then char="All";
if HPV=0 then do; HPVoutcome="Negative"; char="All";end;
if HPV=1 then do; HPVoutcome="Positive"; char="All";end;
if noscr_his=1 then char="Noscr";
if previous_pos_cat^=. then char=cat("Previous_pos-",previous_pos_cat);
if previous_pos_cat=0 then delete;
if noscr_his=0 then delete;
if previous_pos=0 then delete;
if ref=0 then delete;
if previous_pos=1 then char=("Previous_pos");
if ref=1 then char=("Ref");
rename _freq_=N;
ICC_IR_3=cat(compress(put(IR_3,8.1)),' (',compress(put(low_ir_3,8.1)),', ',compress(put(high_ir_3,8.1)),')');
ICC_IR_6=cat(compress(put(IR_6,8.1)),' (',compress(put(low_ir_6,8.1)),', ',compress(put(high_ir_6,8.1)),')');
ICC_IR_over=cat(compress(put(IR_over,8.1)),' (',compress(put(low_ir_over,8.1)),', ',compress(put(high_ir_over,8.1)),')');
drop ref previous_pos previous_pos_cat noscr_his HPV _TYPE_ IR_3 low_ir_3 high_ir_3 IR_6 low_ir_6 high_ir_6 IR_over low_ir_over high_ir_over;
run;


/*organize table*/
data table_b;
format  HPVoutcome $30. char $30.;
set neg_ir;
HPVoutcome="Negative";
if _TYPE_=0 then char="All";
if HPV=0 then char="Negative";
if HPV=1 then char="Positive";
if noscr_his=1 then char="Noscr";
if previous_pos_cat^=. then char=cat("Previous_pos-",previous_pos_cat);
if previous_pos_cat=0 then delete;
if noscr_his=0 then delete;
if previous_pos=0 then delete;
if ref=0 then delete;
if previous_pos=1 then char=("Previous_pos");
if ref=1 then char=("Ref");
rename _freq_=N;
ICC_IR_3=cat(compress(put(IR_3,8.1)),' (',compress(put(low_ir_3,8.1)),', ',compress(put(high_ir_3,8.1)),')');
ICC_IR_6=cat(compress(put(IR_6,8.1)),' (',compress(put(low_ir_6,8.1)),', ',compress(put(high_ir_6,8.1)),')');
ICC_IR_over=cat(compress(put(IR_over,8.1)),' (',compress(put(low_ir_over,8.1)),', ',compress(put(high_ir_over,8.1)),')');
drop ref previous_pos previous_pos_cat noscr_his HPV _TYPE_ IR_3 low_ir_3 high_ir_3 IR_6 low_ir_6 high_ir_6 IR_over low_ir_over high_ir_over;
run;


/*organize table*/
data table_c;
format  HPVoutcome $30. char $30.;
set pos_ir;
HPVoutcome="Positive";
if _TYPE_=0 then char="All";
if HPV=0 then char="Negative";
if HPV=1 then char="Positive";
if noscr_his=1 then char="Noscr";
if previous_pos_cat^=. then char=cat("Previous_pos-",previous_pos_cat);
if previous_pos_cat=0 then delete;
if noscr_his=0 then delete;
if previous_pos=0 then delete;
if ref=0 then delete;
if previous_pos=1 then char=("Previous_pos");
if ref=1 then char=("Ref");
rename _freq_=N;
ICC_IR_3=cat(compress(put(IR_3,8.1)),' (',compress(put(low_ir_3,8.1)),', ',compress(put(high_ir_3,8.1)),')');
ICC_IR_6=cat(compress(put(IR_6,8.1)),' (',compress(put(low_ir_6,8.1)),', ',compress(put(high_ir_6,8.1)),')');
ICC_IR_over=cat(compress(put(IR_over,8.1)),' (',compress(put(low_ir_over,8.1)),', ',compress(put(high_ir_over,8.1)),')');
drop ref previous_pos previous_pos_cat noscr_his HPV _TYPE_ IR_3 low_ir_3 high_ir_3 IR_6 low_ir_6 high_ir_6 IR_over low_ir_over high_ir_over;
run;

data table;
set table_:;
run;

proc sort data=table;
by HPVoutcome char;
run;

proc export data=table
    outfile="&mydir.pre_ abn_chr_table_ICC_IR_time&sysdate..xlsx"
    dbms=xlsx
    replace;
sheet="Sheet1";
run;

title 'Table split IR for ICC';
ods rtf file="&mydir.pre_ abn_chr_table_ICC_IR_time&sysdate..rtf";
proc print data=table noobs;run;
title '';
ods rtf close;
/*export results*/
