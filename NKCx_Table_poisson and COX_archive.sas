/*compare the results of poisson and cox*/
/*poisson regression*/
/* ph assumption was tested using STATA due to faster time split*/
/*run this code immediately after the cumulative incidence code*/

/*poisson*/

%let mydir=P:\ACCES\ACCES_Research\Qingyun\previous abnormality\Doc\NKCx_260604_exnsh\;

proc export data=pop_sur_time outfile="&mydir.\cox_stata.dta" replace dbms=stata replace; run; 
data pop_sur_time;
set pop_sur_time;
if sur_cin2=. then sur_cin2=0.0000001;
run;


data pop_sur_time;
set pop_sur_time;
log_cin2=log(sur_cin2);
log_cancer=log(sur_cancer);
if noscr_his=1 then do; previous_pos=.;previous_pos_cat=.;end;
if previous_pos=1 then  noscr_his=.;
rename fu_cin2=cin2;
run;





proc genmod data=pop_sur_time;
class previous_pos_cat(ref="0") birthco hpv_year/param=glm;
model cin2=previous_pos_cat  / type3 dist=poisson link=log offset=log_cin2;
store p1;
ods output ParameterEstimates=poi_pos_cin2_cat;
where previous_pos_cat in (0 1) and HPVrisk=0;
run;
data poi_pos_cin2_cat;
set poi_pos_cin2_cat;
format char $30.;
keep RR L_RR H_RR poi_p;
char=cat(compress(parameter),'-',compress(level1));
RR=exp(estimate);
L_RR=exp(lowerwaldcl);
H_RR=exp(upperwaldcl);
poi_p=probchisq;
if parameter='previous_pos_cat';
run;

option mprint mlogic;
%macro poisson(outcome, exposure,name,subgroup);
proc genmod data=pop_sur_time;
class &exposure.(ref="0") birthco hpv_year/param=glm;
model &outcome.=&exposure. / type3 dist=poisson link=log offset=log_&outcome.;
%if &subgroup.=POS %then %do;
where HPVrisk>0;
%end;
%else %if &subgroup.=NEG %then %do;
where HPVrisk=0;
%end;
store p1;
ods output ParameterEstimates=poiorg_&outcome._&name._&subgroup.;
run;
data poiorg_&outcome._&name._&subgroup.;
set poiorg_&outcome._&name._&subgroup.;
format char $30.;
keep char RR L_RR H_RR poi_p;
char=cat("&subgroup.",'-',compress(parameter),'-',compress(level1));
RR=exp(estimate);
L_RR=exp(lowerwaldcl);
H_RR=exp(upperwaldcl);
poi_p=probchisq;
if parameter="&exposure.";
run;
proc genmod data=pop_sur_time;
class &exposure.(ref="0") birthco hpv_year/param=glm;
model &outcome.=&exposure. birthco hpv_year/ type3 dist=poisson link=log offset=log_&outcome.;
%if &subgroup.=POS %then %do;
where HPVrisk>0;
%end;
%else %if &subgroup.=NEG %then %do;
where HPVrisk=0;
%end;
store p1;
ods output ParameterEstimates=poiadj_&outcome._&name._&subgroup.;
run;
data poiadj_&outcome._&name._&subgroup.;
set poiadj_&outcome._&name._&subgroup.;
format char $30.;
keep char RR L_RR H_RR poi_p;
char=cat("&subgroup.",'-',compress(parameter),'-',compress(level1));
RR=exp(estimate);
L_RR=exp(lowerwaldcl);
H_RR=exp(upperwaldcl);
poi_p=probchisq;
if parameter="&exposure.";
run;
%mend;
%poisson(cin2,previous_pos,pos,All);
%poisson(cin2,previous_pos,pos,POS);
%poisson(cin2,previous_pos,pos,NEG);
%poisson(cin2,previous_pos_cat,cat,All);
%poisson(cin2,previous_pos_cat,cat,POS);
%poisson(cin2,previous_pos_cat,cat,NEG);
%poisson(cin2,noscr_his,scr,All);
%poisson(cin2,noscr_his,scr,POS);
%poisson(cin2,noscr_his,sct,NEG);



/*organize*/

data poi_cin2_org;
set poiorg_cin2_:;
if index(char,'0')>0 then delete;
run; 
data poi_cin2_org;
set poi_cin2_org;
keep char CIN2_org_RR CIN2_org_p;
CIN2_org_p=poi_p;
CIN2_org_RR=cat(compress(put(RR,8.2)),' (',compress(put(L_RR,8.2)),', ',compress(put(H_RR,8.2)),')');
run;
data poi_cin2_adj;
set poiadj_cin2_:;
if index(char,'0')>0 then delete;
run; 
data poi_cin2_adj;
set poi_cin2_adj;
keep char CIN2_adj_RR CIN2_adj_p;
CIN2_adj_p=poi_p;
CIN2_adj_RR=cat(compress(put(RR,8.2)),' (',compress(put(L_RR,8.2)),', ',compress(put(H_RR,8.2)),')');
run;
/*cancer*/

%poisson(cancer,previous_pos,pos,All);
%poisson(cancer,previous_pos,pos,POS);
%poisson(cancer,previous_pos,pos,NEG);
%poisson(cancer,previous_pos_cat,cat,All);
%poisson(cancer,previous_pos_cat,cat,POS);
%poisson(cancer,previous_pos_cat,cat,NEG);
%poisson(cancer,noscr_his,scr,All);
%poisson(cancer,noscr_his,scr,POS);
%poisson(cancer,noscr_his,sct,NEG);
/*organize*/

data poi_cancer_org;
set poiorg_cancer_:;
if index(char,'0')>0 then delete;
run; 
data poi_cancer_org;
set poi_cancer_org;
keep char ICC_org_RR ICC_org_p;
ICC_org_p=poi_p;
ICC_org_RR=cat(compress(put(RR,8.2)),' (',compress(put(L_RR,8.2)),', ',compress(put(H_RR,8.2)),')');
run;
data poi_cancer_adj;
set poiadj_cancer_:;
if index(char,'0')>0 then delete;
run; 
data poi_cancer_adj;
set poi_cancer_adj;
keep char ICC_adj_RR ICC_adj_p;
ICC_adj_p=poi_p;
ICC_adj_RR=cat(compress(put(RR,8.2)),' (',compress(put(L_RR,8.2)),', ',compress(put(H_RR,8.2)),')');
run;


/*save organize poisson results*/
proc sort data=poi_cin2_org;
by char;
proc sort data=poi_cin2_adj;
by char;
proc sort data=poi_cancer_org;
by char;
proc sort data=poi_cancer_adj;
by char;
run;

data poisson;
merge poi_cin2_org poi_cin2_adj poi_cancer_org poi_cancer_adj;
by char;
run;


title 'poisson results';
ods rtf file="&mydir.cin2 and ICC poisson_&sysdate..rtf";
proc print data=poisson noobs;run;
ods rtf close;
title;



proc datasets library=work nolist;
delete poi:;
quit;



/*cox*/


proc phreg data=pop_sur_time simple;
class previous_pos (param=ref ref="0") hpv_year birthco;
model sur_cancer*cancer(0)=previous_pos /ties=exact rl=pl type3(all);
where HPVrisk=0;
ods output ParameterEstimates=cox_pos_cancer_cat;
run;
data cox_pos_cancer_cat;
set cox_pos_cancer_cat;
keep label estimate stderr hazardratio HRLowerPLCL HRUpperPLCL probchisq ;
rename hazardratio=HR HRLowerPLCL=L_HR HRUpperPLCL=H_HR probchisq=cox_p_org;
L_RR=exp(lowerwaldcl);
H_RR=exp(upperwaldcl);
poi_p=probchisq;
if parameter='previous_pos';
run;

%macro cox(outcome, exposure,name,subgroup);
proc phreg data=pop_sur_time simple;
class &exposure. (param=ref ref="0") hpv_year birthco;
model sur_&outcome.*&outcome.(0)=&exposure. /ties=exact rl=pl type3(all);
%if &subgroup.=POS %then %do;
where HPVrisk>0;
%end;
%else %if &subgroup.=NEG %then %do;
where HPVrisk=0;
%end;
ods output ParameterEstimates=coxorg_&outcome._&name._&subgroup.;
run;
data coxorg_&outcome._&name._&subgroup.;
set coxorg_&outcome._&name._&subgroup.;
format char $30.;
char=cat("&subgroup.",'-',compress(label));
keep char estimate stderr hazardratio HRLowerPLCL HRUpperPLCL probchisq ;
rename hazardratio=HR HRLowerPLCL=L_HR HRUpperPLCL=H_HR probchisq=cox_p_org;
if parameter="&exposure.";
run;

proc phreg data=pop_sur_time simple;
class &exposure. (param=ref ref="0") hpv_year birthco;
model sur_&outcome.*&outcome.(0)=&exposure. hpv_year birthco/ties=exact rl=pl type3(all);
%if &subgroup.=POS %then %do;
where HPVrisk>0;
%end;
%else %if &subgroup.=NEG %then %do;
where HPVrisk=0;
%end;
ods output ParameterEstimates=coxadj_&outcome._&name._&subgroup.;
run;
data coxadj_&outcome._&name._&subgroup.;
set coxadj_&outcome._&name._&subgroup.;
format char $30.;
char=cat("&subgroup.",'-',compress(label));
keep char estimate stderr hazardratio HRLowerPLCL HRUpperPLCL probchisq ;
rename hazardratio=HR HRLowerPLCL=L_HR HRUpperPLCL=H_HR probchisq=cox_p_adj;
if parameter="&exposure.";
run;
%mend;

%cox(cancer,previous_pos,pos,All);
%cox(cancer,previous_pos,pos,POS);
%cox(cancer,previous_pos,pos,NEG);
%cox(cancer,previous_pos_cat,cat,All);
%cox(cancer,previous_pos_cat,cat,POS);
%cox(cancer,previous_pos_cat,cat,NEG);
%cox(cancer,noscr_his,scr,All);
%cox(cancer,noscr_his,scr,POS);
%cox(cancer,noscr_his,sct,NEG);

%cox(cin2,previous_pos,pos,All);
%cox(cin2,previous_pos,pos,POS);
%cox(cin2,previous_pos,pos,NEG);
%cox(cin2,previous_pos_cat,cat,All);
%cox(cin2,previous_pos_cat,cat,POS);
%cox(cin2,previous_pos_cat,cat,NEG);
%cox(cin2,noscr_his,scr,All);
%cox(cin2,noscr_his,scr,POS);
%cox(cin2,noscr_his,sct,NEG);


/*floting zero point so group hpv_year to avoid this problem*/
data pop_sur_time;
set pop_sur_time;
if hpv_year>=2020 then hpv_period=1;
else if hpv_year>=2017 then hpv_period=2;
else hpv_period=3;
run;
proc phreg data=pop_sur_time simple;
class previous_pos (param=ref ref="0") hpv_period birthco;
model sur_cin2*cin2(0)=previous_pos hpv_period birthco/ties=exact rl=pl type3(all);
where HPVrisk=0;
ods output ParameterEstimates=coxadj_cin2_pos_neg;
run;
data coxadj_cin2_pos_neg;
set coxadj_cin2_pos_neg;
format char $30.;
char=cat("NEG",'-',compress(label));
keep char estimate stderr hazardratio HRLowerPLCL HRUpperPLCL probchisq ;
rename hazardratio=HR HRLowerPLCL=L_HR HRUpperPLCL=H_HR probchisq=cox_p_adj;
if parameter="previous_pos";
run;
proc phreg data=pop_sur_time simple;
class previous_pos_cat (param=ref ref="0") hpv_period birthco;
model sur_cin2*cin2(0)=previous_pos_cat hpv_period birthco/ties=exact rl=pl type3(all);
where HPVrisk=0;
ods output ParameterEstimates=coxadj_cin2_cat_neg;
run;
data coxadj_cin2_cat_neg;
set coxadj_cin2_cat_neg;
format char $30.;
char=cat("NEG",'-',compress(label));
keep char estimate stderr hazardratio HRLowerPLCL HRUpperPLCL probchisq ;
rename hazardratio=HR HRLowerPLCL=L_HR HRUpperPLCL=H_HR probchisq=cox_p_adj;
if parameter="previous_pos_cat";
run;


/*organize*/

data cox_cin2_org;
set coxorg_cin2_:;
run; 
data cox_cin2_org;
set cox_cin2_org;
keep char CIN2_org_HR CIN2_org_p;
CIN2_org_p=cox_p_org;
CIN2_org_HR=cat(compress(put(HR,8.2)),' (',compress(put(L_HR,8.2)),', ',compress(put(H_HR,8.2)),')');
run;
data cox_cin2_adj;
set coxadj_cin2_:;
run; 
data cox_cin2_adj;
set cox_cin2_adj;
keep char CIN2_adj_HR CIN2_adj_p;
CIN2_adj_p=cox_p_adj;
CIN2_adj_HR=cat(compress(put(HR,8.2)),' (',compress(put(L_HR,8.2)),', ',compress(put(H_HR,8.2)),')');
run;

data cox_cancer_org;
set coxorg_cancer_:;
run; 
data cox_cancer_org;
set cox_cancer_org;
keep char ICC_org_HR ICC_org_p;
ICC_org_p=cox_p_org;
ICC_org_HR=cat(compress(put(HR,8.2)),' (',compress(put(L_HR,8.2)),', ',compress(put(H_HR,8.2)),')');
run;
data cox_cancer_adj;
set coxadj_cancer_:;
run; 
data cox_cancer_adj;
set cox_cancer_adj;
keep char ICC_adj_HR ICC_adj_p;
ICC_adj_p=cox_p_adj;
ICC_adj_HR=cat(compress(put(HR,8.2)),' (',compress(put(L_HR,8.2)),', ',compress(put(H_HR,8.2)),')');
run;


/*run both poisson and cox but keep cox as the main results*/

proc phreg data=pop_sur_time simple;
class previous_pos (param=ref ref="0") hpv_year birthco;
model sur_cancer*cancer(0)=previous_pos hpv_year birthco/ties=exact rl=pl type3(all);
*where HPVrisk=0;
ods output ParameterEstimates=test;
run;




/*save organize poisson results*/
proc sort data=cox_cin2_org;
by char;
proc sort data=cox_cin2_adj;
by char;
proc sort data=cox_cancer_org;
by char;
proc sort data=cox_cancer_adj;
by char;
run;

data cox;
merge cox_cin2_org cox_cin2_adj cox_cancer_org cox_cancer_adj;
by char;
run;


title 'Cox results';
ods rtf file="&mydir.cin2 and ICC cox&sysdate..rtf";
proc print data=cox noobs;run;
ods rtf close;
title;
