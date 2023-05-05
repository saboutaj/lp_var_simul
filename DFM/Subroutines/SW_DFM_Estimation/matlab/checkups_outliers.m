%%saves into csv files the datasets with and without outlier adjustment
clear all; 

load_data = 1;
levels = 0;
i_demean = 0;
bw_bw = 100; 

[dnobs_m,calvec_m,calds_m] = calendar_make([1959 1],[2014 12],12);  % Monthly Calendar
[dnobs_q,calvec_q,calds_q] = calendar_make([1959 1],[2014 4],4);    % Quarterly Calendar

datain_all;

%saving the series of: corrected data, uncorrected for outliers data,
%variable labels, and dates
csvwrite("bpdata.csv",bpdata);
csvwrite("bpdata_noa.csv",bpdata_noa);
writecell(transpose(bplabvec_long),"bplabvec_long.csv");
csvwrite("calds.csv", calds);