%% objective function for fitting confidence model
function [sse, confpred, r2] = func_conf_HIST_TE(params,df,idx,sub, pr)
if nargin < 5
    pr = 0; %param recovery
end
% params
OFFSET = params(1);
TE_WEIGHT = params(2);
K = params(3);

ntrials = sum(idx==1);
rotation = df.rot(idx);
ha = df.ha(sub,idx); % hand angle data
ha=ha';

te = abs(ha+rotation);

if pr == 1 %param recovery mode
    conf = df.confpred_HIST_TE(sub, :);
else
    conf = df.conf(sub,idx); % confidence ratings
end

confpred = nan(1,length(ha));
est_err = 0;

for t = 2:ntrials
    del_err = te(t-1)-est_err;
    est_err = est_err + K*del_err;
    confpred(t) = OFFSET - TE_WEIGHT * est_err;
end

valid = find(~isnan(conf));
sse = nansum((conf(valid)-confpred(valid)).^2);
re = nansum((conf(valid)-nanmean(conf(valid))).^2); 
r2 = 1-sse/re;