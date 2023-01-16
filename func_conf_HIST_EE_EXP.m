%% objective function for fitting confidence model
function [sse, confpred, r2] = func_conf_HIST_EE_EXP(params,df,idx,sub)
% params
OFFSET = params(1);
EE_WEIGHT = params(2);
EE_EXP = params(3);
K = params(4);

ntrials = sum(idx==1);
rotation = df.rot(idx);
aim = df.aim(sub,idx); % hand angle data
aim=aim';

ee = abs(aim+rotation).^EE_EXP;
conf = df.conf(sub,idx); % confidence ratings

confpred = nan(1,length(aim));
est_err = 0;

for t = 2:ntrials
    del_err = ee(t-1)-est_err;
    est_err = est_err + K*del_err;
    confpred(t) = OFFSET - EE_WEIGHT * est_err;
end

valid = find(~isnan(conf));
sse = nansum((conf(valid)-confpred(valid)).^2);
re = nansum((conf(valid)-nanmean(conf(valid))).^2); 
r2 = 1-sse/re;