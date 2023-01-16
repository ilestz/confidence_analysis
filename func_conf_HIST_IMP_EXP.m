%% objective function for fitting confidence model
function [sse, confpred, r2] = func_conf_HIST_IMP_EXP(params,df,idx,sub)
% params
OFFSET = params(1);
IMP_WEIGHT = params(2);
IMP_EXP = params(3);
K = params(4);

ntrials = sum(idx==1);
rotation = df.rot(idx);
aim = df.aim(sub,idx); % aim angle data
aim=aim';

ha = df.ha(sub,idx); % hand angle data
ha=ha';

imp = abs(ha-aim).^IMP_EXP;
conf = df.conf(sub,idx); % confidence ratings

confpred = nan(1,length(ha));
est_err = 0;

for t = 2:ntrials
    del_err = imp(t-1)-est_err;
    est_err = est_err + K*del_err;
    confpred(t) = OFFSET - IMP_WEIGHT * est_err;
end

valid = find(~isnan(conf));
sse = nansum((conf(valid)-confpred(valid)).^2);
re = nansum((conf(valid)-nanmean(conf(valid))).^2); 
r2 = 1-sse/re;