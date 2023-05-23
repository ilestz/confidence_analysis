%% objective function for fitting confidence model
function [sse, confpred, r2] = func_conf_1TB_TE(params,df,idx,sub, pr)
if nargin < 5
    pr = 0;
end
% params
OFFSET = params(1);
TE_WEIGHT = params(2);

ntrials = sum(idx==1);
rotation = df.rot(sub,idx);
rotation = rotation';
ha = df.ha(sub,idx); % hand angle data
ha=ha';

te = abs(ha+rotation);

if pr == 1 %param recovery mode
    conf = df.confpred_1TB_TE(sub, :);
else
    conf = df.conf(sub,idx); % confidence ratings
end

confpred = nan(1,length(ha));

for t = 2:ntrials
    confpred(t) = OFFSET - TE_WEIGHT * te(t-1);
end

valid = find(~isnan(conf));
sse = nansum((conf(valid)-confpred(valid)).^2);
re = nansum((conf(valid)-nanmean(conf(valid))).^2); 
r2 = 1-sse/re;