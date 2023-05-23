clear; clc;
exp_num = 2;
if exp_num == 1
    load dataConf_EXP1.mat; % load data frame; EXP1 or EXP2
    dataConf = dataConf_EXP1; %EXP1 or EXP2
else
    load dataConf_EXP2.mat; % load data frame; EXP1 or EXP2
    dataConf = dataConf_EXP2; %EXP1 or EXP2
end
[nsubs,N] = size(dataConf.ha); % get dimensions

%ha and aim are always positive if in the opposite direction of rot; to
%compute the errors, we need to subtract the absolute value of rotation
dataConf.rot = -1 * abs(dataConf.rot); %negative abs of rot; err = ha + rot
%additionally, confidence is flipped, with 100 as low confidence and 0 as
%high confidence; we should fix this here
dataConf.conf = 100-dataConf.conf;

aim_conf_corr = nan(2, nsubs);
ha_conf_corr = nan(2, nsubs);
imp_conf_corr = nan(2, nsubs);
aim_sd = nan(1, nsubs);
conf_sd = nan(1,nsubs);
for si = 1:nsubs
    rot_phase = dataConf.conf~=100; % fitting learning phase only (for now)
    rot_phase= rot_phase(1,:);
    sub_ha = dataConf.ha(si, rot_phase);
    sub_aim = dataConf.aim(si, rot_phase);
    sub_imp = dataConf.implicit(si, rot_phase);
    sub_conf = dataConf.conf(si, rot_phase);
    
    [r,p] = corrcoef(diff(sub_aim), sub_conf(1:(end-1)), 'rows', 'complete');
    aim_conf_corr(1,si) = r(1,2);
    aim_conf_corr(2,si) = p(1,2);
    
    [r,p] = corrcoef(diff(sub_ha), sub_conf(1:(end-1)), 'rows', 'complete');
    ha_conf_corr(1,si) = r(1,2);
    ha_conf_corr(2,si) = p(1,2);
    
    [r,p] = corrcoef(diff(sub_imp), sub_conf(1:(end-1)), 'rows', 'complete');
    imp_conf_corr(1,si) = r(1,2);
    imp_conf_corr(2,si) = p(1,2);
    
    aim_sd(si) = nanstd(sub_aim);
    conf_sd(si) = nanstd(sub_conf);
    ha_sd(si) = nanstd(sub_ha);
end