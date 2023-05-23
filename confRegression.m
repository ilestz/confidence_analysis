clearvars -except confReg; clc;
exp_num = 2;
if exp_num == 1
    load dataConf_EXP1.mat; % load data frame; EXP1 or EXP2
    dataConf = dataConf_EXP1;
    rot_phase = 73:312;
    good_subs = 1:18;
elseif exp_num == 2
    load dataConf_EXP2.mat; % load data frame; EXP1 or EXP2
    dataConf = dataConf_EXP2;
    rot_phase = 33:240;
    good_subs = 1:20;
elseif exp_num == 3
    load dataConf_online_gradual.mat; % load data frame; EXP1 or EXP2
    dataConf = dataConf_online_gradual;
    good_subs = [1:17 19:25];
    rot_phase = 31:180;
else
    load dataConf_online_zeromean.mat; % load data frame; EXP1 or EXP2
    dataConf = dataConf_online_zeromean;
    good_subs = [1:13 15:23 25:29];
    rot_phase = 31:180;
end

[nsubs,N] = size(dataConf.ha); % get dimensions
if exp_num == 3
    nsubs = nsubs-1;
elseif exp_num == 4
    nsubs = nsubs-2;
end

%ha and aim are always positive if in the opposite direction of rot; to
%compute the errors, we need to subtract the absolute value of rotation
dataConf.raw_ha = -1*sign(dataConf.rot+eps).*dataConf.ha; %restore raw hand angle
if exp_num < 3
    dataConf.raw_aim = -1*sign(dataConf.rot+ eps).*dataConf.aim;
end
dataConf.raw_rot = dataConf.rot;
dataConf.rot = -1 * abs(dataConf.rot); %negative abs of rot; err = ha + rot
%additionally, confidence is flipped, with 100 as low confidence and 0 as
%high confidence; we should fix this here
if exp_num < 3 %only for experiments 1 and 2; online was correct
    dataConf.conf = 100-dataConf.conf;
    dataConf.conf(isnan(dataConf.conf)) = 50;
end

for si = 1:nsubs
    s = good_subs(si);
    %isolate only the rotation phase
    sub_rot_ha = dataConf.raw_ha(s, rot_phase);
    % figure;
    % plot(sub_rot_ha);
    % pause(0.3)
    if exp_num < 3
        sub_rot_aim = dataConf.raw_aim(s, rot_phase);
    end
    sub_rot = dataConf.rot(s, rot_phase);
    sub_raw_rot = dataConf.raw_rot(s, rot_phase);
    %cursor error
    %te = abs(dataConf.ha(s, rot_phase)+ dataConf.rot(s, rot_phase));
    %te = te(1:(end-1));
    raw_te = sub_rot_ha+ sub_raw_rot;
    raw_te = raw_te(1:(end-1));
    te = abs(raw_te); %magnitude of the error
    %delta, in the compensatory direction
    delta_ha = -diff(sub_rot_ha).*sign(eps+raw_te);
    %delta_ha = diff(sub_rot_ha);
    if exp_num < 3
        delta_aim = -diff(sub_rot_aim).*sign(eps+raw_te);
        %delta_aim = diff(sub_rot_aim);
    end
    %confidence
    conf = dataConf.conf(s, rot_phase);
    conf = conf(1:(end-1));
    if exp_num < 3
        conf(conf == 50) = NaN;
    end
    %conf(conf >98) = NaN;
    rmv(s) = sum(isnan(conf));

    % if(s < 5)
    %     figure; hold on;
    %     plot(sub_rot_ha);
    %     plot(-sub_raw_rot);
    %     plot(raw_te);
    % 
    %     figure; hold on;
    %     plot(te);
    %     plot(delta_ha);
    % end

    %trial number (within a miniblock)
    trial_num = nan(length(rot_phase)-1,1);
    trial_num(1) = 1;
    for i = 2:length(trial_num)
        if sub_rot(i)==sub_rot(i-1) || exp_num > 2
            trial_num(i) = trial_num(i-1)+1;
        else
            trial_num(i) = 1;
        end
    end

    %here I would z-score but I don't understand why we would z-score
    %I am going to keep everything in raw degrees
    %update: Z-scoring doesn't change the results
    %cleaning step + transposition
    delta_ha(delta_ha<-20) = NaN; %wrong way, by mistake
    delta_ha = delta_ha';
    delta_ha(nan_zscore(delta_ha)>3) = NaN;
    %delta_ha = nan_zscore(delta_ha);
    if exp_num < 3
        delta_aim(delta_aim<-20) = NaN; %wrong way, by mistake
        delta_aim = delta_aim';
        delta_aim(nan_zscore(delta_aim)>3) = NaN;
        %delta_aim = nan_zscore(delta_aim);
    end
    te = te';
    raw_te(nan_zscore(te)>3) = NaN;
    te(nan_zscore(te)>3) = NaN;
    %te = nan_zscore(te);
    conf = nan_zscore(conf)'; %scaling on confidence is arbitrary

    %te = raw_te';
    %finally, the regression
    %aim + confidence (ac)
    if exp_num < 3
        X = [trial_num te conf te.*conf ones(length(conf),1)];
        [B_ac(s, :), ~, R,~, stats] = regress(delta_aim, X);
        sse = nansum(R.^2);
        r2_ac(s) = stats(1);
        k = length(X(1,:)) - 1;
        aic_ac(si) = k + length(X) * (log((2*pi)*sse./(length(X)-k)) + 1);

        %aim only (a)
        X = [trial_num te ones(length(conf),1)];
        [B_a(s,:), ~, R,~, stats] = regress(delta_aim, X);
        sse = nansum(R.^2);
        r2_a(s) = stats(1);
        k = length(X(1,:)) - 1;
        aic_a(si) = k + length(X) * (log((2*pi)*sse./(length(X)-k)) + 1);
    end

    %ha + confidence (hc)
    X = [trial_num te conf te.*conf ones(length(conf),1)];
    [B_hc(s,:), ~, R,~, stats] = regress(delta_ha, X);
    sse = nansum(R.^2);
    r2_hc(s) = stats(1);
    k = length(X(1,:)) - 1;
    aic_hc(si) = k + length(X) * (log((2*pi)*sse./(length(X)-k)) + 1);

    %ha only (h)
    X = [trial_num te ones(length(conf),1)];
    [B_h(s,:), ~, R,~, stats] = regress(delta_ha, X);
    sse = nansum(R.^2);
    r2_h(s) = stats(1);
    k = length(X(1,:)) - 1;
    aic_h(si) = k + length(X) * (log((2*pi)*sse./(length(X)-k)) + 1);
end

%results
if exp_num <3
conf_param_aim = mean(B_ac(:, 3))
[~,p_n] = kstest(B_ac(:, 3))
[p,~,stats] = signrank(B_ac(:, 3))
aic_diff = mean(aic_ac-aic_a)
[p_m, ~, stats] = signrank(aic_ac - aic_a)
figure;
bar(1:length(aic_a), sort(aic_ac - aic_a))
end

conf_param_ha = mean(B_hc(:, 3))
[~,p_n] = kstest(B_hc(:, 3))
[p, ~, stats] = signrank(B_hc(:, 3))
aic_diff = mean(aic_hc-aic_h)
[p_m, ~, stats] = signrank(aic_hc- aic_h)
figure;
bar(1:length(aic_hc), sort(aic_hc-aic_h))

%confReg.aicdiff_ha_exp4 = aic_hc - aic_h;
