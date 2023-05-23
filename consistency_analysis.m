clc;
exp_num = 1;
if exp_num == 1
    load dataConf_EXP1.mat; % load data frame; EXP1 or EXP2
    dataConf = dataConf_EXP1; %EXP1 or EXP2
elseif exp_num == 2
    load dataConf_EXP2.mat; % load data frame; EXP1 or EXP2
    dataConf = dataConf_EXP2; %EXP1 or EXP2
elseif exp_num == 3
    load dataConf_online_gradual.mat; % load data frame; EXP1 or EXP2
    dataConf = dataConf_online_gradual; %EXP1 or EXP2
else
    load dataConf_online_zeromean.mat; % load data frame; EXP1 or EXP2
    dataConf = dataConf_online_zeromean; %EXP1 or EXP2
end
[nsubs,N] = size(dataConf.ha); % get dimensions

%ha and aim are always positive if in the opposite direction of rot; to
%compute the errors, we need to subtract the absolute value of rotation
if exp_num<4
    dataConf.rot = -1 * abs(dataConf.rot); %negative abs of rot; err = ha + rot
end
%additionally, confidence is flipped, with 100 as low confidence and 0 as
%high confidence; we should fix this here
if exp_num <3
    dataConf.conf = 100-dataConf.conf;
end
sub_rs = nan(nsubs, 1);
sub_ps = nan(nsubs,1);
sub_wash = nan(nsubs, 1);
sub_ret = nan(nsubs,1);
for s = 1:nsubs
    if exp_num < 3
        %aiming consistency and tracking w/ confidence
        aim_hist = nan(10, 1);
        rot_phase = dataConf.conf~=100; % fitting learning phase only (for now)
        x = 1:360;
        rot_phase= x(rot_phase(1,:));
        sub_aimsd = nan(length(rot_phase), 1);
        sub_conf = dataConf.conf(s, rot_phase);
        confidence(s,:) = sub_conf;
        te(s,:) = abs(dataConf.ha(s, :) + dataConf.rot(s,:));
        ha(s,:) = dataConf.ha(s,:);
        for t = 1:length(rot_phase)
            tt = rot_phase(t);
            aim_hist(2:10) = aim_hist(1:9);
            aim_hist(1) = dataConf.aim(s, tt);
            sub_aimsd(t) = nanstd(aim_hist);
        end
        [R, P] = corrcoef(sub_conf, sub_aimsd, 'rows', 'complete');
        sub_rs(s) = R(1,2);
        sub_ps(s) = P(1,2);
        if exp_num == 1
            sub_wash(s) = nanmean(ha(s,313:360));
            sub_ret(s) = nanmean(ha(s, 313:317))/nanmean(ha(s, 308:312));
        end
        
    end
end

if exp_num == 1

    mean(sub_wash)
    std(sub_wash)
    [h,p,~,stats] = ttest(sub_wash) %is it significantly different from zero

    mean_conf = nanmean(confidence);
    mean_drop = nanmean(nanmean(confidence(:, 52:56)')-nanmean(confidence(:, 44:48)'))
    std_drop = nanstd(nanmean(confidence(:, 52:56)')-nanmean(confidence(:, 44:48)'))
    mean_base_conf = mean(nanmean(confidence(:, 1:48)'))
    std_base_conf = std(nanmean(confidence(:, 1:48)'))

    mean_asym = mean(nanmean(ha(:, 293:312)'))
    std_asym = std(nanmean(ha(:, 293:312)'))

    sub_asym_conf = nanmean(confidence(:, 269:288)');
    sub_base_conf = nanmean(confidence(:, 1:48)')
    [~,p,~,stats] = ttest(sub_base_conf, sub_asym_conf)
elseif exp_num == 2
    for s = 1:nsubs
        rot = dataConf.rot(s, 33:end)';
        trial_pos = 1;
        block = 1;
        for t = 2:length(rot)
            if rot(t) == rot(t-1)
                trial_pos = [trial_pos trial_pos(t-1)+1];
                block = [block block(t-1)];
            else
                trial_pos = [trial_pos 1];
                block = [block block(t-1)+1];
            end
        end
        mean_ha = te(s, 33:end)';
        comp = 1-(-1*mean_ha./rot);
        comp(isinf(comp)) = NaN;
        conf = confidence(s, 17:end)';
        block = block';
        trial_pos = trial_pos';
        % mean_ha = mean_ha(abs(rot)>0);
        % mean_conf = mean_conf(abs(rot)>0);
        trial_pos = trial_pos(abs(rot)>0);
        block = block(abs(rot)>0);
        rot = rot(abs(rot)>0);
        X = [trial_pos block trial_pos.*block ones(length(block),1)];
        Y = comp;
        [B_ha(s,:), ~,~,~, stats] = regress(Y, X);
        ha_rs(s) = stats(1);
        Y = conf;
        [B_conf(s,:), ~,~,~, stats] = regress(Y, X);
        conf_rs(s) = stats(1);
        % conf_lm = fitlm([rot, trial_pos, block], mean_conf, [0 0 0 0; 1 0 0 0; 0 1 0 0; 0 0 1 0; 0 1 1 0]);
    end
    mean_ha = nanmean(te(:, 33:end))';
    comp = 1-(-1*mean_ha./rot);


elseif exp_num == 3
    good_subs3 = [1:17 19:25];
    dataConf.ha = dataConf.ha(good_subs3, :);
    dataConf.conf = dataConf.conf(good_subs3, :);
    dataConf.rot = dataConf.rot(good_subs3, :);

    dataConf.te = abs(dataConf.ha - abs(dataConf.rot));
    %asymptote
    sub_asymptote = nanmean(dataConf.ha(:, 161:180), 2);
    mean_asym = mean(sub_asymptote)
    sd_asym = std(sub_asymptote)
    % %compensation
    mean_comp = mean(sub_asymptote)/15
    sd_comp = std(sub_asymptote)/15

    %washout
    sub_wash = nanmean(dataConf.ha(:, 181:200), 2);
    mean_wash = mean(sub_wash)
    sd_wash = std(sub_wash)
    [~,p,~,stats] = ttest(sub_asymptote)

    %retention ratio
    sub_wo = nanmean(dataConf.ha(:, 181:185), 2);
    sub_asym = nanmean(dataConf.ha(:, 176:180), 2);
    sub_rr = sub_wo./sub_asym;
    mean_ret = mean(sub_rr)
    sd_ret = std(sub_rr)
    [~,p,~,stats] = ttest(sub_rr)

    %confidence analyses
    base_conf = nanmean(dataConf.conf(:, 11:30), 2);
    late_conf = nanmean(dataConf.conf(:, 161:180), 2);

    mean_base_conf = mean(base_conf)
    std_base_conf = std(base_conf)

    mean_late_conf = mean(late_conf)
    std_late_conf = std(late_conf)

    [~,p,~,stats] = ttest(base_conf, late_conf)

    %similar stuff for te?
    base_te = nanmean(dataConf.te(:, 11:30), 2);
    late_te = nanmean(dataConf.te(:, 161:180), 2);

    mean_base_te = mean(base_te)
    std_base_te = std(base_te)

    mean_late_te = mean(late_te)
    std_late_te = std(late_te)

    [~,p,~,stats] = ttest(base_te, late_te)

elseif exp_num == 4
    for s = 1:nsubs
        sub_ha = -dataConf.ha(s, 31:180).*sign(dataConf.rot(s, 31:180)); %raw ha
        sub_delta = diff(sub_ha);
        sub_rot = dataConf.rot(s, 31:180);
        sub_te = sub_ha + sub_rot;
        dataConf.te(s,:) = abs(sub_te);
        sub_stl_trace = -sub_delta.*sign(sub_te(1:(end-1))); %compensatory direction for previous rotation
        sub_stl(s) = nanmean(sub_stl_trace);
        early_stl(s) = nanmean(sub_stl_trace(1:20));
        late_stl(s) = nanmean(sub_stl_trace((end-20+1):end));
    end
    %sub_stl = nanmean(-diff(dataConf.ha(:,31:180)')'.*sign(dataConf.rot(:, 31:179)));
    good_subs4 = [1:13 15:23 25:29];
    sub_stl = sub_stl(good_subs4);
    mean_stl = mean(sub_stl)
    sd_stl = std(sub_stl)
    [~,p,~,stats] = ttest(sub_stl)

    %compare across blocks
    early_stl = early_stl(good_subs4);
    late_stl = late_stl(good_subs4);
    mean_early_stl = mean(early_stl)
    mean_late_stl = mean(late_stl)
    [~,p,~,stats] = ttest(early_stl, late_stl)

    %confidence analyses
    base_conf = nanmean(dataConf.conf(good_subs4, 11:30), 2);
    late_conf = nanmean(dataConf.conf(good_subs4, 161:180), 2);

    mean_base_conf = mean(base_conf)
    std_base_conf = std(base_conf)

    mean_late_conf = mean(late_conf)
    std_late_conf = std(late_conf)

    [~,p,~,stats] = ttest(base_conf, late_conf)

    %TE over time?
    base_te = nanmean(dataConf.te(good_subs4, 1:20), 2);
    late_te = nanmean(dataConf.te(good_subs4, 131:150), 2);

    mean_base_te = mean(base_te)
    std_base_te = std(base_te)

    mean_late_te = mean(late_te)
    std_late_te = std(late_te)

    [~,p,~,stats] = ttest(base_te, late_te)
end


function shadedErrorPlot(data, color, alpha)
    std_dat = nanstd(data);
    sem_dat = std_dat/sqrt(size(data,1));
    mean_dat = nanmean(data);
    x_vec = [(1:length(data)) fliplr(1:length(data))];
    patch = fill(x_vec, [mean_dat-sem_dat,fliplr(mean_dat+sem_dat)], color);
    set(patch, 'edgecolor', 'none');
    set(patch, 'FaceAlpha', alpha);
end

