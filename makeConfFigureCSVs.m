%so first things first, lets do hand angle curves
load dataConf_EXP1.mat dataConf_EXP1
dataConf1 = dataConf_EXP1;
load dataConf_EXP2.mat dataConf_EXP2
dataConf2 = dataConf_EXP2;
load dataConf_online_gradual.mat dataConf_online_gradual
dataConf3 = dataConf_online_gradual;
good_subs3 = [1:17 19:25];
load dataConf_online_zeromean.mat dataConf_online_zeromean
dataConf4 = dataConf_online_zeromean;
good_subs4 = [1:13 15:23 25:29];

numsubs1 = size(dataConf1.ha, 1);
numsubs2 = size(dataConf2.ha, 1);
numsubs3 = length(good_subs3);
numsubs4 = length(good_subs4);

conf1.ha = nanmean(dataConf1.ha)';
conf1.ha_sem = nanstd(dataConf1.ha)'/sqrt(numsubs1);
conf2.ha = nanmean(dataConf2.ha)';
conf2.ha_sem = nanstd(dataConf2.ha)'/sqrt(numsubs2);
conf3.ha = nanmean(dataConf3.ha(good_subs3, :))';
conf3.ha_sem = nanstd(dataConf3.ha(good_subs3, :))'/sqrt(numsubs3);
conf4.ha = nanmean(dataConf4.ha(good_subs4, :))';
conf4.ha_sem = nanstd(dataConf4.ha(good_subs4, :))'/sqrt(numsubs4);

%again for confidence
conf1.conf = nanmean(dataConf1.conf)';
conf1.conf_sem = nanstd(dataConf1.conf)'/sqrt(numsubs1);
conf2.conf = nanmean(dataConf2.conf)';
conf2.conf_sem = nanstd(dataConf2.conf)'/sqrt(numsubs2);
conf3.conf = nanmean(dataConf3.conf(good_subs3, :))';
conf3.conf_sem = nanstd(dataConf3.conf(good_subs3, :))'/sqrt(numsubs3);
conf4.conf = nanmean(dataConf4.conf(good_subs4, :))';
conf4.conf_sem = nanstd(dataConf4.conf(good_subs4, :))'/sqrt(numsubs4);

%finally, do it for rot
conf1.rot = dataConf1.rot(1, :)';
conf2.rot = nanmean(dataConf2.rot)';
conf3.rot = [zeros(1, 30) 0.33*(1:45) repmat(15, 1, 105) zeros(1, 20)]'; %idealized
conf4.rot = zeros(1,200)'; %idealized

%now we need to add in the confidence model predictions
load modelFit_HIST_TE_EXP.mat modelFit_HIST_TE_EXP
model1 = modelFit_HIST_TE_EXP;
load modelFit_HIST_TE_EXP2.mat modelFit_HIST_TE_EXP
model2 = modelFit_HIST_TE_EXP;
load modelFit_HIST_TE_EXP3.mat modelFit_HIST_TE_EXP
model3 = modelFit_HIST_TE_EXP;
load modelFit_HIST_TE_EXP4.mat modelFit_HIST_TE_EXP
model4 = modelFit_HIST_TE_EXP;

conf1.confpred = [nan(1, 24) nanmean(model1.confpred) nan(1,48)]';
conf1.confpred_sem = [nan(1, 24) nanstd(model1.confpred)/sqrt(numsubs1) nan(1,48)]';
conf2.confpred = [nan(1,16) nanmean(model2.confpred)]';
conf2.confpred_sem = [nan(1,16) nanstd(model2.confpred)/sqrt(numsubs2)]';
conf3.confpred = [nan(1,10) nanmean(model3.confpred(good_subs3, :)) nan(1,20)]';
conf3.confpred_sem = [nan(1,10) nanstd(model3.confpred(good_subs3, :))/sqrt(numsubs3) nan(1,20)]';
conf4.confpred = [nan(1,10) nanmean(model4.confpred(good_subs4, :)) nan(1,20)]';
conf4.confpred_sem = [nan(1,10) nanstd(model4.confpred(good_subs4, :))/sqrt(numsubs4) nan(1,20)]';

writetable(struct2table(conf1), "conf_exp1.csv");
writetable(struct2table(conf2), "conf_exp2.csv");
writetable(struct2table(conf3), "conf_exp3.csv");
writetable(struct2table(conf4), "conf_exp4.csv");

%ok so now we want to aggregate winning model params, and other necessary
%things for the figures

%we will pad everything to have the same length, with the max being 29
numsubs4 = 29;

%params
conf_params.offset_exp1 = fillLength(model1.offset, numsubs4)';
conf_params.weight_exp1 = fillLength(model1.weight, numsubs4)';
conf_params.exp_exp1 = fillLength(model1.exponent, numsubs4)';
conf_params.k_exp1 = fillLength(model1.K, numsubs4)';
conf_params.offset_exp2 = fillLength(model2.offset, numsubs4)';
conf_params.weight_exp2 = fillLength(model2.weight, numsubs4)';
conf_params.exp_exp2 = fillLength(model2.exponent, numsubs4)';
conf_params.k_exp2 = fillLength(model2.K, numsubs4)';
conf_params.offset_exp3 = fillLength(model3.offset, numsubs4)';
conf_params.weight_exp3 = fillLength(model3.weight, numsubs4)';
conf_params.exp_exp3 = fillLength(model3.exponent, numsubs4)';
conf_params.k_exp3 = fillLength(model3.K, numsubs4)';
conf_params.offset_exp4 = fillLength(model4.offset, numsubs4)';
conf_params.weight_exp4 = fillLength(model4.weight, numsubs4)';
conf_params.exp_exp4 = fillLength(model4.exponent, numsubs4)';
conf_params.k_exp4 = fillLength(model4.K, numsubs4)';

%rs and true_rs
conf_params.rs_exp1 = fillLength(model1.r2, numsubs4)';
conf_params.truers_exp1 = fillLength(model1.true_r2, numsubs4)';
conf_params.rs_exp2 = fillLength(model2.r2, numsubs4)';
conf_params.truers_exp2 = fillLength(model2.true_r2, numsubs4)';
conf_params.rs_exp3 = fillLength(model3.r2, numsubs4)';
conf_params.truers_exp3 = fillLength(model3.true_r2, numsubs4)';
conf_params.rs_exp4 = fillLength(model4.r2, numsubs4)';
conf_params.truers_exp4 = fillLength(model4.true_r2, numsubs4)';

%subject wise delta aics
conf_params.aic1TB_exp1 = fillLength(model1.aic_1TB_TE, numsubs4)';
conf_params.aic1TBe_exp1 = fillLength(model1.aic_1TB_TE_EXP, numsubs4)';
conf_params.aicHIST_exp1 = fillLength(model1.aic_HIST_TE, numsubs4)';
conf_params.aicHISTe_exp1 = fillLength(model1.aic_HIST_TE_EXP, numsubs4)';
conf_params.aic1TB_exp2 = fillLength(model2.aic_1TB_TE, numsubs4)';
conf_params.aic1TBe_exp2 = fillLength(model2.aic_1TB_TE_EXP, numsubs4)';
conf_params.aicHIST_exp2 = fillLength(model2.aic_HIST_TE, numsubs4)';
conf_params.aicHISTe_exp2 = fillLength(model2.aic_HIST_TE_EXP, numsubs4)';
conf_params.aic1TB_exp3 = fillLength(model3.aic_1TB_TE, numsubs4)';
conf_params.aic1TBe_exp3 = fillLength(model3.aic_1TB_TE_EXP, numsubs4)';
conf_params.aicHIST_exp3 = fillLength(model3.aic_HIST_TE, numsubs4)';
conf_params.aicHISTe_exp3 = fillLength(model3.aic_HIST_TE_EXP, numsubs4)';
conf_params.aic1TB_exp4 = fillLength(model4.aic_1TB_TE, numsubs4)';
conf_params.aic1TBe_exp4 = fillLength(model4.aic_1TB_TE_EXP, numsubs4)';
conf_params.aicHIST_exp4 = fillLength(model4.aic_HIST_TE, numsubs4)';
conf_params.aicHISTe_exp4 = fillLength(model4.aic_HIST_TE_EXP, numsubs4)';

%finally the winning model vs. best alternative
conf_params.aicdiff_exp1 = fillLength(model1.aic_diff, numsubs4)';
conf_params.aicdiff_exp2 = fillLength(model2.aic_diff, numsubs4)';
conf_params.aicdiff_exp3 = fillLength(model3.aic_diff, numsubs4)';
conf_params.aicdiff_exp4 = fillLength(model4.aic_diff, numsubs4)';

writetable(struct2table(conf_params), "conf_params.csv");

%% a summary table of model params and statistical comparisons
param_smry.c0_exp1 = mean(model1.offset);
param_smry.c0sd_exp1 = std(model1.offset);
param_smry.c0_exp2 = mean(model2.offset);
param_smry.c0sd_exp2 = std(model2.offset);
param_smry.c0_exp3 = mean(model3.offset(good_subs3));
param_smry.c0sd_exp3 = std(model3.offset(good_subs3));
param_smry.c0_exp4 = mean(model4.offset(good_subs4));
param_smry.c0sd_exp4 = std(model4.offset(good_subs4));

param_smry.w_exp1 = mean(model1.weight);
param_smry.wsd_exp1 = std(model1.weight);
param_smry.w_exp2 = mean(model2.weight);
param_smry.wsd_exp2 = std(model2.weight);
param_smry.w_exp3 = mean(model3.weight(good_subs3));
param_smry.wsd_exp3 = std(model3.weight(good_subs3));
param_smry.w_exp4 = mean(model4.weight(good_subs4));
param_smry.wsd_exp4 = std(model4.weight(good_subs4));

param_smry.exp_exp1 = mean(model1.exponent);
param_smry.expsd_exp1 = std(model1.exponent);
param_smry.exp_exp2 = mean(model2.exponent);
param_smry.expsd_exp2 = std(model2.exponent);
param_smry.exp_exp3 = mean(model3.exponent(good_subs3));
param_smry.expsd_exp3 = std(model3.exponent(good_subs3));
param_smry.exp_exp4 = mean(model4.exponent(good_subs4));
param_smry.expsd_exp4 = std(model4.exponent(good_subs4));

param_smry.k_exp1 = mean(model1.K);
param_smry.ksd_exp1 = std(model1.K);
param_smry.k_exp2 = mean(model2.K);
param_smry.ksd_exp2 = std(model2.K);
param_smry.k_exp3 = mean(model3.K(good_subs3));
param_smry.ksd_exp3 = std(model3.K(good_subs3));
param_smry.k_exp4 = mean(model4.K(good_subs4));
param_smry.ksd_exp4 = std(model4.K(good_subs4));

[p, ~,stat] = ranksum(model1.offset, model2.offset);
param_comp.p_exp12(1) = p;
param_comp.z_exp12(1) = stat.zval;
[p, ~,stat] = ranksum(model1.weight, model2.weight);
param_comp.p_exp12(2) = p;
param_comp.z_exp12(2) = stat.zval;
[p, ~,stat] = ranksum(model1.exponent, model2.exponent);
param_comp.p_exp12(3) = p;
param_comp.z_exp12(3) = stat.zval;
[p, ~,stat] = ranksum(model1.K, model2.K);
param_comp.p_exp12(4) = p;
param_comp.z_exp12(4) = stat.zval;

[p, ~,stat] = ranksum(model1.offset, model3.offset(good_subs3));
param_comp.p_exp13(1) = p;
param_comp.z_exp13(1) = stat.zval;
[p, ~,stat] = ranksum(model1.weight, model3.weight(good_subs3));
param_comp.p_exp13(2) = p;
param_comp.z_exp13(2) = stat.zval;
[p, ~,stat] = ranksum(model1.exponent, model3.exponent(good_subs3));
param_comp.p_exp13(3) = p;
param_comp.z_exp13(3) = stat.zval;
[p, ~,stat] = ranksum(model1.K, model3.K(good_subs3));
param_comp.p_exp13(4) = p;
param_comp.z_exp13(4) = stat.zval;

[p, ~,stat] = ranksum(model1.offset, model4.offset(good_subs4));
param_comp.p_exp14(1) = p;
param_comp.z_exp14(1) = stat.zval;
[p, ~,stat] = ranksum(model1.weight, model4.weight(good_subs4));
param_comp.p_exp14(2) = p;
param_comp.z_exp14(2) = stat.zval;
[p, ~,stat] = ranksum(model1.exponent, model4.exponent(good_subs4));
param_comp.p_exp14(3) = p;
param_comp.z_exp14(3) = stat.zval;
[p, ~,stat] = ranksum(model1.K, model4.K(good_subs4));
param_comp.p_exp14(4) = p;
param_comp.z_exp14(4) = stat.zval;

[p, ~,stat] = ranksum(model2.offset, model3.offset(good_subs3));
param_comp.p_exp23(1) = p;
param_comp.z_exp23(1) = stat.zval;
[p, ~,stat] = ranksum(model2.weight, model3.weight(good_subs3));
param_comp.p_exp23(2) = p;
param_comp.z_exp23(2) = stat.zval;
[p, ~,stat] = ranksum(model2.exponent, model3.exponent(good_subs3));
param_comp.p_exp23(3) = p;
param_comp.z_exp23(3) = stat.zval;
[p, ~,stat] = ranksum(model2.K, model3.K(good_subs3));
param_comp.p_exp23(4) = p;
param_comp.z_exp23(4) = stat.zval;

[p, ~,stat] = ranksum(model2.offset, model4.offset(good_subs4));
param_comp.p_exp24(1) = p;
param_comp.z_exp24(1) = stat.zval;
[p, ~,stat] = ranksum(model2.weight, model4.weight(good_subs4));
param_comp.p_exp24(2) = p;
param_comp.z_exp24(2) = stat.zval;
[p, ~,stat] = ranksum(model2.exponent, model4.exponent(good_subs4));
param_comp.p_exp24(3) = p;
param_comp.z_exp24(3) = stat.zval;
[p, ~,stat] = ranksum(model2.K, model4.K(good_subs4));
param_comp.p_exp24(4) = p;
param_comp.z_exp24(4) = stat.zval;

[p, ~,stat] = ranksum(model3.offset(good_subs3), model4.offset(good_subs4));
param_comp.p_exp34(1) = p;
param_comp.z_exp34(1) = stat.zval;
[p, ~,stat] = ranksum(model3.weight(good_subs3), model4.weight(good_subs4));
param_comp.p_exp34(2) = p;
param_comp.z_exp34(2) = stat.zval;
[p, ~,stat] = ranksum(model3.exponent(good_subs3), model4.exponent(good_subs4));
param_comp.p_exp34(3) = p;
param_comp.z_exp34(3) = stat.zval;
[p, ~,stat] = ranksum(model3.K(good_subs3), model4.K(good_subs4));
param_comp.p_exp34(4) = p;
param_comp.z_exp34(4) = stat.zval;

%false discovery rate corrections
ps(1,:) = param_comp.p_exp12;
ps(2,:) = param_comp.p_exp13;
ps(3,:) = param_comp.p_exp14;
ps(4,:) = param_comp.p_exp23;
ps(5,:) = param_comp.p_exp24;
ps(6,:) = param_comp.p_exp34;

[h, ~,~,adj_p] = fdr_bh(ps, 0.05, 'pdep');

param_comp.h_exp12 = h(1,:);
param_comp.h_exp13 = h(2,:);
param_comp.h_exp14 = h(3,:);
param_comp.h_exp23 = h(4,:);
param_comp.h_exp24 = h(5,:);
param_comp.h_exp34 = h(6,:);

param_comp.adjp_exp12 = adj_p(1,:);
param_comp.adjp_exp13 = adj_p(2,:);
param_comp.adjp_exp14 = adj_p(3,:);
param_comp.adjp_exp23 = adj_p(4,:);
param_comp.adjp_exp24 = adj_p(5,:);
param_comp.adjp_exp34 = adj_p(6,:);

save param_comp.mat param_comp
save param_smry.mat param_smry

%for fun; does larger average TE mean smaller exponent?
dataConf1.te = abs(dataConf1.ha - abs(dataConf1.rot));
dataConf2.te = abs(dataConf2.ha - abs(dataConf2.rot));
dataConf3.te = abs(dataConf3.ha - abs(dataConf3.rot));
dataConf4.te = abs(dataConf4.ha - abs(dataConf4.rot));

avgte1 = nanmean(dataConf1.te(:, 73:312), 2);
avgte2 = nanmean(dataConf2.te(:, 33:240), 2);
avgte3 = nanmean(dataConf3.te(:, 31:180), 2);
avgte4 = nanmean(dataConf4.te(:, 31:180), 2);

avgte = [avgte1; avgte2; avgte3(good_subs3); avgte4(good_subs4)];

allexps = [model1.exponent'; model2.exponent'; model3.exponent(good_subs3)'; model4.exponent(good_subs4)'];
err_exp = figure; hold on;
plot(avgte1, model1.exponent, 'b.');
plot(avgte2, model2.exponent, 'r.');
plot(avgte3(good_subs3), model3.exponent(good_subs3), 'g.');
plot(avgte4(good_subs4), model4.exponent(good_subs4), 'm.');
%plot(avgte, allexps, '.');
legend({'Exp. 1', 'Exp. 2', 'Exp. 3', 'Exp. 4'});
coeffs = polyfit(avgte, allexps, 1);
besty = polyval(coeffs, avgte);
plot(avgte, besty, '-');
xlabel("Average TE"); ylabel("Exponent Parameter Value");
[r, p] = corr(avgte, allexps)
exportgraphics(err_exp, 'te_vs_exp.svg', 'BackgroundColor', 'none','ContentType', 'vector')
%% functions
%pads end of array with nan to 29
function out = fillLength(arr, len)
    padding = len - length(arr);
    out =  [arr nan(1, padding)];
end