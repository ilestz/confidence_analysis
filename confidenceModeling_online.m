%% quick pass modeling of Chris' confidence task
%% SDM; New Haven CT; 08/08/2022

clear; clc;
exp_num = 4;
if exp_num == 3
    load dataConf_online_gradual.mat; % load data frame; EXP1 or EXP2
    dataConf = dataConf_online_gradual; %EXP1 or EXP2
    good_subs = [1:17 19:25];
else
    load dataConf_online_zeromean.mat; % load data frame; EXP1 or EXP2
    dataConf = dataConf_online_zeromean; %EXP1 or EXP2
    good_subs = [1:13 15:23 25:29];
end
[nsubs,N] = size(dataConf.ha); % get dimensions

%ha and aim are always positive if in the opposite direction of rot; to
%compute the errors, we need to subtract the absolute value of rotation
dataConf.rot = -1 * abs(dataConf.rot); %negative abs of rot; err = ha + rot
%additionally, confidence is flipped, with 100 as low confidence and 0 as
%high confidence; we should fix this here
%dataConf.conf = 100-dataConf.conf;
%dataConf.conf(isnan(dataConf.conf)) = 50;

%quick performance metrics
% TODO: calculate percent compensation
% TODO: calculate STL on zero mean

%need a trial column to select only the reporting phase
dataConf.trial = repmat(1:200, 25, 1);
%% MAIN LOOP %%
for s = 1:length(good_subs)
    si = good_subs(s);
    subject = si; % get sub number
    rot_phase = dataConf.trial > 10 & dataConf.trial<=180; % fitting learning phase only (for now)
    rot_phase= rot_phase(1,:); 
    n_iter = 10; % fit iterations with X different random start values to avoid local minima
    
   confidence(si,:) = dataConf.conf(subject,rot_phase); % make a matrix just of confidence ratings
   
    %% model 1: confidence tracks last task error (aka task error one trial back); objective %%       
    disp(['now fitting subject ',num2str(si),' and 1TB TE']);
    
    for k = 1:n_iter
        
        offset = rand*100; %baseline confidence
        weight = rand*100; %weight/sensitivity to error
        
        params = [offset,weight];
        options=optimset('display','off');
        LB = [0   0];
        UB = [100 100];
        [params, error] = fmincon(@func_conf_1TB_TE,params,[],[],[],[],LB,UB,[],options,dataConf,rot_phase,subject);
        
        model1.p(k,:) = params;
        model1.error(k) = error;
    end     
    [modelFit_1TB_TE.ll(si),best] = min(model1.error);    
    modelFit_1TB_TE.offset(si) = model1.p(best,1);
    modelFit_1TB_TE.weight(si) = model1.p(best,2);
    modelFit_1TB_TE.n_params(si) = length(model1.p(1,:));
    [modelFit_1TB_TE.sse(si), modelFit_1TB_TE.confpred(si,:), modelFit_1TB_TE.r2(si)] = func_conf_1TB_TE(model1.p(best,:),dataConf,rot_phase,subject);
    modelFit_1TB_TE.bic(si) = sum(rot_phase)*log(modelFit_1TB_TE.sse(si)/sum(rot_phase)) + modelFit_1TB_TE.n_params(si)*log(sum(rot_phase));
    modelFit_1TB_TE.aic(si) = sum(rot_phase)*log(modelFit_1TB_TE.sse(si)/sum(rot_phase)) + 2*modelFit_1TB_TE.n_params(si);
        
    
%     plot(confidence(si,:),'k','linewidth',2); hold on; % plot
%     plot(modelFit_1TB_TE.confpred(si,:),'g');
%     pause(0.5);clf; % examine figure

%% model 2: confidence tracks last task error (aka task error one trial back); subjective %%       
    disp(['now fitting subject ',num2str(si),' and 1TB TE EXP']);
    
    for k = 1:n_iter
        
        offset = rand*100; %baseline confidence
        weight = rand*100; %weight/sensitivity to error
        te_exp = rand*4; %exponent on errors
        
        params = [offset,weight, te_exp];
        options=optimset('display','off');
        LB = [0   0    0];
        UB = [100 100  4];
        [params, error] = fmincon(@func_conf_1TB_TE_EXP,params,[],[],[],[],LB,UB,[],options,dataConf,rot_phase,subject);
        
        model2.p(k,:) = params;
        model2.error(k) = error;
    end     
    [modelFit_1TB_TE_EXP.ll(si),best] = min(model2.error);    
    modelFit_1TB_TE_EXP.offset(si) = model2.p(best,1);
    modelFit_1TB_TE_EXP.weight(si) = model2.p(best,2);
    modelFit_1TB_TE_EXP.exponent(si) = model2.p(best,3);
    modelFit_1TB_TE_EXP.n_params(si) = length(model2.p(1,:));
    [modelFit_1TB_TE_EXP.sse(si), modelFit_1TB_TE_EXP.confpred(si,:), modelFit_1TB_TE_EXP.r2(si)] = func_conf_1TB_TE_EXP(model2.p(best,:),dataConf,rot_phase,subject);
    modelFit_1TB_TE_EXP.bic(si) = sum(rot_phase)*log(modelFit_1TB_TE_EXP.sse(si)/sum(rot_phase)) + modelFit_1TB_TE_EXP.n_params(si)*log(sum(rot_phase));
    modelFit_1TB_TE_EXP.aic(si) = sum(rot_phase)*log(modelFit_1TB_TE_EXP.sse(si)/sum(rot_phase)) + 2*modelFit_1TB_TE_EXP.n_params(si);
    
%     plot(confidence(si,:),'k','linewidth',2); hold on; % plot
%     plot(modelFit_1TB_TE_EXP.confpred(si,:),'g');
%     pause(0.5);clf; % examine figure

    %% model 3: confidence tracks history of errors (HIST TE); objective %%       
    disp(['now fitting subject ',num2str(si),' and HIST TE']);
    
    for k = 1:n_iter
        
        offset = rand*100; %baseline confidence
        weight = rand*100; %weight/sensitivity to error
        K = rand; %learning rate for state estimate of errors
        
        params = [offset,weight, K];
        options=optimset('display','off');
        LB = [0   0    0];
        UB = [100 100  1];
        [params, error] = fmincon(@func_conf_HIST_TE,params,[],[],[],[],LB,UB,[],options,dataConf,rot_phase,subject);
        
        model3.p(k,:) = params;
        model3.error(k) = error;
    end     
    [modelFit_HIST_TE.ll(si),best] = min(model3.error);    
    modelFit_HIST_TE.offset(si) = model3.p(best,1);
    modelFit_HIST_TE.weight(si) = model3.p(best,2);
    modelFit_HIST_TE.K(si) = model3.p(best,3);
    modelFit_HIST_TE.n_params(si) = length(model3.p(1,:));
    [modelFit_HIST_TE.sse(si), modelFit_HIST_TE.confpred(si,:), modelFit_HIST_TE.r2(si)] = func_conf_HIST_TE(model3.p(best,:),dataConf,rot_phase,subject);
    modelFit_HIST_TE.bic(si) = sum(rot_phase)*log(modelFit_HIST_TE.sse(si)/sum(rot_phase)) + modelFit_HIST_TE.n_params(si)*log(sum(rot_phase));
    modelFit_HIST_TE.aic(si) = sum(rot_phase)*log(modelFit_HIST_TE.sse(si)/sum(rot_phase)) + 2*modelFit_HIST_TE.n_params(si);
    
%     plot(confidence(si,:),'k','linewidth',2); hold on; % plot
%     plot(modelFit_HIST_TE.confpred(si,:),'g');
%     pause(0.5);clf; % examine figure

%% model 4: confidence tracks history of errors (HIST TE); subjective %%       
    disp(['now fitting subject ',num2str(si),' and HIST TE EXP']);
    
    for k = 1:n_iter
        
        offset = rand*100; %baseline confidence
        weight = rand*100; %weight/sensitivity to error
        te_exp = rand*4; %exponent on errors
        K = rand; %learning rate for state estimate of errors
        
        params = [offset,weight, te_exp, K];
        options=optimset('display','off');
        LB = [0   0    0  0];
        UB = [100 100  4  1];
        [params, error] = fmincon(@func_conf_HIST_TE_EXP,params,[],[],[],[],LB,UB,[],options,dataConf,rot_phase,subject);
        
        model4.p(k,:) = params;
        model4.error(k) = error;
    end
    [modelFit_HIST_TE_EXP.ll(si),best] = min(model4.error);    
    modelFit_HIST_TE_EXP.offset(si) = model4.p(best,1);
    modelFit_HIST_TE_EXP.weight(si) = model4.p(best,2);
    modelFit_HIST_TE_EXP.exponent(si) = model4.p(best,3);
    modelFit_HIST_TE_EXP.K(si) = model4.p(best,4);
    modelFit_HIST_TE_EXP.n_params(si) = length(model4.p(1,:));
    [modelFit_HIST_TE_EXP.sse(si), modelFit_HIST_TE_EXP.confpred(si,:), modelFit_HIST_TE_EXP.r2(si)] = func_conf_HIST_TE_EXP(model4.p(best,:),dataConf,rot_phase,subject);
    modelFit_HIST_TE_EXP.bic(si) = sum(rot_phase)*log(modelFit_HIST_TE_EXP.sse(si)/sum(rot_phase)) + modelFit_HIST_TE_EXP.n_params(si)*log(sum(rot_phase));
    modelFit_HIST_TE_EXP.aic(si) = sum(rot_phase)*log(modelFit_HIST_TE_EXP.sse(si)/sum(rot_phase)) + 2*modelFit_HIST_TE_EXP.n_params(si);

    %set up for param recovery
    dataConf.confpred_HIST_TE_EXP(si, :) = modelFit_HIST_TE_EXP.confpred(si, :);
    % plot(confidence(si,:),'k','linewidth',2); hold on; % plot
    % plot(modelFit_HIST_TE_EXP.confpred(si,:),'g');
    % waitforbuttonpress;clf; % examine figure


end
% 
% %parameter recovery
% pr_rhos = nan(10, 4);
% best_errs = zeros(nsubs, 1);
% for k = 1:10
%     for si = 1:nsubs
%         offset = rand*100; %baseline confidence
%         weight = rand*100; %weight/sensitivity to error
%         te_exp = rand*4;
%         K = rand;
%         params = [offset,weight, te_exp, K];
%         [params, error] = fmincon(@func_conf_HIST_TE_EXP,params,[],[],[],[],LB,UB,[],options,dataConf,rot_phase,si, 1);
%         if error<best_errs(si) || k == 1 %only keep the best fitting params so far
%             best_errs(si) = error;
%             modelFit_HIST_TE_EXP.pr_offset(si) = params(1);
%             modelFit_HIST_TE_EXP.pr_weight(si) = params(2);
%             modelFit_HIST_TE_EXP.pr_exp(si) = params(3);
%             modelFit_HIST_TE_EXP.pr_K(si) = params(4);
%         end
%     end
%     pr_rhos(k, 1) = corr(modelFit_HIST_TE_EXP.pr_offset', modelFit_HIST_TE_EXP.offset');
%     pr_rhos(k, 2) = corr(modelFit_HIST_TE_EXP.pr_weight', modelFit_HIST_TE_EXP.weight');
%     pr_rhos(k, 3) = corr(modelFit_HIST_TE_EXP.pr_exp', modelFit_HIST_TE_EXP.exponent');
%     pr_rhos(k, 4) = corr(modelFit_HIST_TE_EXP.pr_K', modelFit_HIST_TE_EXP.K');
% end
% param_names = ["Conf_{max}" "\eta" "\gamma" "\alpha"];
% for p = 1:4
%     figure;
%     plot(pr_rhos(:, p), 'o');
%     xlabel("Iterations");
%     ylabel("Pearson correlation");
%     title(param_names(p));
% end

%%
%compute the true R2 (above is merely the pseudo R2)
for si = 1:length(good_subs)
    s = good_subs(si);
    true_r = corrcoef(modelFit_1TB_TE.confpred(s,:)', dataConf.conf(s, rot_phase)', 'rows', 'complete');
    modelFit_1TB_TE.true_r2(s) = true_r(1,2)^2;
    true_r = corrcoef(modelFit_1TB_TE_EXP.confpred(s,:)', dataConf.conf(s, rot_phase)', 'rows', 'complete');
    modelFit_1TB_TE_EXP.true_r2(s) = true_r(1,2)^2;
    true_r = corrcoef(modelFit_HIST_TE.confpred(s,:)', dataConf.conf(s, rot_phase)', 'rows', 'complete');
    modelFit_HIST_TE.true_r2(s) = true_r(1,2)^2;
    true_r = corrcoef(modelFit_HIST_TE_EXP.confpred(s,:)', dataConf.conf(s, rot_phase)', 'rows', 'complete');
    modelFit_HIST_TE_EXP.true_r2(s) = true_r(1,2)^2;
end
%%
%add stuff into the model fitting that is useful for figures
modelFit_HIST_TE_EXP.aic_1TB_TE = modelFit_1TB_TE.aic - modelFit_HIST_TE_EXP.aic;
modelFit_HIST_TE_EXP.aic_1TB_TE_EXP = modelFit_1TB_TE_EXP.aic - modelFit_HIST_TE_EXP.aic;
modelFit_HIST_TE_EXP.aic_HIST_TE = modelFit_HIST_TE.aic - modelFit_HIST_TE_EXP.aic;
modelFit_HIST_TE_EXP.aic_HIST_TE_EXP = zeros(1, length(modelFit_1TB_TE.aic));
modelFit_HIST_TE_EXP.aic_diff = modelFit_HIST_TE_EXP.aic - min(modelFit_HIST_TE.aic, min(modelFit_1TB_TE_EXP.aic,  modelFit_1TB_TE.aic));

% 
% if exp_num == 3
%     save modelFit_1TB_TE3.mat modelFit_1TB_TE
%     save modelFit_1TB_TE_EXP3.mat modelFit_1TB_TE_EXP
%     save modelFit_HIST3.mat modelFit_HIST_TE
      %save modelFit_HIST_TE_EXP3.mat modelFit_HIST_TE_EXP
% else
%     save modelFit_1TB_TE4.mat modelFit_1TB_TE
%     save modelFit_1TB_TE_EXP4.mat modelFit_1TB_TE_EXP
%     save modelFit_HIST4.mat modelFit_HIST_TE
%     save modelFit_HIST_TE_EXP4.mat modelFit_HIST_TE_EXP
% end
% 
% if ispc
% !shutdown -s -f -t 0
% else
% !shutdown -h now
% end
%%
%show avg model fits
gcf_fit = figure;
c=plot(nanmean(confidence),'k','linewidth',3); hold on; % plot
shadedErrorPlot(confidence, 'k', 0.5,0);
%c1=plot(nanmean(modelFit_1TB_TE.confpred),'g', 'linewidth', 2); % plot fit
%c2=plot(nanmean(modelFit_1TB_TE_EXP.confpred), 'b', 'linewidth', 2); %plot fit
%c3=plot(nanmean(modelFit_HIST_TE.confpred),'r', 'linewidth', 2); % plot fit
model_conf = modelFit_HIST_TE_EXP.confpred(:, 2:end);
shadedErrorPlot(model_conf, 'm', 0.5, 1);
c4=plot([NaN nanmean(modelFit_HIST_TE_EXP.confpred(:, 2:end))], 'm', 'linewidth', 2); %plot fit
xlabel('Trial Number'); ylabel('Confidence (1-100)');
legend([c c4], {"Average Confidence",...
    sprintf('ESS_{Subj} (R^{2}= %.2f)', nanmean(modelFit_HIST_TE_EXP.true_r2))});
%exportgraphics(gcf_fit, 'modelFit_exp2.eps', 'ContentType', 'vector');
%%
%aic comparison
gcf_aic = figure;
aic_1TB_TE = sum(modelFit_1TB_TE.aic) - sum(modelFit_HIST_TE_EXP.aic);
aic_1TB_TE_EXP = sum(modelFit_1TB_TE_EXP.aic) - sum(modelFit_HIST_TE_EXP.aic);
aic_HIST_TE = sum(modelFit_HIST_TE.aic) - sum(modelFit_HIST_TE_EXP.aic);
aic_HIST_TE_EXP = 0;
bar(1:4, [aic_HIST_TE_EXP aic_HIST_TE aic_1TB_TE_EXP aic_1TB_TE]); xticks(1:4);
xticklabels({'ESS_{Subj}', 'ESS_{Obj}', 'OTB_{Subj}', 'OTB_{Obj}'}); ylabel("\Delta AIC");
xlabel("Model");
%exportgraphics(gcf_aic,sprintf('modelfitbic_exp%d_vec.eps', exp_num),'BackgroundColor','none','ContentType','vector')

%aic second-best
gcf_aic2 = figure;
aic_diff = modelFit_HIST_TE_EXP.aic - min(modelFit_HIST_TE.aic, min(modelFit_1TB_TE_EXP.aic,  modelFit_1TB_TE.aic));
aic_diff = sort(aic_diff);
bar(aic_diff); xlabel("Subject"); ylabel("\Delta AIC");
%exportgraphics(gcf_aic2, sprintf('aic_deltas_exp%d_vec.eps', exp_num), 'BackgroundColor', 'none','ContentType', 'vector');


%% across experiment parameter comparisons
%load the data frames
% load modelFit_HIST_TE_EXP.mat
model1 = modelFit_HIST_TE_EXP;
figure;
violinplot(model1.offset);
ylim([0 100]);
ylabel("Offset");

figure;
violinplot(model1.exponent);
ylabel("Exponent");

figure;
violinplot(model1.weight);
ylabel("Weight");

figure;
violinplot(model1.K);
ylabel("Learning Rate");



function shadedErrorPlot(data, color, alpha, offset)
    std_dat = nanstd(data);
    sem_dat = std_dat/sqrt(size(data,1));
    mean_dat = nanmean(data);
    x_vec = [offset + (1:length(data)) offset+fliplr(1:length(data))];
    patch = fill(x_vec, [mean_dat+sem_dat,fliplr(mean_dat-sem_dat)], color);
    set(patch, 'edgecolor', 'none');
    set(patch, 'FaceAlpha', alpha);
end








