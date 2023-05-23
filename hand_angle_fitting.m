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

%% MAIN LOOP %%
for si = 1:nsubs
    
    subject = si; % get sub number
    rot_phase = dataConf.conf~=100; % fitting learning phase only (for now)
    rot_phase= rot_phase(1,:); 
    %rot_phase(17:218) = 0;
    %rot_phase = [zeros(218,1);ones(92,1); zeros(50,1)];
    n_iter = 5; % fit iterations with X different random start values to avoid local minima
    
   hand_angle(si,:) = dataConf.ha(subject,rot_phase); % make a matrix just of confidence ratings
   
    %% model 1: basic state space model %%       
    disp(['now fitting subject ',num2str(si),' and SS']);
    
    for k = 1:n_iter
        
        A = rand; %retention
        B = rand; %error sensitivity/ learning rate
        
        params = [A,B];
        options=optimset('display','off');
        LB = [0 0];
        UB = [1 1];
        [params, error] = fmincon(@fitSS,params,[],[],[],[],LB,UB,[],options,dataConf,rot_phase,subject);
        
        model1.p(k,:) = params;
        model1.error(k) = error;
    end     
    [modelSS.ll(si),best] = min(model1.error);    
    modelSS.A(si) = model1.p(best,1);
    modelSS.B(si) = model1.p(best,2);
    modelSS.n_params(si) = length(model1.p(1,:));
    [modelSS.sse(si), modelSS.hapred(si,:), modelSS.r2(si)] = fitSS(model1.p(best,:),dataConf,rot_phase,subject);
    modelSS.bic(si) = sum(rot_phase)*log(modelSS.sse(si)/sum(rot_phase)) + modelSS.n_params(si)*log(sum(rot_phase));
    modelSS.aic(si) = sum(rot_phase)*log(modelSS.sse(si)/sum(rot_phase)) + 2*modelSS.n_params(si);
    
%      plot(hand_angle(si,:),'k','linewidth',2); hold on; % plot
%      plot(modelSS.hapred(si,:),'b');
%      pause(0.5);clf; % examine figure
     
    %% model 2: high confidence results in greater error sensitivity %%       
    disp(['now fitting subject ',num2str(si),' and SS_pos']);
    
    for k = 1:n_iter
        
        A = rand; %retention
        B_m = rand; %learning rate slope
        B_b = rand; %learning rate intercept
        
        params = [A,B_m, B_b];
        options=optimset('display','off');
        LB = [0 -50 -50];
        UB = [1 50 50];
        [params, error] = fmincon(@fitSS_pos, params, [], [], [], [], LB, UB, [], options, dataConf, rot_phase, subject);
        
        model2.p(k,:) = params;
        model2.error(k) = error;
    end     
    [modelSS_pos.ll(si),best] = min(model2.error);    
    modelSS_pos.A(si) = model2.p(best,1);
    modelSS_pos.B_m(si) = model2.p(best,2);
    modelSS_pos.B_b(si) = model2.p(best,3);
    modelSS_pos.n_params(si) = length(model2.p(1,:));
    [modelSS_pos.sse(si), modelSS_pos.hapred(si,:), modelSS_pos.r2(si)] = fitSS_pos(model2.p(best,:),dataConf,rot_phase,subject);
    modelSS_pos.bic(si) = sum(rot_phase)*log(modelSS_pos.sse(si)/sum(rot_phase)) + modelSS_pos.n_params(si)*log(sum(rot_phase));
    modelSS_pos.aic(si) = sum(rot_phase)*log(modelSS_pos.sse(si)/sum(rot_phase)) + 2*modelSS_pos.n_params(si);
    
     % plot(hand_angle(si,:),'k','linewidth',2); hold on; % plot
     % plot(modelSS_pos.hapred(si,:),'r');
     % pause(0.5);clf; % examine figure
    
    %% model 3: high uncertainty results in greater error sensitivity %%       
    disp(['now fitting subject ',num2str(si),' and SS_neg']);
    
    for k = 1:n_iter
        
        A = rand; %retention
        B_m = rand; %learning rate slope
        B_b = rand; %learning rate intercept
        
        params = [A,B_m, B_b];
        options=optimset('display','off');
        LB = [0 0 0];
        UB = [1 1 1];
        [params, error] = fmincon(@fitSS_neg,params,[],[],[],[],LB,UB,[],options,dataConf,rot_phase,subject);
        
        model3.p(k,:) = params;
        model3.error(k) = error;
    end     
    [modelSS_neg.ll(si),best] = min(model3.error);    
    modelSS_neg.A(si) = model3.p(best,1);
    modelSS_neg.B_m(si) = model3.p(best,2);
    modelSS_neg.B_b(si) = model3.p(best,3);
    modelSS_neg.n_params(si) = length(model3.p(1,:));
    [modelSS_neg.sse(si), modelSS_neg.hapred(si,:), modelSS_neg.r2(si)] = fitSS_neg(model3.p(best,:),dataConf,rot_phase,subject);
    modelSS_neg.bic(si) = sum(rot_phase)*log(modelSS_neg.sse(si)/sum(rot_phase)) + modelSS_neg.n_params(si)*log(sum(rot_phase));
    modelSS_neg.aic(si) = sum(rot_phase)*log(modelSS_neg.sse(si)/sum(rot_phase)) + 2*modelSS_neg.n_params(si);
    
%      plot(hand_angle(si,:),'k','linewidth',2); hold on; % plot
%      plot(modelSS_neg.hapred(si,:),'g');
%      pause(0.5);clf; % examine figure
    
    %% model 4: low and high confidence results in greater error sensitivity %%       
    disp(['now fitting subject ',num2str(si),' and SS_abs']);
    
    for k = 1:n_iter
        
        A_m = rand; %retention
        A_b = rand;
        B_m = rand; %learning rate slope
        B_b = rand; %learning rate intercept
        
        params = [A_m, A_b,B_m, B_b];
        options=optimset('display','off');
        LB = [-50 -50 -50 -50];
        UB = [50 50 50 50];
        [params, error] = fmincon(@fitSS_abs,params,[],[],[],[],LB,UB,[],options,dataConf,rot_phase,subject);
        
        model4.p(k,:) = params;
        model4.error(k) = error;
    end     
    [modelSS_abs.ll(si),best] = min(model4.error);    
    modelSS_abs.A_m(si) = model4.p(best,1);
    modelSS_abs.A_b(si) = model4.p(best,2);
    modelSS_abs.B_m(si) = model4.p(best,3);
    modelSS_abs.B_b(si) = model4.p(best,4);
    modelSS_abs.n_params(si) = length(model4.p(1,:));
    [modelSS_abs.sse(si), modelSS_abs.hapred(si,:), modelSS_abs.r2(si)] = fitSS_abs(model4.p(best,:),dataConf,rot_phase,subject);
    modelSS_abs.bic(si) = sum(rot_phase)*log(modelSS_abs.sse(si)/sum(rot_phase)) + modelSS_abs.n_params(si)*log(sum(rot_phase));
    modelSS_abs.aic(si) = sum(rot_phase)*log(modelSS_abs.sse(si)/sum(rot_phase)) + 2*modelSS_abs.n_params(si);
    
%      plot(hand_angle(si,:),'k','linewidth',2); hold on; % plot
%      plot(modelSS_abs.hapred(si,:),'m');
%      pause(0.5);clf; % examine figure

    %% model 5: high confidence tunes up adaptation %%       
    disp(['now fitting subject ',num2str(si),' and SS_bias']);
    
    for k = 1:n_iter
        
        A = rand; %retention
        B = rand; %learning rate
        B_high = rand;
        
        params = [A,B, B_high];
        options=optimset('display','off');
        LB = [0 0 0];
        UB = [1 1 1];
        [params, error] = fmincon(@fitSS_bias,params,[],[],[],[],LB,UB,[],options,dataConf,rot_phase,subject);
        
        model5.p(k,:) = params;
        model5.error(k) = error;
    end     
    [modelSS_bias.ll(si),best] = min(model5.error);    
    modelSS_bias.A(si) = model5.p(best,1);
    modelSS_bias.B(si) = model5.p(best,2);
    modelSS_bias.B_high(si) = model5.p(best, 3);
    modelSS_bias.n_params(si) = length(model5.p(1,:));
    [modelSS_bias.sse(si), modelSS_bias.hapred(si,:), modelSS_bias.r2(si)] = fitSS_bias(model5.p(best,:),dataConf,rot_phase,subject);
    modelSS_bias.bic(si) = sum(rot_phase)*log(modelSS_bias.sse(si)/sum(rot_phase)) + modelSS_bias.n_params(si)*log(sum(rot_phase));
    modelSS_bias.aic(si) = sum(rot_phase)*log(modelSS_bias.sse(si)/sum(rot_phase)) + 2*modelSS_bias.n_params(si);
    
%      plot(hand_angle(si,:),'k','linewidth',2); hold on; % plot
%      plot(modelSS_abs.hapred(si,:),'m');
%      pause(0.5);clf; % examine figure
    
end

gcf_fit = figure;
h=plot(nanmean(hand_angle),'k','linewidth',3); hold on; % plot
%c1=plot(nanmean(modelFit_1TB_TE.confpred),'g', 'linewidth', 2); % plot fit
%c2=plot(nanmean(modelFit_1TB_TE_EXP.confpred), 'b', 'linewidth', 2); %plot fit
%c3=plot(nanmean(modelFit_HIST_TE.confpred),'r', 'linewidth', 2); % plot fit
m1 = plot(nanmean(modelSS.hapred), 'b', 'linewidth', 1); %plot fit
m2 = plot(nanmean(modelSS_pos.hapred), 'r', 'linewidth', 1); %plot fit
m3 = plot(nanmean(modelSS_neg.hapred), 'g', 'linewidth', 1); %plot fit
m4 = plot(nanmean(modelSS_bias.hapred), 'm', 'linewidth', 1); %plot fit

xlabel('Trial Number'); ylabel('HA');
%legend([h m1 m2 m3 m4], {"Average HA", 'SS'});

%heatmap of best fitting model
modelaics = [modelSS.aic; modelSS_pos.aic; modelSS_neg.aic; modelSS_bias.aic];
modelaics = modelaics - repmat(min(modelaics),4,1);
modelaics(modelaics>0) = 1;
figure; heatmap(modelaics);