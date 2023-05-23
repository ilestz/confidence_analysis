%across experiment parameter comparisons; not yet set up to load online exps
%load the data frames
% load modelFit_HIST_TE_EXP.mat
% model1 = modelFit_HIST_TE_EXP;
% load modelFit_HIST_TE_EXP2.mat
% model2 = modelFit_HIST_TE_EXP;
%offset
gcf_offset = figure;
exp = [repmat(cellstr("Exp. 1"), 18,1);...
    repmat(cellstr("Exp. 2"), 20,1)];
expOffsets = [model1.offset model2.offset];
violinplot(expOffsets, exp, 'GroupOrder', {'Exp. 1', 'Exp. 2'});
ylabel("Conf_{Max} (1-100)");
%exportgraphics(gcf_offset,'expoffsets_vec.eps','ContentType','vector');
%weight
gcf_weight = figure;
expWeights = [model1.weight model2.weight];
violinplot(expWeights, exp, 'GroupOrder', {'Exp. 1', 'Exp. 2'});
ylabel("\eta", 'FontSize', 20);
%exportgraphics(gcf_weight,'expweights_vec.eps','ContentType','vector');
%exponent
gcf_exponent = figure;
expExp = [model1.exponent model2.exponent];
violinplot(expExp, exp, 'GroupOrder', {'Exp. 1', 'Exp. 2'});
ylabel("\gamma", 'FontSize', 20);
%exportgraphics(gcf_exponent,'expexponent_vec.eps','ContentType','vector');
%learning rate
gcf_K = figure;
expKs = [model1.K model2.K];
violinplot(expKs, exp, 'GroupOrder', {'Exp. 1', 'Exp. 2'});
ylabel("\alpha", 'FontSize', 20);
%exportgraphics(gcf_K,'expalphas_vec.eps','ContentType','vector');



%parameter comparison
%offset
figure;
modelOffsets = [modelFit_1TB_TE.offset modelFit_1TB_TE_EXP.offset modelFit_HIST_TE.offset modelFit_HIST_TE_EXP.offset]';
violinplot(modelOffsets, model, 'GroupOrder', {'1TB TE', '1TB TE Subj', 'Hist TE', 'Hist TE Subj'}, 'ViolinColor', ...
    [0 1 0; 0 0 1; 1 0 0; 1 0 1]);
xlabel("Model"); ylabel("Offset");

%weight
figure;
modelWeights = [modelFit_1TB_TE.weight modelFit_1TB_TE_EXP.weight modelFit_HIST_TE.weight modelFit_HIST_TE_EXP.weight]';
violinplot(modelWeights, model, 'GroupOrder', {'1TB TE', '1TB TE Subj', 'Hist TE', 'Hist TE Subj'}, 'ViolinColor', ...
    [0 1 0; 0 0 1; 1 0 0; 1 0 1]);
xlabel("Model"); ylabel("Weight");

%exponent
figure;
model_exp = [repmat(cellstr('1TB TE Subj'), nsubs, 1); ...
    repmat(cellstr('Hist TE Subj'), nsubs, 1)];
modelExps = [modelFit_1TB_TE_EXP.exponent modelFit_HIST_TE_EXP.exponent]';
violinplot(modelExps, model_exp, 'GroupOrder', {'1TB TE Subj', 'Hist TE Subj'}, 'ViolinColor', ...
    [0 0 1; 1 0 1]);
xlabel("Model"); ylabel("Exponent");

%learning rate
figure;
model_k = [repmat(cellstr('Hist TE'), nsubs, 1); ...
    repmat(cellstr('Hist TE Subj'), nsubs, 1)];
modelKs = [modelFit_HIST_TE.K modelFit_HIST_TE_EXP.K]';
violinplot(modelKs, model_k, 'GroupOrder', {'Hist TE', 'Hist TE Subj'}, 'ViolinColor', ...
    [1 0 0; 1 0 1]);
xlabel("Model"); ylabel("Learning Rate");


function shadedErrorPlot(data, color, alpha, offset)
    std_dat = nanstd(data);
    sem_dat = std_dat/sqrt(size(data,1));
    mean_dat = nanmean(data);
    x_vec = [offset + (1:length(data)) offset+fliplr(1:length(data))];
    patch = fill(x_vec, [mean_dat+sem_dat,fliplr(mean_dat-sem_dat)], color);
    set(patch, 'edgecolor', 'none');
    set(patch, 'FaceAlpha', alpha);
end