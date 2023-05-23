%% sdm quick pass on confidence vmr task
%% May 17 2022; New Haven CT
clear;close all;clc;

%df = readtable('confidence-EXP1_IMV.csv');
df = readtable('confidence_online_zeromean.csv');

subs = unique(df.id);
nsubs = length(subs);

num_excluded = 0;

for j = 1:nsubs
    sidx = df.id == subs(j);
    signs = -sign((df.manipulation_angle(sidx)));
    signs(signs==0) = signs(100);
    ha(j,:) = df.diff_angle(sidx);
    for i = 1:length(ha(j,:))
        if ha(j,i) > 180
            ha(j,i) = ha(j,i) - 360;
        end
    end
    ha(j,:) = ha(j,:).*signs';
    %ha(j,:) = ha(j,:)*-1;
    confidence(j,:) = df.judgement_width(sidx);
    ttype = df.is_judged(sidx);    
    
    num_excluded = num_excluded + sum(confidence(j,:)==50);
    %confidence(j,confidence(j,:)==50) = nan;
   
    rot(j,:) = df.manipulation_angle(sidx);
    
    
end

% save data
dataConf_online_zeromean.ha = ha;
dataConf_online_zeromean.conf = confidence;
dataConf_online_zeromean.rot = rot;

%rot_exp1 = rot;
save dataConf_online_zeromean dataConf_online_zeromean

% figure;
% for k = 1:nsubs
%     subplot(4,4,k);hold on;
%     plot(ha(k,:));
%     plot(aim(k,:));
%     plot(implicit(k,:));
%     plot(confidence(k,:)/2,'linewidth',4);
%     plot([0 length(tmp)],[30 30],'-k');       
%     xlabel('trial');
%     ylabel('angle; confidence/60');
%     legend('ha','aim','imp','conf');
% end
% figure;
% %plot(1:length(rot_exp1), rot_exp1, 'linewidth', 2); hold on;
% plot(1:length(rot), rot, 'linewidth', 2);
% xlabel("Trial Number"); ylabel("Perturbation (\circ)");
% %ylim([-65 65]);
% %ylim([-5 35]);
% %legend({"Experiment 1", "Experiment 2"});
% 
% 
% gcf_ha_rot = figure; %error bars!!
% ha_plot = plot(nanmean(ha), 'b', 'linewidth', 2);hold on;
% shadedErrorPlot(ha, 'b', 0.5, 0);
% rot_plot = plot(rot, 'k', 'linewidth', 2);
% %plot(nanmean(implicit), 'linewidth', 2);
% legend([ha_plot rot_plot], {'Hand Angle (\circ)', 'Perturbation (\circ)'});
% xlabel("Trial Number"); %ylim([-65 65]);
% exportgraphics(gcf_ha_rot, 'ha_rot_exp1.eps', 'ContentType', 'vector');
% 
% figure; % error bars!!
% conf_plot = plot(17:240, 100-nanmean(confidence(:, 17:240)), 'k', 'linewidth', 2);
% xlabel("Trial Number"); ylabel("Confidence (1-100)");
% legend('Confidence (1-100)'); ylim([0 100]);

function shadedErrorPlot(data, color, alpha, offset)
    std_dat = nanstd(data);
    sem_dat = std_dat/sqrt(size(data,1));
    mean_dat = nanmean(data);
    x_vec = [offset + (1:length(data)) offset+fliplr(1:length(data))];
    patch = fill(x_vec, [mean_dat-sem_dat,fliplr(mean_dat+sem_dat)], color);
    set(patch, 'edgecolor', 'none');
    set(patch, 'FaceAlpha', alpha);
end

