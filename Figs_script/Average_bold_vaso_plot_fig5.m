%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Script to create plots of average BOLD and VASO layer profiles 
% for each hippocampal subfield

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define paths to layer-wise GLM results (outputs of GLM_layers_normalized.m)
fileContent = fileread('unique_no_PA.txt');
folders = splitlines(fileContent);
folders = folders(1:end-1); % remove last empty entry
folders = regexp(folders, '^(.*?/func)/', 'tokens', 'once');
folders = vertcat(folders{:});

%% load constant variables which do not change across subjects and subregions 
Subfield = 4; % number of subfields
N = 30; % number of layers 
depths = 1:N;
names = {'SRLM','inner','outer'};
x_SRLM(1,1:30) = 10; % 10th bin is where SRLM ends
x_inner(1,1:30) = 20; % 20th bin where inner surface ends
N_sub = 6;
%% Compute the mean and SEM values for BOLD and VASO for the contrast memory vs math
%starting with BOLD

for n = 1:N_sub
    cd(folders{n})
    cd('first_level/bold')
    for m = 1:Subfield
        RESULT_Z{m,1} = strcat('results_Ztransformed','_',string(m),'.mat');
        ARRAYS_Z{n,m} = load(RESULT_Z{m,1});                
    end
    cd ../../../../     
end

for m = 1:length(ARRAYS_Z) % Sub = 1; ca1 = 2; ca2 = 3; ca3 = 4, subfields labels
    ALL_cont_sub_Z(m,1:30) = ARRAYS_Z{m,1}.con_array;
    ALL_cont_ca1_Z(m,1:30) = ARRAYS_Z{m,2}.con_array;
    ALL_cont_ca2_Z(m,1:30) = ARRAYS_Z{m,3}.con_array;  
    ALL_cont_ca3_Z(m,1:20) = ARRAYS_Z{m,4}.con_array(11:30); % ca3 has 20 bins due to excluding srlm section
end

mean_sub_Z = mean(ALL_cont_sub_Z);
sem_sub_Z = std(ALL_cont_sub_Z)./sqrt(N_sub);
upper_sub = mean_sub_Z + sem_sub_Z;
lower_sub = mean_sub_Z - sem_sub_Z;


mean_ca1_Z = mean(ALL_cont_ca1_Z);
sem_ca1_Z = std(ALL_cont_ca1_Z)./sqrt(N_sub);
upper_ca1 = mean_ca1_Z + sem_ca1_Z;
lower_ca1 = mean_ca1_Z - sem_ca1_Z;


mean_ca2_Z = mean(ALL_cont_ca2_Z);
sem_ca2_Z = std(ALL_cont_ca2_Z)./sqrt(N_sub);
upper_ca2 = mean_ca2_Z + sem_ca2_Z;
lower_ca2 = mean_ca2_Z - sem_ca2_Z;


mean_ca3_Z = mean(ALL_cont_ca3_Z);
sem_ca3_Z = std(ALL_cont_ca3_Z)./sqrt(N_sub);
upper_ca3 = mean_ca3_Z + sem_ca3_Z;
lower_ca3 = mean_ca3_Z - sem_ca3_Z;

% now repeat the same for VASO

for n = 1:N_sub
    cd(folders{n})
    cd('first_level/vaso')
    for m = 1:Subfield
        RESULT_Z_v{m,1} = strcat('results_Ztransformed','_',string(m),'.mat');
        ARRAYS_Z_v{n,m} = load(RESULT_Z_v{m,1});                   
    end
    cd ../../../../ 
end

       
for m = 1:length(ARRAYS_Z_v)    
    ALL_cont_sub_Z_v(m,1:30) = ARRAYS_Z_v{m,1}.con_array;    
    ALL_cont_ca1_Z_v(m,1:30) = ARRAYS_Z_v{m,2}.con_array;    
    ALL_cont_ca2_Z_V(m,1:30) = ARRAYS_Z_v{m,3}.con_array;    
    ALL_cont_ca3_Z_v(m,1:20) = ARRAYS_Z_v{m,4}.con_array(11:30);
end

mean_sub_Z_v = mean(ALL_cont_sub_Z_v);
sem_sub_Z_v = std(ALL_cont_sub_Z_v)./sqrt(N_sub);
upper_sub_v = mean_sub_Z_v + sem_sub_Z_v;
lower_sub_v = mean_sub_Z_v - sem_sub_Z_v;


mean_ca1_Z_v = mean(ALL_cont_ca1_Z_v);
sem_ca1_Z_v = std(ALL_cont_ca1_Z_v)./sqrt(N_sub);
upper_ca1_v = mean_ca1_Z_v + sem_ca1_Z_v;
lower_ca1_v = mean_ca1_Z_v - sem_ca1_Z_v;


mean_ca2_Z_v = mean(ALL_cont_ca2_Z_V);
sem_ca2_Z_v = std(ALL_cont_ca2_Z_V)./sqrt(N_sub);
upper_ca2_v = mean_ca2_Z_v + sem_ca2_Z_v;
lower_ca2_v = mean_ca2_Z_v - sem_ca2_Z_v;


mean_ca3_Z_v = mean(ALL_cont_ca3_Z_v);
sem_ca3_Z_v = std(ALL_cont_ca3_Z_v)./sqrt(N_sub);
upper_ca3_v = mean_ca3_Z_v + sem_ca3_Z_v;
lower_ca3_v = mean_ca3_Z_v - sem_ca3_Z_v;

%% plot the average profiles

y_threshold_Z = linspace(-0.2,0.6,30);

figure('Color','w');
set(gcf, 'GraphicsSmoothing', 'off');
set(gcf, 'Renderer', 'painters');  
subplot(2,2,1)
xlim([1 30]);hold on
plot(x_SRLM,y_threshold_Z,'LineWidth',2,'Color','k','LineStyle','--');hold on
set(gca,'xtick',[5,15,25],'xticklabel',names,'Color','w')
plot(x_inner,y_threshold_Z,'LineWidth',2,'Color','k','LineStyle','--');

ylabel('memory-math contrast [z]');
%axis('square');

box on
ax = gca;
ax.LineWidth = 2;
ax.FontSize = 15; 
fill([depths, fliplr(depths)], [upper_sub, fliplr(lower_sub)], [0 0 0.8], ...
     'FaceAlpha', 0.2, 'EdgeColor', 'none');
H1 = plot(depths,mean_sub_Z,'LineWidth',3,'Color','b'); 
fill([depths, fliplr(depths)], [upper_sub_v, fliplr(lower_sub_v)], [0.8 0 0], ...
     'FaceAlpha', 0.2, 'EdgeColor', 'none');
H2 = plot(depths,mean_sub_Z_v,'LineWidth',3,'Color','r'); 
title('Subiculum'); 
hleglines = [H1(1),H2(1)];
Hleg = legend(hleglines,'BOLD','VASO','Location','northeast');hold on

subplot(2,2,2)
xlim([1 30]);hold on
plot(x_SRLM,y_threshold_Z,'LineWidth',2,'Color','k','LineStyle','--');hold on
set(gca,'xtick',[5,15,25],'xticklabel',names,'Color','w')
plot(x_inner,y_threshold_Z,'LineWidth',2,'Color','k','LineStyle','--');

ylabel('memory-math contrast [z]');
%axis('square');

box on
ax = gca;
ax.LineWidth = 2;
ax.FontSize = 15; 
fill([depths, fliplr(depths)], [upper_ca1, fliplr(lower_ca1)], [0 0 0.8], ...
     'FaceAlpha', 0.2, 'EdgeColor', 'none');
H3 = plot(depths,mean_ca1_Z,'LineWidth',3,'Color','b'); 
fill([depths, fliplr(depths)], [upper_ca1_v, fliplr(lower_ca1_v)], [0.8 0 0], ...
     'FaceAlpha', 0.2, 'EdgeColor', 'none');
H4 = plot(depths,mean_ca1_Z_v,'LineWidth',3,'Color','r'); 
title('CA1')
hleglines = [H3(1),H4(1)];
Hleg = legend(hleglines,'BOLD','VASO','Location','northeast');hold on


subplot(2,2,3)
xlim([1 30]);hold on
plot(x_SRLM,y_threshold_Z,'LineWidth',2,'Color','k','LineStyle','--');hold on
set(gca,'xtick',[5,15,25],'xticklabel',names,'Color','w')
plot(x_inner,y_threshold_Z,'LineWidth',2,'Color','k','LineStyle','--');

ylabel('memory-math contrast [z]');
%axis('square');

box on
ax = gca;
ax.LineWidth = 2;
ax.FontSize = 15; 
fill([depths, fliplr(depths)], [upper_ca2, fliplr(lower_ca2)], [0 0 0.8], ...
     'FaceAlpha', 0.2, 'EdgeColor', 'none');
H5 = plot(depths,mean_ca2_Z,'LineWidth',3,'Color','b'); 
fill([depths, fliplr(depths)], [upper_ca2_v, fliplr(lower_ca2_v)], [0.8 0 0], ...
     'FaceAlpha', 0.2, 'EdgeColor', 'none');
H6 = plot(depths,mean_ca2_Z_v,'LineWidth',3,'Color','r'); 
title('CA2'); 
hleglines = [H5(1),H6(1)];
Hleg = legend(hleglines,'BOLD','VASO','Location','northeast');hold on

subplot(2,2,4)
xlim([1 30]);hold on
plot(x_SRLM,y_threshold_Z,'LineWidth',2,'Color','k','LineStyle','--');hold on
set(gca,'xtick',[5,15,25],'xticklabel',names,'Color','w')
plot(x_inner,y_threshold_Z,'LineWidth',2,'Color','k','LineStyle','--');

ylabel('memory-math contrast [z]');
%axis('square');
box on
ax = gca;
ax.LineWidth = 2;
ax.FontSize = 15; 
fill([depths(11:30), fliplr(depths(11:30))], [upper_ca3, fliplr(lower_ca3)], [0 0 0.8], ...
     'FaceAlpha', 0.2, 'EdgeColor', 'none');
H7 = plot(depths(11:30),mean_ca3_Z,'LineWidth',3,'Color','b'); 
fill([depths(11:30), fliplr(depths(11:30))], [upper_ca3_v, fliplr(lower_ca3_v)], [0.8 0 0], ...
     'FaceAlpha', 0.2, 'EdgeColor', 'none');
H8 = plot(depths(11:30),mean_ca3_Z_v,'LineWidth',3,'Color','r'); 
title('CA3')
hleglines = [H7(1),H8(1)];
Hleg = legend(hleglines,'BOLD','VASO','Location','northeast');hold on

% save the figure
set(gcf, 'Renderer', 'opengl');
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [2, 2, 29.7, 21]); % Set figure units and size (A4 landscape: 29.7 cm × 21 cm)

exportgraphics(gcf, 'profiles_memory_vs_math.pdf', ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'none', ...
    'Resolution', 600);
