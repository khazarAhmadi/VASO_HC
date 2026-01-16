%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute and plot laminar tSNR for BOLD and VASO across hippocampal subfields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% define path and load layers.mat files

folders = readlines('folders.txt');
folders = strtrim(folders);
folders(folders == "") = [];   % remove empty lines

nSub = numel(folders);
LAYERS = {};   % will grow dynamically

for m = 1:nSub

    baseDir = folders(m);

    % Look for session_* subfolders
    sessDirs = dir(fullfile(baseDir, 'session_*'));
    sessDirs = sessDirs([sessDirs.isdir]);

    if isempty(sessDirs)
        % -------- single-session case --------
        fname = fullfile(baseDir, 'layers_bold.mat');

        if ~exist(fname,'file')
            warning('Missing layers_bold.mat in %s', baseDir);
            continue
        end

        tmp = load(fname);   % loads variable "layers"
        LAYERS{m,1} = tmp.layers;

    else
        % -------- multi-session case --------
        for n = 1:numel(sessDirs)
            sessPath = fullfile(baseDir, sessDirs(n).name);
            fname = fullfile(sessPath, 'layers_bold.mat');

            if ~exist(fname,'file')
                warning('Missing layers_bold.mat in %s', sessPath);
                continue
            end

            tmp = load(fname);
            LAYERS{m,n} = tmp.layers;
        end
    end
end
%% concatenate the layers, split them in each subfields 
[nSub, nSess] = size(LAYERS); %nSess = maximum number of sessions which is 3.
allLayers = struct('sub_b', [], 'ca1_b', [], 'ca2_b', [], 'ca3_b', []);

for m = 1:nSub
    for n = 1:nSess

        if isempty(LAYERS{m,n})
            continue
        end

        L = LAYERS{m,n}{1,1};

        allLayers.sub_b = [allLayers.sub_b, L(:,:,1)];
        allLayers.ca1_b = [allLayers.ca1_b, L(:,:,2)];
        allLayers.ca2_b = [allLayers.ca2_b, L(:,:,3)];
        allLayers.ca3_b = [allLayers.ca3_b, L(:,:,4)];

    end
end
%% get the tSNR for every run of every subject (BOLD)
N_vol = 150; % number of volumes per run
N_run = length(allLayers.sub_b)./N_vol; % total acquired runs in all participants
N_layers = 30;

for n = 1:N_run
    cols = (n-1)*N_vol + 1 : n*N_vol;  
    sub_b{n} = allLayers.sub_b(:, cols);
    ca1_b{n} = allLayers.ca1_b(:, cols);
    ca2_b{n} = allLayers.ca2_b(:, cols);
    ca3_b{n} = allLayers.ca3_b(:, cols);
end


for n = 1:N_run
    for i = 1:N_layers
        mean_sub_b(n,i) = mean(sub_b{1,n}(i,:));
        sd_sub_b(n,i) = std(sub_b{1,n}(i,:));
        tsnr_sub_b(n,i) = mean_sub_b(n,i)./sd_sub_b(n,i);
        mean_ca1_b(n,i) = mean(ca1_b{1,n}(i,:));
        sd_ca1_b(n,i) = std(ca1_b{1,n}(i,:));
        tsnr_ca1_b(n,i) = mean_ca1_b(n,i)./sd_ca1_b(n,i);
        mean_ca2_b(n,i) = mean(ca2_b{1,n}(i,:));
        sd_ca2_b(n,i) = std(ca2_b{1,n}(i,:));
        tsnr_ca2_b(n,i) = mean_ca2_b(n,i)./sd_ca2_b(n,i);
        mean_ca3_b(n,i) = mean(ca3_b{1,n}(i,:));
        sd_ca3_b(n,i) = std(ca3_b{1,n}(i,:));
        tsnr_ca3_b(n,i) = mean_ca3_b(n,i)./sd_ca3_b(n,i);
    end
end 

%% caluclate the mean tsnr across runs 

avg_tsnr_sub_b = mean(tsnr_sub_b);
sem_tsnr_sub_b = std(tsnr_sub_b)./sqrt(N_run); 
upper_sub_b = avg_tsnr_sub_b + sem_tsnr_sub_b;
lower_sub_b = avg_tsnr_sub_b - sem_tsnr_sub_b;

avg_tsnr_ca1_b = mean(tsnr_ca1_b);
sem_tsnr_ca1_b = std(tsnr_ca1_b)./sqrt(N_run); 
upper_ca1_b = avg_tsnr_ca1_b + sem_tsnr_ca1_b;
lower_ca1_b = avg_tsnr_ca1_b - sem_tsnr_ca1_b;

avg_tsnr_ca2_b = mean(tsnr_ca2_b);
sem_tsnr_ca2_b = std(tsnr_ca2_b)./sqrt(N_run); 
upper_ca2_b = avg_tsnr_ca2_b + sem_tsnr_ca2_b;
lower_ca2_b = avg_tsnr_ca2_b - sem_tsnr_ca2_b;

avg_tsnr_ca3_b = mean(tsnr_ca3_b);
sem_tsnr_ca3_b = std(tsnr_ca3_b)./sqrt(N_run); 
upper_ca3_b = avg_tsnr_ca3_b + sem_tsnr_ca3_b;
lower_ca3_b = avg_tsnr_ca3_b - sem_tsnr_ca3_b;

%% Repeat the same for VASO 
LAYERS_V = {};   % will grow dynamically

for m = 1:nSub

    baseDir = folders(m);

    % Look for session_* subfolders
    sessDirs = dir(fullfile(baseDir, 'session_*'));
    sessDirs = sessDirs([sessDirs.isdir]);

    if isempty(sessDirs)
        % -------- single-session case --------
        fname = fullfile(baseDir, 'layers_vaso.mat');

        if ~exist(fname,'file')
            warning('Missing layers_vaso.mat in %s', baseDir);
            continue
        end

        tmp = load(fname);   % loads variable "LAYERS"
        LAYERS_V{m,1} = tmp.layers;

        else
        % -------- multi-session case --------
        for n = 1:numel(sessDirs)
            sessPath = fullfile(baseDir, sessDirs(n).name);
            fname = fullfile(sessPath, 'layers_vaso.mat');

            if ~exist(fname,'file')
                warning('Missing layers_vaso.mat in %s', sessPath);
                continue
            end

            tmp = load(fname);
            LAYERS_V{m,n} = tmp.layers;
        end
    end
end

%% concatenate the layers, split them in each subfields (for VASO)
[nSub, nSess] = size(LAYERS_V);
allLayers_v = struct('sub_v', [], 'ca1_v', [], 'ca2_v', [], 'ca3_v', []);

for m = 1:nSub
    for n = 1:nSess

        if isempty(LAYERS_V{m,n})
            continue
        end

        L = LAYERS_V{m,n}{1,1};

        allLayers_v.sub_v = [allLayers_v.sub_v, L(:,:,1)];
        allLayers_v.ca1_v = [allLayers_v.ca1_v, L(:,:,2)];
        allLayers_v.ca2_v = [allLayers_v.ca2_v, L(:,:,3)];
        allLayers_v.ca3_v = [allLayers_v.ca3_v, L(:,:,4)];

    end
end

%% get the tSNR for every run of every subject (VASO)

for n = 1:N_run
    cols = (n-1)*N_vol + 1 : n*N_vol;  
    sub_v{n} = allLayers_v.sub_v(:, cols);
    ca1_v{n} = allLayers_v.ca1_v(:, cols);
    ca2_v{n} = allLayers_v.ca2_v(:, cols);
    ca3_v{n} = allLayers_v.ca3_v(:, cols);
end


for n = 1:N_run
    for i = 1:N_layers
        mean_sub_v(n,i) = mean(sub_v{1,n}(i,:));
        sd_sub_v(n,i) = std(sub_v{1,n}(i,:));
        tsnr_sub_v(n,i) = mean_sub_v(n,i)./sd_sub_v(n,i);
        mean_ca1_v(n,i) = mean(ca1_v{1,n}(i,:));
        sd_ca1_v(n,i) = std(ca1_v{1,n}(i,:));
        tsnr_ca1_v(n,i) = mean_ca1_v(n,i)./sd_ca1_v(n,i);
        mean_ca2_v(n,i) = mean(ca2_v{1,n}(i,:));
        sd_ca2_v(n,i) = std(ca2_v{1,n}(i,:));
        tsnr_ca2_v(n,i) = mean_ca2_v(n,i)./sd_ca2_v(n,i);
        mean_ca3_v(n,i) = mean(ca3_v{1,n}(i,:));
        sd_ca3_v(n,i) = std(ca3_v{1,n}(i,:));
        tsnr_ca3_v(n,i) = mean_ca3_v(n,i)./sd_ca3_v(n,i);
    end
end 

%% caluclate the mean tsnr across runs (VASO)

avg_tsnr_sub_v = mean(tsnr_sub_v);
sem_tsnr_sub_v = std(tsnr_sub_v)./sqrt(N_run); 
upper_sub_v = avg_tsnr_sub_v + sem_tsnr_sub_v;
lower_sub_v = avg_tsnr_sub_v - sem_tsnr_sub_v;

avg_tsnr_ca1_v = mean(tsnr_ca1_v);
sem_tsnr_ca1_v = std(tsnr_ca1_v)./sqrt(N_run); 
upper_ca1_v = avg_tsnr_ca1_v + sem_tsnr_ca1_v;
lower_ca1_v = avg_tsnr_ca1_v - sem_tsnr_ca1_v;

avg_tsnr_ca2_v = mean(tsnr_ca2_v);
sem_tsnr_ca2_v = std(tsnr_ca2_v)./sqrt(N_run); 
upper_ca2_v = avg_tsnr_ca2_v + sem_tsnr_ca2_v;
lower_ca2_v = avg_tsnr_ca2_v - sem_tsnr_ca2_v;

avg_tsnr_ca3_v = mean(tsnr_ca3_v);
sem_tsnr_ca3_v = std(tsnr_ca3_v)./sqrt(N_run); 
upper_ca3_v = avg_tsnr_ca3_v + sem_tsnr_ca3_v;
lower_ca3_v = avg_tsnr_ca3_v - sem_tsnr_ca3_v;

%% Create the plots 

depths = (1:N_layers);
names = {'SRLM','inner','outer'};
x_SRLM(1,1:30) = 10; % 10th bin is where SRLM ends
x_inner(1,1:30) = 20; % 20th bin where inner surface ends
y_threshold = linspace(0,120,30);

figure('Color','w');
set(gcf, 'GraphicsSmoothing', 'off');
set(gcf, 'Renderer', 'painters');  
subplot(2,2,1)
xlim([1 30]);hold on
plot(x_SRLM,y_threshold,'LineWidth',2,'Color','k','LineStyle','--');hold on
set(gca,'xtick',[5,15,25],'xticklabel',names,'Color','w')
plot(x_inner,y_threshold,'LineWidth',2,'Color','k','LineStyle','--');

ylabel('tSNR [a.u]');
box on
ax = gca;
ax.LineWidth = 2;
ax.FontSize = 15; 
fill([depths, fliplr(depths)], [upper_sub_b, fliplr(lower_sub_b)], [0 0 0.8], ...
     'FaceAlpha', 0.2, 'EdgeColor', 'none');
H1 = plot(depths,avg_tsnr_sub_b,'LineWidth',3,'Color','b'); 
fill([depths, fliplr(depths)], [upper_sub_v, fliplr(lower_sub_v)], [0.8 0 0], ...
     'FaceAlpha', 0.2, 'EdgeColor', 'none');
H2 = plot(depths,avg_tsnr_sub_v,'LineWidth',3,'Color','r'); 

title('Sub'); 
hleglines = [H1(1),H2(1)];
Hleg = legend(hleglines,'BOLD','VASO','Location','southwest');hold on

subplot(2,2,2)
xlim([1 30]);hold on
plot(x_SRLM,y_threshold,'LineWidth',2,'Color','k','LineStyle','--');hold on
set(gca,'xtick',[5,15,25],'xticklabel',names,'Color','w')
plot(x_inner,y_threshold,'LineWidth',2,'Color','k','LineStyle','--');

ylabel('tSNR [a.u]');
box on
ax = gca;
ax.LineWidth = 2;
ax.FontSize = 15; 
fill([depths, fliplr(depths)], [upper_ca1_b, fliplr(lower_ca1_b)], [0 0 0.8], ...
     'FaceAlpha', 0.2, 'EdgeColor', 'none');
H3 = plot(depths,avg_tsnr_ca1_b,'LineWidth',3,'Color','b'); 
fill([depths, fliplr(depths)], [upper_ca1_v, fliplr(lower_ca1_v)], [0.8 0 0], ...
     'FaceAlpha', 0.2, 'EdgeColor', 'none');
H4 = plot(depths,avg_tsnr_ca1_v,'LineWidth',3,'Color','r'); 

title('CA1'); 
hleglines = [H3(1),H4(1)];
Hleg = legend(hleglines,'BOLD','VASO','Location','southwest');hold on

subplot(2,2,3)
xlim([1 30]);hold on
plot(x_SRLM,y_threshold,'LineWidth',2,'Color','k','LineStyle','--');hold on
set(gca,'xtick',[5,15,25],'xticklabel',names,'Color','w')
plot(x_inner,y_threshold,'LineWidth',2,'Color','k','LineStyle','--');

ylabel('tSNR [a.u]');
box on
ax = gca;
ax.LineWidth = 2;
ax.FontSize = 15; 
fill([depths, fliplr(depths)], [upper_ca2_b, fliplr(lower_ca2_b)], [0 0 0.8], ...
     'FaceAlpha', 0.2, 'EdgeColor', 'none');
H5 = plot(depths,avg_tsnr_ca2_b,'LineWidth',3,'Color','b'); 
fill([depths, fliplr(depths)], [upper_ca2_v, fliplr(lower_ca2_v)], [0.8 0 0], ...
     'FaceAlpha', 0.2, 'EdgeColor', 'none');
H6 = plot(depths,avg_tsnr_ca2_v,'LineWidth',3,'Color','r'); 

title('CA2'); 
hleglines = [H5(1),H6(1)];
Hleg = legend(hleglines,'BOLD','VASO','Location','southwest');hold on

subplot(2,2,4)
xlim([1 30]);hold on
plot(x_SRLM,y_threshold,'LineWidth',2,'Color','k','LineStyle','--');hold on
set(gca,'xtick',[5,15,25],'xticklabel',names,'Color','w')
plot(x_inner,y_threshold,'LineWidth',2,'Color','k','LineStyle','--');

ylabel('tSNR [a.u]');
box on
ax = gca;
ax.LineWidth = 2;
ax.FontSize = 15; 
fill([depths(11:30), fliplr(depths(11:30))], [upper_ca3_b(11:30), fliplr(lower_ca3_b(11:30))], [0 0 0.8], ...
     'FaceAlpha', 0.2, 'EdgeColor', 'none');
H7 = plot(depths(11:30),avg_tsnr_ca3_b(11:30),'LineWidth',3,'Color','b'); 
fill([depths(11:30), fliplr(depths(11:30))], [upper_ca3_v(11:30), fliplr(lower_ca3_v(11:30))], [0.8 0 0], ...
     'FaceAlpha', 0.2, 'EdgeColor', 'none');
H8 = plot(depths(11:30),avg_tsnr_ca3_v(11:30),'LineWidth',3,'Color','r'); 

title('CA3'); 
hleglines = [H7(1),H8(1)];
Hleg = legend(hleglines,'BOLD','VASO','Location','southwest');hold on

set(gcf, 'Renderer', 'painters');
print(gcf, '-depsc2', 'Laminar_tSNR_BOLD_VASO.eps'); % save as eps

