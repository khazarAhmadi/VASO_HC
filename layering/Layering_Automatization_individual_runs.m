%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Script for sampling hippocampal subfields
%  For each functional run  (BOLD and VASO), it calls VPF_create_hippocampus_layers (must be added to MATLAB path)
% and creates a structure with 30 bins in each hippocampal subregion and each time point 
% requires SPM to be added to MATLAB path.
% We thank Viktor Pfaffenrot for sharing the sampling script with us. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run this script if you want to sample subfileds in each run separately.
%% Locate all functional and anatomical data
fileContent = fileread('FOLDERS_no_PA.txt');
folders = splitlines(fileContent);
folders = folders(1:end-1); % remove last empty entry

fileContent_anat = fileread('FOLDERS-ANAT.txt');
folders_anat = splitlines(fileContent_anat);
folders_anat = folders_anat(~cellfun('isempty', folders_anat)); % remove empty last row

%% prepare the input to be ready to run the layering script 

for m = 1:length(folders)
    bold_data{m,1} = strcat(folders{m,1},'/',...
    'NORDIC_bold_coreg_resampled.nii.gz');
    vaso_data{m,1} = strcat(folders{m,1},'/vaso_split/',...
    'VASO_LN_coreg.nii.gz');
end

rule = '1 2 3 4'; % 1 = subiculum, 2 = ca1, 3 = ca2, 4 = ca3 
N_layers = [20 10]; % 20 bins between inner and outer surfaces and then extende another 10 bins beyond the inner surface to cover SRLM
N_sub = 6; % number of participants

subject_id = 1: N_sub;

subject_bold = cell(length(subject_id), 1);
subject_vaso = cell(length(subject_id), 1);
subjects = unique(cellfun(@(x) x(1:3), bold_data, 'UniformOutput', false)); % extract unique subject IDs

% group functional data of each subject into one cell
for n = 1:N_sub
    subject_bold{n} = bold_data(contains(bold_data, subjects{n}));
    subject_vaso{n} = vaso_data(contains(vaso_data, subjects{n}));
end

layers_bold = cell(N_sub, 1);
layers_vaso = cell(N_sub, 1);

%% run the sampling algorithm 
for n = 1:N_sub
    layers_bold{n} = cell(length(subject_bold{n}), 1);
    layers_vaso{n} = cell(length(subject_vaso{n}), 1);

    % Ensure folders_anat is absolute
    if ~startsWith(folders_anat{n}, '/')
        base_anat_path = fullfile(pwd, folders_anat{n});
    else
        base_anat_path = folders_anat{n};
    end

    surf_path = fullfile(base_anat_path, 'HUoutput_T1', 'hippunfold', ...
                         sprintf('Sub%04d', n), 'surf');

    T1_path = fullfile(base_anat_path, 'HUoutput_T1', 'hippunfold', ...
                       sprintf('Sub%04d', n), 'anat', ...
                       sprintf('sub-%04d_desc-preproc_T1w.nii', n));

    for m = 1:length(subject_bold{n})   

        sampled_img = subject_bold{n}(m);

        layers_bold{n}{m} = VPF_create_hippocampus_layers(sampled_img,subject_id(n),surf_path,...
            T1_path,N_layers,rule);

        clear sampled_img
        sampled_img = subject_vaso{n}(m);
        layers_vaso{n}{m} = VPF_create_hippocampus_layers(sampled_img,subject_id(n),surf_path,...
            T1_path,N_layers,rule);
        clear sampled_img  
    end
end

%% save the output as .mat file
run_idx = 1;  % index to track folders

for n = 1:N_sub
    for m = 1: length(layers_bold{n})
        layers_bold_per_run = layers_bold{n}{m};
        layers_vaso_per_run = layers_vaso{n}{m};
        folder_path = folders{run_idx};
        save(fullfile(folder_path,'layers_BOLD.mat'),'layers_bold_per_run');
        save(fullfile(folder_path,'layers_VASO.mat'),'layers_vaso_per_run');
        run_idx = run_idx + 1;
    end 
end 
