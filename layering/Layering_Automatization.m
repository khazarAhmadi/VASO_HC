%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Script for sampling hippocampal subfields
%  For each functional run  (BOLD and VASO), it calls VPF_create_hippocampus_layers (must be added to MATLAB path)
% and creates a structure with 30 bins in each hippocampal subregion and each time point 
% requires SPM to be added to MATLAB path.
% We thank Viktor Pfaffenrot for sharing the sampling script with us. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Locate merged coregistered functional and anatomical data
fileContent = fileread('unique_no_PA.txt');
folders = splitlines(fileContent);
folders = folders(1:end-1); % remove last empty entry
folders = regexp(folders, '^(.*?/func)/', 'tokens', 'once');
folders = vertcat(folders{:});

fileContent_anat = fileread('FOLDERS-ANAT.txt');
folders_anat = splitlines(fileContent_anat);
folders_anat = folders_anat(~cellfun('isempty', folders_anat)); % remove empty last row
%% prepare the input to be ready to run the layering script 

for m = 1:length(folders)
    bold_data{m,1} = strcat(folders{m,1},'/',...
    'NORDIC_BOLD_coreg_merged_all.nii.gz');
    vaso_data{m,1} = strcat(folders{m,1},'/',...
    'NORDIC_VASO_coreg_merged_all.nii.gz');
end

rule = '1 2 3 4'; % 1 = subiculum, 2 = ca1, 3 = ca2, 4 = ca3 
N_layers = [20 10]; % 20 bins between inner and outer surfaces and then extende another 10 bins beyond the inner surface to cover SRLM
N_sub = 6; % number of participants

subject_id = 1: N_sub;

layers_bold = cell(N_sub, 1);
layers_vaso = cell(N_sub, 1);

%% run the sampling algorithm 
for n = 1:N_sub
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

     sampled_img = {bold_data{n}};

     layers_bold{n} = VPF_create_hippocampus_layers(sampled_img,subject_id(n),surf_path,...
            T1_path,N_layers,rule);

     clear sampled_img
     sampled_img = {vaso_data{n}};
     layers_vaso{n} = VPF_create_hippocampus_layers(sampled_img,subject_id(n),surf_path,...
            T1_path,N_layers,rule);
     clear sampled_img  
end

%% save the output as .mat file
run_idx = 1;  % index to track folders

for n = 1:N_sub
    layers_bold_per_sub = layers_bold{n};
    layers_vaso_per_sub = layers_vaso{n};

    folder_path = folders{run_idx};
    save(fullfile(folder_path,'layers_BOLD_merged.mat'),'layers_bold_per_sub');
    save(fullfile(folder_path,'layers_VASO_merged.mat'),'layers_vaso_per_sub');
    run_idx = run_idx + 1;

end 
% Note: Since the acquisition slab does not cover the right hemisphere, 
% the second row of the output .mat file will contain only zeros. 
% Additionally, the dentate gyrus (DG) and CA4 subfields are not layerized 
% in this analysis, so the second column of the output file will be empty.