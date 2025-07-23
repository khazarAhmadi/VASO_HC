%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Script to extract 10 components for nuisance regressors 
%  (5 for CSF and 5 for WM) from all co-registered BOLD and VASO data
% requires "fmri_compcor.m" to be added to MATLAB path. check https://github.com/dmascali/fmri_denoising
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Locate all functional data

fileContent = fileread('FOLDERS_no_PA.txt');
folders = splitlines(fileContent);
folders = folders(1:end-1);

for m = 1:length(folders)
    data_bold{m,1} = strcat(folders{m,1},'/',...
    'NORDIC_bold_coreg_resampled.nii.gz');
    data_vaso{m,1} = strcat(folders{m,1},'/vaso_split/',...
    'VASO_LN_coreg.nii.gz');
end

%% locate WM and CSF mask in each subject

fileContent_anat = fileread('FOLDERS-ANAT.txt');
folders_anat = splitlines(fileContent_anat);
% remove empty last row
folders_anat = folders_anat(~cellfun('isempty', folders_anat));

for n = 1:length(folders_anat)
    rois{n,1} = strcat(folders_anat{n,1},'/FSL_fast/csf_new_cropped_eroded.nii.gz');
    rois{n,2} = strcat(folders_anat{n,1},'/FSL_fast/WM_new_bin_cropped_eroded.nii.gz');
end

%% Run aCompcor

dime = [5 5];
subjects = unique(cellfun(@(x) x(1:3), data_bold, 'UniformOutput', false)); % extract unique subject IDs
num_subjects = length(subjects);

subject_bold = cell(num_subjects, 1);
subject_vaso = cell(num_subjects, 1);

% group functional data of each subject into one cell
for n = 1:num_subjects
    subject_bold{n} = data_bold(contains(data_bold, subjects{n}));
    subject_vaso{n} = data_vaso(contains(data_vaso, subjects{n}));
end

x_bold = cell(num_subjects, 1);
x_vaso = cell(num_subjects, 1);

for n = 1:num_subjects
    x_bold{n} = cell(length(subject_bold{n}), 1);
    x_vaso{n} = cell(length(subject_vaso{n}), 1);
    for m = 1:length(subject_bold{n})
        x_bold{n}{m} = fmri_compcor(subject_bold{n}{m}, rois{n, :}, dime, 'PolOrder', 1);
        x_vaso{n}{m} = fmri_compcor(subject_vaso{n}{m}, rois{n, :}, dime, 'PolOrder', 1);
    end
end

%% save the output as txt file 

run_idx = 1;  % index to track folders

for n = 1:num_subjects
    for m = 1:length(x_bold{n})
        % Get current output
        comp_bold = x_bold{n}{m}; 
        comp_vaso = x_vaso{n}{m}; 

        % Get corresponding folder path
        folder_path = folders{run_idx};

        % Define output file path
        out_path_bold = fullfile(folder_path, 'compcor_bold.txt');
        out_path_vaso = fullfile(folder_path, 'compcor_vaso.txt');

        % Save as plain text
        writematrix(comp_bold, out_path_bold,'Delimiter','\t');
        writematrix(comp_vaso, out_path_vaso,'Delimiter','\t');
        run_idx = run_idx + 1;
    end
end
