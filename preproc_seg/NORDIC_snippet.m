%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Snippet to run NORDIC denoising on BOLD and VASO data 

% requires VASO-adaptation of NORDIC to be added to MATLAB path (see
% Knudsen et al., 2025. Frontiers in Neuroscience paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% denoise bold (none-nulled) image
pathtofile = pwd;
ARG.NORDIC= 1;
ARG.noise_volume_last = 0; % consider the last volume to inform about noise distribution
ARG.save_add_info = 1;
ARG.magnitude_only= 1;
ARG.factor_error= 1.2;
ARG.kernel_size_PCA = [10, 10, 10];
fn_magn_in='bold_slice_Vol_Remove.nii.gz';
fn_out=['NORDIC_' fn_magn_in];
fn_phase_in=fn_magn_in;
NIFTI_NORDIC(fn_magn_in,fn_phase_in,fn_out,ARG); 

clearvars -except ARG pathtofile RG

%% Apply the same for VASO
fn_magn_in='vaso_slice_Vol_Remove.nii.gz';
fn_out=['NORDIC_' fn_magn_in];
fn_phase_in=fn_magn_in;
NIFTI_NORDIC(fn_magn_in,fn_phase_in,fn_out,ARG); 
