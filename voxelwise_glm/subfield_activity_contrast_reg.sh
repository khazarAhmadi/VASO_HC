#!/bin/bash 

# This script should be run from the parent directory (i.e., the "Part2" folder) after the SPM first-level analyses have been completed.

#--------------------------------------------------
# Section 1: Compute the mean activity within each subfield and subject, based on a conjunction ROI of significant clusters for BOLD + VASO from the memory vs. math contrast.
#--------------------------------------------------

# first get a binary mask of each subfield from HippUnfold outputs. 
# Subfield label IDs based on HippUnfold's BigBrain atlas i.e., 1 = Subiculum, 2 = CA1, 3 = CA2, 4 = CA3, 5 = CA4, 6 = DG, 7 = SRLM  
# Note DG and CA4 are merged into one ROI.

while IFS= read -r line1 && IFS= read -r line2 <&3; do 
fslmaths $line1/HUoutput_T1/hippunfold/Sub$line2/anat/sub-"$line2"_hemi-L_space-T1w_desc-subfields_atlas-bigbrain_dseg.nii.gz -uthr 1 -thr 1 $line1/HUoutput_T1/hippunfold/left-sub.nii.gz;
fslmaths $line1/HUoutput_T1/hippunfold/Sub$line2/anat/sub-"$line2"_hemi-L_space-T1w_desc-subfields_atlas-bigbrain_dseg.nii.gz -uthr 2 -thr 2 -bin $line1/HUoutput_T1/hippunfold/left-ca1.nii.gz;
fslmaths $line1/HUoutput_T1/hippunfold/Sub$line2/anat/sub-"$line2"_hemi-L_space-T1w_desc-subfields_atlas-bigbrain_dseg.nii.gz -uthr 3 -thr 3 -bin $line1/HUoutput_T1/hippunfold/left-ca2.nii.gz;
fslmaths $line1/HUoutput_T1/hippunfold/Sub$line2/anat/sub-"$line2"_hemi-L_space-T1w_desc-subfields_atlas-bigbrain_dseg.nii.gz -uthr 4 -thr 4 -bin $line1/HUoutput_T1/hippunfold/left-ca3.nii.gz;
fslmaths $line1/HUoutput_T1/hippunfold/Sub$line2/anat/sub-"$line2"_hemi-L_space-T1w_desc-subfields_atlas-bigbrain_dseg.nii.gz -uthr 6 -thr 5 -bin $line1/HUoutput_T1/hippunfold/left-DG.nii.gz;
fslmaths $line1/HUoutput_T1/hippunfold/Sub$line2/anat/sub-"$line2"_hemi-L_space-T1w_desc-subfields_atlas-bigbrain_dseg.nii.gz -uthr 7 -thr 7 -bin $line1/HUoutput_T1/hippunfold/left-SRLM.nii.gz;
done <  FOLDERS-ANAT.txt 3< Hippunfold_indices.txt


subjects=(001 002 003 004 005 006)

# bold
paste -d' ' FOLDERS-ANAT.txt <(printf "%s\n" "${subjects[@]}") | while read anat_dir subj; do
  subj_dir_b="${subj}/func/first_level/bold"
  fslmaths "${subj_dir_b}/mem_math_sig_clust.nii" -mul "${anat_dir}/HUoutput_T1/hippunfold/left-sub.nii.gz" -nan "${subj_dir_b}/Sub_sig_clust.nii.gz"
  fslmaths "${subj_dir_b}/mem_math_sig_clust.nii" -mul "${anat_dir}/HUoutput_T1/hippunfold/left-ca1.nii.gz" -nan "${subj_dir_b}/ca1_sig_clust.nii.gz"
  fslmaths "${subj_dir_b}/mem_math_sig_clust.nii" -mul "${anat_dir}/HUoutput_T1/hippunfold/left-ca2.nii.gz" -nan "${subj_dir_b}/ca2_sig_clust.nii.gz"
  fslmaths "${subj_dir_b}/mem_math_sig_clust.nii" -mul "${anat_dir}/HUoutput_T1/hippunfold/left-ca3.nii.gz" -nan "${subj_dir_b}/ca3_sig_clust.nii.gz"
  fslmaths "${subj_dir_b}/mem_math_sig_clust.nii" -mul "${anat_dir}/HUoutput_T1/hippunfold/left-DG.nii.gz" -nan "${subj_dir_b}/DG_sig_clust.nii.gz"
  fslmaths "${subj_dir_b}/mem_math_sig_clust.nii" -mul "${anat_dir}/HUoutput_T1/hippunfold/left-SRLM.nii.gz" -nan "${subj_dir_b}/SRLM_sig_clust.nii.gz"
done

# repeat it for vaso
paste -d' ' FOLDERS-ANAT.txt <(printf "%s\n" "${subjects[@]}") | while read anat_dir subj; do
  subj_dir_v="${subj}/func/first_level/vaso"
  fslmaths "${subj_dir_v}/mem_math_sig_clust.nii" -mul "${anat_dir}/HUoutput_T1/hippunfold/left-sub.nii.gz" -nan "${subj_dir_v}/Sub_sig_clust.nii.gz"
  fslmaths "${subj_dir_v}/mem_math_sig_clust.nii" -mul "${anat_dir}/HUoutput_T1/hippunfold/left-ca1.nii.gz" -nan "${subj_dir_v}/ca1_sig_clust.nii.gz"
  fslmaths "${subj_dir_v}/mem_math_sig_clust.nii" -mul "${anat_dir}/HUoutput_T1/hippunfold/left-ca2.nii.gz" -nan "${subj_dir_v}/ca2_sig_clust.nii.gz"
  fslmaths "${subj_dir_v}/mem_math_sig_clust.nii" -mul "${anat_dir}/HUoutput_T1/hippunfold/left-ca3.nii.gz" -nan "${subj_dir_v}/ca3_sig_clust.nii.gz"
  fslmaths "${subj_dir_v}/mem_math_sig_clust.nii" -mul "${anat_dir}/HUoutput_T1/hippunfold/left-DG.nii.gz" -nan "${subj_dir_v}/DG_sig_clust.nii.gz"
  fslmaths "${subj_dir_v}/mem_math_sig_clust.nii" -mul "${anat_dir}/HUoutput_T1/hippunfold/left-SRLM.nii.gz" -nan "${subj_dir_v}/SRLM_sig_clust.nii.gz"
done


# concatenate the significant subfields  

for subj in "${subjects[@]}"; do
	subj_dir_b="${subj}/func/first_level/bold"
	subj_dir_v="${subj}/func/first_level/vaso"
	fslmaths "${subj_dir_b}/Sub_sig_clust.nii.gz" -bin -add "${subj_dir_v}/Sub_sig_clust.nii.gz" -bin "${subj}/func/first_level/Sub_sig_both_contrasts.nii.gz"
	fslmaths "${subj_dir_b}/ca1_sig_clust.nii.gz" -bin -add "${subj_dir_v}/ca1_sig_clust.nii.gz" -bin "${subj}/func/first_level/CA1_sig_both_contrasts.nii.gz"
	fslmaths "${subj_dir_b}/ca2_sig_clust.nii.gz" -bin -add "${subj_dir_v}/ca2_sig_clust.nii.gz" -bin "${subj}/func/first_level/CA2_sig_both_contrasts.nii.gz"
	fslmaths "${subj_dir_b}/ca3_sig_clust.nii.gz" -bin -add "${subj_dir_v}/ca3_sig_clust.nii.gz" -bin "${subj}/func/first_level/CA3_sig_both_contrasts.nii.gz"
	fslmaths "${subj_dir_b}/DG_sig_clust.nii.gz" -bin -add "${subj_dir_v}/DG_sig_clust.nii.gz" -bin "${subj}/func/first_level/DG_sig_both_contrasts.nii.gz"
	fslmaths "${subj_dir_b}/SRLM_sig_clust.nii.gz" -bin -add "${subj_dir_v}/SRLM_sig_clust.nii.gz" -bin "${subj}/func/first_level/SRLM_sig_both_contrasts.nii.gz"
done 


# Define the subfields and their corresponding output names

subfields=("Sub" "ca1" "ca2" "ca3" "DG" "SRLM")
output_files_bold=("bold_Sub_activity.txt" "bold_CA1_activity.txt" "bold_CA2_activity.txt" "bold_CA3_activity.txt" "bold_DG_activity.txt" "bold_SRLM_activity.txt")
output_files_vaso=("vaso_Sub_activity.txt" "vaso_CA1_activity.txt" "vaso_CA2_activity.txt" "vaso_CA3_activity.txt" "vaso_DG_activity.txt" "vaso_SRLM_activity.txt")

# Loop over subjects and subfields
for i in "${!subjects[@]}"; do
  subj="${subjects[$i]}"
  subj_dir="${subj}/func/first_level"

  # Loop over the subfields and compute fslstats for each one
  for j in "${!subfields[@]}"; do
    subfield="${subfields[$j]}"
    output_file_bold="${output_files_bold[$j]}"
    output_file_vaso="${output_files_vaso[$j]}"

    # Compute the mean activity for each subfield and append to the corresponding output file
    fslstats "${subj_dir}/bold/con_0001.nii" -k  "${subj_dir}/${subfield}_sig_both_contrasts.nii.gz" -M -S >> "${output_file_bold}"
    fslstats "${subj_dir}/vaso/con_0001.nii" -k "${subj_dir}/${subfield}_sig_both_contrasts.nii.gz" -M -S  >> "${output_file_vaso}"
  done
done


# Now get the mean activity across entire anatomically defined ROIs regardless of significant clusters 

anat_masks=("left-sub" "left-ca1" "left-ca2" "left-ca3" "left-DG" "left-SRLM")

output_files_bold_anat=(
  "bold_Sub_anat_activity.txt"
  "bold_CA1_anat_activity.txt"
  "bold_CA2_anat_activity.txt"
  "bold_CA3_anat_activity.txt"
  "bold_DG_anat_activity.txt"
  "bold_SRLM_anat_activity.txt"
)

output_files_vaso_anat=(
  "vaso_Sub_anat_activity.txt"
  "vaso_CA1_anat_activity.txt"
  "vaso_CA2_anat_activity.txt"
  "vaso_CA3_anat_activity.txt"
  "vaso_DG_anat_activity.txt"
  "vaso_SRLM_anat_activity.txt"
)


for i in "${!subjects[@]}"; do
  subj="${subjects[$i]}"
  subj_dir="${subj}/func/first_level"

  # you need the correct anat_dir for each subject
  anat_dir=$(sed -n "$((i+1))p" FOLDERS-ANAT.txt)

  for j in "${!anat_masks[@]}"; do
    mask="${anat_dir}/HUoutput_T1/hippunfold/${anat_masks[$j]}.nii.gz"

    # BOLD
    fslstats "${subj_dir}/bold/con_0001.nii" \
      -k "$mask" -M -S >> "${output_files_bold_anat[$j]}"

    # VASO
    fslstats "${subj_dir}/vaso/con_0001.nii" \
      -k "$mask" -M -S >> "${output_files_vaso_anat[$j]}"

  done
done

#--------------------------------------------------
# Section 2: Register individual contrast maps to template space, required for second-level analysis
#--------------------------------------------------



affines=(
  "T1_template/SST_T1_sub00100GenericAffine.mat"
  "T1_template/SST_T1_sub00210GenericAffine.mat"
  "T1_template/SST_T1_sub00320GenericAffine.mat"
  "T1_template/SST_T1_sub00430GenericAffine.mat"
  "T1_template/SST_T1_sub00540GenericAffine.mat"
  "T1_template/SST_T1_sub00650GenericAffine.mat"
)

warps=(
  "T1_template/SST_T1_sub00101Warp.nii.gz"
  "T1_template/SST_T1_sub00211Warp.nii.gz"
  "T1_template/SST_T1_sub00321Warp.nii.gz"
  "T1_template/SST_T1_sub00431Warp.nii.gz"
  "T1_template/SST_T1_sub00541Warp.nii.gz"
  "T1_template/SST_T1_sub00651Warp.nii.gz"
)

for i in "${!subjects[@]}"; do
  subj="${subjects[$i]}"
  subj_dir="${subj}/func/first_level"
  affine="${affines[$i]}"
  warp="${warps[$i]}"

  for contrast in {1..4}; do
    con_file=$(printf "con_000%d.nii" "$contrast")
    out_file=$(printf "con_000%d_reg_Temp.nii" "$contrast")

    for modality in bold vaso; do
      antsApplyTransforms -d 3 \
        -i "${subj_dir}/${modality}/${con_file}" \
        -r T1_template/SST_template0.nii.gz \
        -o "${subj_dir}/${modality}/${out_file}" \
        -t "${affine}" \
        -t "${warp}"
    done
  done
done


#--------------------------------------------------
# Section 3: Crop template GM mask to be used later during cluster-level inference
#--------------------------------------------------

cd T1_template/

# to get an estimate of functional slab covergae, take mean vaso image of a random subjec e.g. 002 and align it with template space

aantsApplyTransforms -d 3 -i ../002/func/Run2/QC_vaso/vaso_mean.nii.gz -r SST_template0.nii.gz -o ../002/func/Run2/QC_vaso/vaso_mean_coreg_template.nii.gz -t SST_T1_sub00210GenericAffine.mat -t SST_T1_sub00211Warp.nii.gz 
 
fslmaths FSL_fast/GM_new_bin.nii.gz -mas ../002/func/Run2/QC_vaso/vaso_mean_coreg_template.nii.gz FSL_fast/GM_new_bin_cropped.nii.gz 
