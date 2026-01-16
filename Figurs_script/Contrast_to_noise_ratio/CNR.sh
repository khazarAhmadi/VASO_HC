#!/bin/bash

# -------------------------------
# CNR calculation for entire hippocampus 
# -------------------------------

# Output text files
output_bold="BOLD_HC.txt"
output_vaso="VASO_HC.txt"

subjects=(001 002 003 004 005 006)


# Loop over subjects
for subj in "${subjects[@]}"; do

    # Path to subject's HippUnfold folder and Hippocampus mask
    hu_index=$(printf "%04d" "$((10#$subj))")  # converts 001 → 0001 and so forth
    HIPPO_MASK="${subj}/anat/HUoutput_T1/hippunfold/sub"$hu_index"/anat/sub-"$hu_index"_hemi-L_space-T1w_desc-subfields_atlas-bigbrain_dseg_bin.nii.gz"
   
    

    # ---------------- BOLD ----------------
    BOLD_DIR="${subj}/func/first_level/bold"
    CON_FILE="${BOLD_DIR}/con_0001.nii"
    RES_FILE="${BOLD_DIR}/ResMS.nii"

    mean_signal=$(fslstats "$CON_FILE" -a -n -k "$HIPPO_MASK" -M)
    mean_resvar=$(fslstats "$RES_FILE" -n -k "$HIPPO_MASK" -M)

    # Compute CNR
    CNR=$(echo "scale=4; $mean_signal / sqrt($mean_resvar)" | bc -l)

    # Save to BOLD text file
    echo "$subj $CNR" >> "$output_bold"

    # ---------------- VASO ----------------
    VASO_DIR="${subj}/func/first_level/vaso"
    CON_FILE_V="${VASO_DIR}/con_0001.nii"
    RES_FILE_V="${VASO_DIR}/ResMS.nii"

    mean_signal_v=$(fslstats "$CON_FILE_V" -a -n -k "$HIPPO_MASK" -M)
    mean_resvar_v=$(fslstats "$RES_FILE_V" -n -k "$HIPPO_MASK" -M)

    # Compute CNR
    CNR_v=$(echo "scale=4; $mean_signal_v / sqrt($mean_resvar_v)" | bc -l)

    # Save to VASO text file
    echo "$subj $CNR_v" >> "$output_vaso"

done

echo "CNR calculation complete!"
echo "Results saved in $output_bold and $output_vaso"

