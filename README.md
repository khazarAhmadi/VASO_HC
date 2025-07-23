This repository contains scripts for the analysis of submillimeter-resolution VASO and BOLD fMRI data in the human hippocampus.

## ⚙️ Dependencies

Before running the scripts, the following software must be installed:

* [ANTs](http://stnava.github.io/ANTs/)
* [AFNI](https://afni.nimh.nih.gov/)
* [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki)
* [LayNii](https://layerfmri.com/laynii/)
* [HippUnfold](https://github.com/khanlab/hippunfold)
* [SPM](https://www.fil.ion.ucl.ac.uk/spm/)
* [ITK-SNAP](http://www.itksnap.org/) *(for manual co-registration)*

Additional code dependencies:

* [`fmri_denoising`](https://github.com/dmascali/fmri_denoising) (for physiological denoising, aka, 'aCompCor')
* [`NORDIC_Raw`](https://github.com/SteenMoeller/NORDIC_Raw) (adapted for VASO)

## 📁 Folder Structure & Main Scripts

* `preproc_seg/`

  * **`Pipeline.sh`**: Main pipeline script.

    * Performs NORDIC denoising (VASO adaptation)
    * Motion and distortion correction
    * Correction for BOLD contamination in VASO
    * Co-registration to T1-weighted image
    * Hippocampal subfield segmentation 
    * Generates study-specific T1 template
    * ⚠️ Memory intensive analysis : ~500 GB RAM recommended*

    * **`nonGM_CompCor.m`**: 
    * Used to regress out physiological noise i.e., non-gray matter signals during GLM analysis

* `voxelwise_glm/`

  * **`subfield_activity_contrast_reg.sh`**:

    * Computes mean signal changes for BOLD and VASO per hippocampal subfield for the contrast memory vs math 
    * Aligns SPM first-level contrast maps with template space
    * Assumes SPM first-level analysis is already completed for all subjects

* `layering/`

  * MATLAB scripts for gray matter sampling in hippocampal subfields
  * Performs layer-wise GLM analyses

* `Figs_script/`

  * MATLAB and R scripts to generate figures from the manuscript

---

## 📚 References
1. Pfaffenrot, V., Bouyeure, A., Gomes, C. A., Kashyap, S., Axmacher, N., & Norris, D. G. (2025). Characterizing BOLD activation patterns in the human hippocampus with laminar fMRI. Imaging Neuroscience, 3, imag_a_00532.

2. DeKraker, J., Haast, R. A., Yousif, M. D., Karat, B., Lau, J. C., Köhler, S., & Khan, A. R. (2022). Automated hippocampal unfolding for morphometry and subfield segmentation with HippUnfold. Elife, 11, e77945.

3. Behzadi, Y., Restom, K., Liau, J., & Liu, T. T. (2007). A component based noise correction method (CompCor) for BOLD and perfusion based fMRI. Neuroimage, 37(1), 90–101.

4. Knudsen, L., Vizioli, L., De Martino, F., Faes, L. K., Handwerker, D. A., Moeller, S., Bandettini, P. A., & Huber, L. (2025). NORDIC denoising on VASO data. Frontiers in Neuroscience, 18, 1499762.
