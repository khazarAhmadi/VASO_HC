#!/bin/bash 

# This script should be run from the parent directory i.e., Part2 folder

#--------------------------------------------------
# Section1: entire preprocessing pipeline which consists of the following steps:
# 1) Remove slices at the outermost edge of the acquisition slab
# 2) Remove dummy volumes
# 3) Apply NORDIC denoising
# 4) Perform motion and distortion correction (if reverse PE data is available)
# 5) MPRAGE-ize the anatomy (requires `presurfer`)
# 6) Co-register functional data to anatomical reference
# 7) Generate basic QC measures (e.g., mean images and tSNR maps)
#--------------------------------------------------


# Remove 3 slices (1 from bottom and 2 from top at the edge of acquisition)
while IFS= read -r line1 && IFS= read -r line2 <&3; do  fslroi $line1 $line2 0 -1 0 -1 1 33; done <  Files.txt 3< Files_Slice_rm.txt

# Discard two dummy volumes and the last two noise volumes at the end 
cat FOLDERS.txt | while read line; do NumVol=`3dinfo -nv "$line/vaso_sliceRemove.nii.gz"`; b="$((NumVol - 3))";  3dTcat -prefix $line/vaso_slice_Vol_Remove.nii.gz $line/vaso_sliceRemove.nii.gz[2.."$b"]; done
cat FOLDERS.txt | while read line; do NumVol=`3dinfo -nv "$line/bold_sliceRemove.nii.gz"`; b="$((NumVol - 3))";  3dTcat -prefix $line/bold_slice_Vol_Remove.nii.gz $line/bold_sliceRemove.nii.gz[2.."$b"]; done

# Run NORDIC, and gzip the outputs 
START_DIR=$(pwd)
cat FOLDERS.txt | while read line; do cd $line/ || continue ; matlab -nodisplay -nodesktop -r "run('$START_DIR/NORDIC_snippet.m'); quit"; cd "$START_DIR"; done  
cat FOLDERS.txt | while read line; do cd $line/ || continue ; gzip NORDIC_bold_slice_Vol_Remove.nii; gzip NORDIC_vaso_slice_Vol_Remove.nii;  cd "$START_DIR"; done 


# Do motion and distortion correction using ANTS pipeline 

cat FOLDERS_no_PA.txt | while read line; do cd $line/ || continue 
if [ -d ../reverse_phase ]; then
  "$START_DIR"/sk_ants_Realign_Estimate_KA.sh -n 28 -t 3 -a NORDIC_bold_slice_Vol_Remove.nii.gz -b ../reverse_phase/NORDIC_bold_slice_Vol_Remove.nii.gz
  "$START_DIR"/sk_ants_Realign_Estimate_KA.sh -n 28 -t 3 -a NORDIC_vaso_slice_Vol_Remove.nii.gz -b ../reverse_phase/NORDIC_vaso_slice_Vol_Remove.nii.gz
else
   echo "Skipping $line: reverse_phase directory not found."
fi
cd "$START_DIR"; done


# merge motion and distortion corrected data 
cat FOLDERS_no_PA.txt | while read line; do cd $line/ || continue 
if [ -d ../reverse_phase ]; then
   cd NORDIC_bold_slice_Vol_Remove_mats
   M=$(ls *_0GenericAffine.mat | wc -l)
   F=$(($M-1))
   D=$((1$F))
   for n in $(seq 1000 $D); do
    antsApplyTransforms \
      --output-data-type int \
      --dimensionality 3 \
      --interpolation LanczosWindowedSinc \
      --transform NORDIC_bold_slice_Vol_Remove_"$n"_0GenericAffine.mat \  
      --transform ../NORDIC_bold_slice_Vol_Remove_DistCorr_01Warp.nii.gz \
      --input ../NORDIC_bold_slice_Vol_Remove_split/NORDIC_bold_slice_Vol_Remove_"$n".nii.gz \
      --output ../NORDIC_bold_slice_Vol_Remove_split/NORDIC_bold_slice_Vol_Remove_"$n"_MoCo_DisCor.nii.gz \
      --reference-image ../NORDIC_bold_slice_Vol_Remove_DistCorr_template0.nii.gz \
      -v
  done
  cd ../NORDIC_bold_slice_Vol_Remove_split
  ImageMath 4 NORDIC_bold_MoCo_DisCor_merged.nii.gz TimeSeriesAssemble 3 0 *_MoCo_DisCor.nii.gz
  
  # do the same for vaso
  cd ../NORDIC_vaso_slice_Vol_Remove_mats
  for n in $(seq 1000 $D); do
    antsApplyTransforms \
      --output-data-type int \
      --dimensionality 3 \
      --interpolation LanczosWindowedSinc \
      --transform NORDIC_vaso_slice_Vol_Remove_"$n"_0GenericAffine.mat \ 
      --transform ../NORDIC_vaso_slice_Vol_Remove_DistCorr_01Warp.nii.gz \
      --input ../NORDIC_vaso_slice_Vol_Remove_split/NORDIC_vaso_slice_Vol_Remove_"$n".nii.gz \
      --output ../NORDIC_vaso_slice_Vol_Remove_split/NORDIC_vaso_slice_Vol_Remove_"$n"_MoCo_DisCor.nii.gz \
      --reference-image ../NORDIC_vaso_slice_Vol_Remove_DistCorr_template0.nii.gz \
      -v
  done
  cd ../NORDIC_vaso_slice_Vol_Remove_split
  ImageMath 4 NORDIC_vaso_MoCo_DisCor_merged.nii.gz TimeSeriesAssemble 3 0 *_MoCo_DisCor.nii.gz
  else
  echo "Skipping $line: reverse_phase directory not found."
fi
cd "$START_DIR"; done  



# In Sub 002 and session_3 in sub 001 where there are no reverse-phase data, adapted version of k_ants_Realign_Estimate_KA.sh is performed for motion correction

cd 002/func/Run2/ 
mkdir NORDIC_bold_slice_Vol_Remove_split NORDIC_vaso_slice_Vol_Remove_split

ImageMath 4 NORDIC_bold_slice_Vol_Remove_split/NORDIC_bold_slice_Vol_Remove_.nii.gz TimeSeriesDisassemble NORDIC_bold_slice_Vol_Remove.nii.gz

ImageMath 4 NORDIC_vaso_slice_Vol_Remove_split/NORDIC_vaso_slice_Vol_Remove_.nii.gz TimeSeriesDisassemble NORDIC_vaso_slice_Vol_Remove.nii.gz

cd NORDIC_bold_slice_Vol_Remove_split/;

# 1000-1149 indicate volume indices, counting from 0, making 150 volumes in total.
VOLUMES=$(seq 1000 1149)
for n in $VOLUMES; do 
  antsRegistration \
    --verbose 1 \
    --float 1 \
    --dimensionality 3 \
    --use-histogram-matching 1 \
    --interpolation LanczosWindowedSinc \
    --collapse-output-transforms 1 \
    --output [ NORDIC_bold_$n, NORDIC_bold_"$n"_slice_Vol_Remove_Warped.nii.gz , 1 ] \
    --winsorize-image-intensities [ 0.005 , 0.995 ] \
    --initial-moving-transform [ NORDIC_bold_slice_Vol_Remove_1000.nii.gz , NORDIC_bold_slice_Vol_Remove_$n.nii.gz , 1 ] \
    --transform Rigid[0.1] \
    --metric MI[ NORDIC_bold_slice_Vol_Remove_1000.nii.gz , NORDIC_bold_slice_Vol_Remove_$n.nii.gz , 1 , 64 , Regular , 0.25] \
    --convergence [ 500x250 , 1e-6 , 10 ] \
    --shrink-factors 2x1 \
    --smoothing-sigmas 1x0vox
   done
           
                                   
for n in $VOLUMES; do
    ConvertTransformFile 3 NORDIC_bold_"$n"0GenericAffine.mat NORDIC_bold_"$n"_ants2itk.mat --hm --ras
    c3d_affine_tool -ref NORDIC_bold_slice_Vol_Remove_1000.nii.gz \
      -src NORDIC_bold_slice_Vol_Remove_$n.nii.gz \
      NORDIC_bold_"$n"_ants2itk.mat -ras2fsl \
      -o NORDIC_bold_"$n"_itk2fsl.mat
done                    
                    
                                        
# $FSLDIR refers to fsl root directory which is /usr/local/fsl in Linux systems. The below commands save motion estimates into text files and shows the motion plots                      
for n in $VOLUMES; do ${FSLDIR}/bin/avscale --allparams NORDIC_bold_"$n"_itk2fsl.mat NORDIC_bold_slice_Vol_Remove_1000.nii.gz | grep "Translations" | awk '{print $5 " " $6 " " $7}';done > translation.txt          
for n in $VOLUMES; do ${FSLDIR}/bin/avscale --allparams NORDIC_bold_"$n"_itk2fsl.mat NORDIC_bold_slice_Vol_Remove_1000.nii.gz | grep "Rotation Angles" | awk '{print $6 " " $7 " " $8}';done > Rotation.txt  

paste translation.txt Rotation.txt > Motion.params

${FSLDIR}/bin/fsl_tsplot -i Motion.params -t 'Translations (mm)' -u 1 --start=1 --finish=3 -a x,y,z -w 640 -h 144 -o Motion_translations.png
${FSLDIR}/bin/fsl_tsplot -i Motion.params -t 'Rotations (deg)' -u 1 --start=4 --finish=6 -a x,y,z -w 640 -h 144 -o Motion_rotations.png

for n in $VOLUMES; do ${FSLDIR}/bin/rmsdiff NORDIC_bold_1000_itk2fsl.mat NORDIC_bold_"$n"_itk2fsl.mat NORDIC_bold_slice_Vol_Remove_1000.nii.gz; done >> Motion.rmsabs
for n in $(seq 1001 1149); do F=$(($n-1)); ${FSLDIR}/bin/rmsdiff NORDIC_bold_"$F"_itk2fsl.mat  NORDIC_bold_"$n"_itk2fsl.mat NORDIC_bold_slice_Vol_Remove_1000.nii.gz; done >> Motion.rmsrel

${FSLDIR}/bin/fsl_tsplot -i Motion.rmsabs,Motion.rmsrel  -t 'Mean Displacements (mm)' -u 1 -a absolute,relative  -w 640 -h 144 -o Motion_rms.png

ImageMath 4 NORDIC_bold_MoCo.nii.gz TimeSeriesAssemble 3 0 *_Warped.nii.gz 

# now repeat the above for vaso
cd ../NORDIC_vaso_slice_Vol_Remove_split/
 
for n in $VOLUMES; do
  antsRegistration \
    --verbose 1 \
    --float 1 \
    --dimensionality 3 \
    --use-histogram-matching 1 \
    --interpolation LanczosWindowedSinc \
    --collapse-output-transforms 1 \
    --output [ NORDIC_vaso_$n, NORDIC_vaso_"$n"_slice_Vol_Remove_Warped.nii.gz , 1 ] \
    --winsorize-image-intensities [ 0.005 , 0.995 ] \
    --initial-moving-transform [ NORDIC_vaso_slice_Vol_Remove_1000.nii.gz , NORDIC_vaso_slice_Vol_Remove_$n.nii.gz , 1 ] \
    --transform Rigid[0.1] \
    --metric MI[ NORDIC_vaso_slice_Vol_Remove_1000.nii.gz , NORDIC_vaso_slice_Vol_Remove_$n.nii.gz , 1 , 64 , Regular , 0.25] \
    --convergence [ 500x250 , 1e-6 , 10 ] \
    --shrink-factors 2x1 \
    --smoothing-sigmas 1x0vox
done
           
for n in $VOLUMES; do
    ConvertTransformFile 3 NORDIC_vaso_"$n"0GenericAffine.mat NORDIC_vaso_"$n"_ants2itk.mat --hm --ras
    c3d_affine_tool -ref NORDIC_vaso_slice_Vol_Remove_1000.nii.gz \
      -src NORDIC_vaso_slice_Vol_Remove_$n.nii.gz \
      NORDIC_vaso_"$n"_ants2itk.mat -ras2fsl \
      -o NORDIC_vaso_"$n"_itk2fsl.mat
done               
                    
for n in $VOLUMES; do ${FSLDIR}/bin/avscale --allparams NORDIC_vaso_"$n"_itk2fsl.mat NORDIC_vaso_slice_Vol_Remove_1000.nii.gz | grep "Translations" | awk '{print $5 " " $6 " " $7}';done > translation.txt         
for n in $VOLUMES; do ${FSLDIR}/bin/avscale --allparams NORDIC_vaso_"$n"_itk2fsl.mat NORDIC_vaso_slice_Vol_Remove_1000.nii.gz | grep "Rotation Angles" | awk '{print $6 " " $7 " " $8}';done > Rotation.txt  

paste translation.txt Rotation.txt > Motion.params

${FSLDIR}/bin/fsl_tsplot -i Motion.params -t 'Translations (mm)' -u 1 --start=1 --finish=3 -a x,y,z -w 640 -h 144 -o Motion_translations.png
${FSLDIR}/bin/fsl_tsplot -i Motion.params -t 'Rotations (deg)' -u 1 --start=4 --finish=6 -a x,y,z -w 640 -h 144 -o Motion_rotations.png

for n in $VOLUMES; do ${FSLDIR}/bin/rmsdiff NORDIC_vaso_1000_itk2fsl.mat NORDIC_vaso_"$n"_itk2fsl.mat NORDIC_vaso_slice_Vol_Remove_1000.nii.gz; done >> Motion.rmsabs
for n in $(seq 1001 1149); do F=$(($n-1)); ${FSLDIR}/bin/rmsdiff NORDIC_vaso_"$F"_itk2fsl.mat  NORDIC_vaso_"$n"_itk2fsl.mat NORDIC_vaso_slice_Vol_Remove_1000.nii.gz; done >> Motion.rmsrel

${FSLDIR}/bin/fsl_tsplot -i Motion.rmsabs,Motion.rmsrel  -t 'Mean Displacements (mm)' -u 1 -a absolute,relative  -w 640 -h 144 -o Motion_rms.png

ImageMath 4 NORDIC_vaso_MoCo.nii.gz TimeSeriesAssemble 3 0 *_Warped.nii.gz 
cd ../../../../       
 
# In session_3 of sub 001, there are 3 runs. Creat a for loop to run the same motion correction as in sub 002. 

cd 001/func/session_3/

RUNS=("Run7" "Run8" "Run9")
for RUN in "${RUNS[@]}"; do
  RUN_DIR="${RUN}"
  cd "$RUN_DIR" || continue
  echo "Processing $RUN_DIR..."
  mkdir NORDIC_bold_slice_Vol_Remove_split NORDIC_vaso_slice_Vol_Remove_split
  
  ImageMath 4 NORDIC_bold_slice_Vol_Remove_split/NORDIC_bold_slice_Vol_Remove_.nii.gz TimeSeriesDisassemble NORDIC_bold_slice_Vol_Remove.nii.gz
  ImageMath 4 NORDIC_vaso_slice_Vol_Remove_split/NORDIC_vaso_slice_Vol_Remove_.nii.gz TimeSeriesDisassemble NORDIC_vaso_slice_Vol_Remove.nii.gz
  
  cd NORDIC_bold_slice_Vol_Remove_split || continue
  for n in $VOLUMES; do 
  antsRegistration \
    --verbose 1 \
    --float 1 \
    --dimensionality 3 \
    --use-histogram-matching 1 \
    --interpolation LanczosWindowedSinc \
    --collapse-output-transforms 1 \
    --output [ NORDIC_bold_$n, NORDIC_bold_"$n"_slice_Vol_Remove_Warped.nii.gz , 1 ] \
    --winsorize-image-intensities [ 0.005 , 0.995 ] \
    --initial-moving-transform [ NORDIC_bold_slice_Vol_Remove_1000.nii.gz , NORDIC_bold_slice_Vol_Remove_$n.nii.gz , 1 ] \
    --transform Rigid[0.1] \
    --metric MI[ NORDIC_bold_slice_Vol_Remove_1000.nii.gz , NORDIC_bold_slice_Vol_Remove_$n.nii.gz , 1 , 64 , Regular , 0.25] \
    --convergence [ 500x250 , 1e-6 , 10 ] \
    --shrink-factors 2x1 \
    --smoothing-sigmas 1x0vox
   done
   
   for n in $VOLUMES; do
    ConvertTransformFile 3 NORDIC_bold_"$n"0GenericAffine.mat NORDIC_bold_"$n"_ants2itk.mat --hm --ras
    c3d_affine_tool -ref NORDIC_bold_slice_Vol_Remove_1000.nii.gz \
      -src NORDIC_bold_slice_Vol_Remove_$n.nii.gz \
      NORDIC_bold_"$n"_ants2itk.mat -ras2fsl \
      -o NORDIC_bold_"$n"_itk2fsl.mat
  done  
         
  for n in $VOLUMES; do ${FSLDIR}/bin/avscale --allparams NORDIC_bold_"$n"_itk2fsl.mat NORDIC_bold_slice_Vol_Remove_1000.nii.gz | grep "Translations" | awk '{print $5 " " $6 " " $7}';done >    translation.txt          
  for n in $VOLUMES; do ${FSLDIR}/bin/avscale --allparams NORDIC_bold_"$n"_itk2fsl.mat NORDIC_bold_slice_Vol_Remove_1000.nii.gz | grep "Rotation Angles" | awk '{print $6 " " $7 " " $8}';done >  Rotation.txt  

  paste translation.txt Rotation.txt > Motion.params
  
  ${FSLDIR}/bin/fsl_tsplot -i Motion.params -t 'Translations (mm)' -u 1 --start=1 --finish=3 -a x,y,z -w 640 -h 144 -o Motion_translations.png
  ${FSLDIR}/bin/fsl_tsplot -i Motion.params -t 'Rotations (deg)' -u 1 --start=4 --finish=6 -a x,y,z -w 640 -h 144 -o Motion_rotations.png
  
  for n in $VOLUMES; do ${FSLDIR}/bin/rmsdiff NORDIC_bold_1000_itk2fsl.mat NORDIC_bold_"$n"_itk2fsl.mat NORDIC_bold_slice_Vol_Remove_1000.nii.gz; done >> Motion.rmsabs
  for n in $(seq 1001 1149); do F=$(($n-1)); ${FSLDIR}/bin/rmsdiff NORDIC_bold_"$F"_itk2fsl.mat  NORDIC_bold_"$n"_itk2fsl.mat NORDIC_bold_slice_Vol_Remove_1000.nii.gz; done >> Motion.rmsrel

  ${FSLDIR}/bin/fsl_tsplot -i Motion.rmsabs,Motion.rmsrel  -t 'Mean Displacements (mm)' -u 1 -a absolute,relative  -w 640 -h 144 -o Motion_rms.png

  ImageMath 4 NORDIC_bold_MoCo.nii.gz TimeSeriesAssemble 3 0 *_Warped.nii.gz 

# now repeat the above for vaso

  cd ../NORDIC_vaso_slice_Vol_Remove_split || continue
  for n in $VOLUMES; do
  antsRegistration \
    --verbose 1 \
    --float 1 \
    --dimensionality 3 \
    --use-histogram-matching 1 \
    --interpolation LanczosWindowedSinc \
    --collapse-output-transforms 1 \
    --output [ NORDIC_vaso_$n, NORDIC_vaso_"$n"_slice_Vol_Remove_Warped.nii.gz , 1 ] \
    --winsorize-image-intensities [ 0.005 , 0.995 ] \
    --initial-moving-transform [ NORDIC_vaso_slice_Vol_Remove_1000.nii.gz , NORDIC_vaso_slice_Vol_Remove_$n.nii.gz , 1 ] \
    --transform Rigid[0.1] \
    --metric MI[ NORDIC_vaso_slice_Vol_Remove_1000.nii.gz , NORDIC_vaso_slice_Vol_Remove_$n.nii.gz , 1 , 64 , Regular , 0.25] \
    --convergence [ 500x250 , 1e-6 , 10 ] \
    --shrink-factors 2x1 \
    --smoothing-sigmas 1x0vox
  done
  
  for n in $VOLUMES; do
    ConvertTransformFile 3 NORDIC_vaso_"$n"0GenericAffine.mat NORDIC_vaso_"$n"_ants2itk.mat --hm --ras
    c3d_affine_tool -ref NORDIC_vaso_slice_Vol_Remove_1000.nii.gz \
      -src NORDIC_vaso_slice_Vol_Remove_$n.nii.gz \
      NORDIC_vaso_"$n"_ants2itk.mat -ras2fsl \
      -o NORDIC_vaso_"$n"_itk2fsl.mat
  done 
  
  for n in $VOLUMES; do ${FSLDIR}/bin/avscale --allparams NORDIC_vaso_"$n"_itk2fsl.mat NORDIC_vaso_slice_Vol_Remove_1000.nii.gz | grep "Translations" | awk '{print $5 " " $6 " " $7}';done >  translation.txt      
     
  for n in $VOLUMES; do ${FSLDIR}/bin/avscale --allparams NORDIC_vaso_"$n"_itk2fsl.mat NORDIC_vaso_slice_Vol_Remove_1000.nii.gz | grep "Rotation Angles" | awk '{print $6 " " $7 " " $8}';done >  Rotation.txt  

  paste translation.txt Rotation.txt > Motion.params

  ${FSLDIR}/bin/fsl_tsplot -i Motion.params -t 'Translations (mm)' -u 1 --start=1 --finish=3 -a x,y,z -w 640 -h 144 -o Motion_translations.png
  ${FSLDIR}/bin/fsl_tsplot -i Motion.params -t 'Rotations (deg)' -u 1 --start=4 --finish=6 -a x,y,z -w 640 -h 144 -o Motion_rotations.png

  for n in $VOLUMES; do ${FSLDIR}/bin/rmsdiff NORDIC_vaso_1000_itk2fsl.mat NORDIC_vaso_"$n"_itk2fsl.mat NORDIC_vaso_slice_Vol_Remove_1000.nii.gz; done >> Motion.rmsabs
  for n in $(seq 1001 1149); do F=$(($n-1)); ${FSLDIR}/bin/rmsdiff NORDIC_vaso_"$F"_itk2fsl.mat  NORDIC_vaso_"$n"_itk2fsl.mat NORDIC_vaso_slice_Vol_Remove_1000.nii.gz; done >> Motion.rmsrel

  ${FSLDIR}/bin/fsl_tsplot -i Motion.rmsabs,Motion.rmsrel  -t 'Mean Displacements (mm)' -u 1 -a absolute,relative  -w 640 -h 144 -o Motion_rms.png 
  
  ImageMath 4 NORDIC_vaso_MoCo.nii.gz TimeSeriesAssemble 3 0 *_Warped.nii.gz 
  
  # Return to Session_3 for next run
  cd ../../
done 
       

# preparation to run BOCO function from Laynii package
 
cat FOLDERS_no_PA.txt | while read line; do
  if [ -d ../reverse_phase ]; then
    3dTcat -prefix $line/vaso_bold_combined.nii.gz \
      $line/NORDIC_vaso_slice_Vol_Remove_split/NORDIC_vaso_MoCo_DisCor_merged.nii.gz \
      $line/NORDIC_bold_slice_Vol_Remove_split/NORDIC_bold_MoCo_DisCor_merged.nii.gz
  else
    3dTcat -prefix $line/vaso_bold_combined.nii.gz \
      $line/NORDIC_vaso_slice_Vol_Remove_split/NORDIC_vaso_MoCo.nii.gz \
      $line/NORDIC_bold_slice_Vol_Remove_split/NORDIC_bold_MoCo.nii.gz
  fi
done
 
cat FFOLDERS_no_PA.txt | while read line; do 3dTstat -cvarinv -prefix $line/T1_weighted_vaso_bold.nii.gz $line/vaso_bold_combined.nii.gz; done 


cat FFOLDERS_no_PA.txt | while read line; do
  if [ -d ../reverse_phase ]; then
    3dTstat -mean -prefix $line/vaso_mean.nii.gz $line/NORDIC_vaso_slice_Vol_Remove_split/NORDIC_vaso_MoCo_DisCor_merged.nii.gz
    3dTstat -mean -prefix $line/bold_mean.nii.gz $line/NORDIC_bold_slice_Vol_Remove_split/NORDIC_bold_MoCo_DisCor_merged.nii.gz
    3dcalc -a $line/vaso_mean.nii.gz -b $line/bold_mean.nii.gz -prefix $line/T1_weighted_vaso_bold_divided.nii.gz -expr 'within(((a-b)/(a+b)+1),0,1.2)*(a-b)/(a+b)'
    3dUpsample -datum short -prefix $line/vaso_upsamp.nii.gz -n 2 -input $line/NORDIC_vaso_slice_Vol_Remove_split/NORDIC_vaso_MoCo_DisCor_merged.nii.gz
    3dUpsample -datum short -prefix $line/bold_upsamp.nii.gz -n 2 -input $line/NORDIC_bold_slice_Vol_Remove_split/NORDIC_bold_MoCo_DisCor_merged.nii.gz
  else
    3dTstat -mean -prefix $line/vaso_mean.nii.gz $line/NORDIC_vaso_slice_Vol_Remove_split/NORDIC_vaso_MoCo.nii.gz
    3dTstat -mean -prefix $line/bold_mean.nii.gz $line/NORDIC_bold_slice_Vol_Remove_split/NORDIC_bold_MoCo.nii.gz
    3dcalc -a $line/vaso_mean.nii.gz -b $line/bold_mean.nii.gz -prefix $line/T1_weighted_vaso_bold_divided.nii.gz -expr 'within(((a-b)/(a+b)+1),0,1.2)*(a-b)/(a+b)'
    3dUpsample -datum short -prefix $line/vaso_upsamp.nii.gz -n 2 -input $line/NORDIC_vaso_slice_Vol_Remove_split/NORDIC_vaso_MoCo.nii.gz
    3dUpsample -datum short -prefix $line/bold_upsamp.nii.gz -n 2 -input $line/NORDIC_bold_slice_Vol_Remove_split/NORDIC_bold_MoCo.nii.gz
  fi
done

cat FFOLDERS_no_PA.txt | while read line; do
   NumVol=`3dinfo -nv "$line/bold_upsamp.nii.gz"`
   b="$((NumVol - 2))"
   3dTcat -prefix $line/bold_upsamp_shifted.nii.gz \
   $line/bold_upsamp.nii.gz'[0]' \
   $line/bold_DisCor_MoCor_upsamp.nii.gz[0.."b"]    
done

cat FFOLDERS_no_PA.txt | while read line; do
  cd $line/ || continue
  LN_BOCO -Nulled vaso_upsamp.nii.gz -BOLD bold_upsamp_shifted.nii.gz
  gzip VASO_LN.nii
  rm bold_upsamp_shifted.nii.gz
  cd "$START_DIR"
done

cat FFOLDERS_no_PA.txt | while read line; do 3drefit -TR 3 $line/VASO_LN.nii.gz; done 

cat FFOLDERS_no_PA.txt | while read line; do 3dcalc -a $line/VASO_LN.nii.gz'[0..$(2)]' -expr 'a' -prefix $line/VASO_LN_origVols.nii.gz; done 

# remove the noisy background of the BOLD-corrected VASO, by generating a mask of BOLD data and then applying this mask to VASO 
cat FFOLDERS_no_PA.txt | while read line; do
  if [ -d ../reverse_phase ]; then
    3dAutomask -prefix $line/bold_mask.nii.gz -apply_prefix $line/bold_MoCo_DisCor_maskApplied.nii.gz -dilate 2 $line/NORDIC_bold_slice_Vol_Remove_split/NORDIC_bold_MoCo_DisCor_merged.nii.gz
    fslmaths $line/VASO_LN_origVols.nii.gz -mul $line/bold_mask.nii.gz $line/VASO_LN_origVols_maskApplied.nii.gz
  else
  3dAutomask -prefix $line/bold_mask.nii.gz -apply_prefix $line/bold_MoCo_maskApplied.nii.gz -dilate 2 $line/NORDIC_bold_slice_Vol_Remove_split/NORDIC_bold_MoCo.nii.gz
  fslmaths $line/VASO_LN_origVols.nii.gz -mul $line/bold_mask.nii.gz $line/VASO_LN_origVols_maskApplied.nii.gz
  fi
done

# MPRAGIZE the anatomy to remove the noisy background of UNI image. Requires Presurfer to be added to MATLAB path in advance.   
cat FOLDERS-ANAT.txt | while read line; do cd $line/
UNI='MP2RAGE-UNI_defaced.nii.gz'
INV2='MP2RAGE-INV2_defaced.nii.gz'
matlab -nodisplay -nodesktop -r "[mprageised_im, wmseg_im] = MPRAGEise('$UNI', '$INV2');quit"; cd "$START_DIR"; done 

# Co-registration of functional data to corresponding anatomical image

while IFS= read -r line1 && IFS= read -r line2 <&3; do
  cd $line1/ || continue
  if [ -d ../reverse_phase ]; then
    align_epi_anat.py \
      -anat "$START_DIR"/$line2/MP2RAGE-UNI_defaced_MPRAGEised.nii.gz \
      -anat_has_skull yes \
      -epi NORDIC_bold_slice_Vol_Remove_split/NORDIC_bold_MoCo_DisCor_merged.nii.gz \
      -epi_base 0 \
      -partial_coverage \
      -epi2anat \
      -prep_off \
      -epi_strip 3dAutomask \
      -cmass cmass \
      -giant_move \
      -suffix _coreg
   
    3dAFNItoNIFTI -prefix NORDIC_bold_MoCo_DisCor_merged_coreg.nii.gz NORDIC_bold_MoCo_DisCor_merged_coreg+orig.
    3dresample -master "$START_DIR"/$line2/MP2RAGE-UNI_defaced_MPRAGEised.nii.gz -prefix NORDIC_bold_coreg_resampled.nii.gz -input NORDIC_bold_MoCo_DisCor_merged_coreg.nii.gz
  else 
    align_epi_anat.py \
      -anat "$START_DIR"/$line2/MP2RAGE-UNI_defaced_MPRAGEised.nii.gz \
      -anat_has_skull yes \
      -epi NORDIC_bold_slice_Vol_Remove_split/NORDIC_bold_MoCo.nii.gz \
      -epi_base 0 \
      -partial_coverage \
      -epi2anat \
      -prep_off \
      -epi_strip 3dAutomask \
      -cmass cmass \
      -giant_move \
      -suffix _coreg
      
    3dAFNItoNIFTI -prefix NORDIC_bold_MoCo_coreg.nii.gz NORDIC_bold_MoCo_coreg+orig.  
    3dresample -master "$START_DIR"/$line2/MP2RAGE-UNI_defaced_MPRAGEised.nii.gz -prefix NORDIC_bold_coreg_resampled.nii.gz -input NORDIC_bold_MoCo_DisCor_merged_coreg.nii.gz
  fi
  cd "$START_DIR"
done <  FOLDERS_no_PA.txt 3< FOLDERS-ANAT.txt 
   
# coregistration of vaso to anatomy involves manual approach using ITK-snap tools, to do so, first get the mean VASO image, split the volumes and then once the transformation matrix is generated from ITK-snap alignment 

cat FOLDERS_no_PA.txt | while read line; do
  mkdir $line/vaso_split
  3dTstat -mean -prefix $line/VASO_LN_origVols_masked_mean.nii.gz $line/VASO_LN_origVols_maskApplied.nii.gz
  cp $line/VASO_LN_origVols_maskApplied.nii.gz $line/vaso_split/
  cd $line/vaso_split/ || continue
  fslsplit VASO_LN_origVols_maskApplied.nii.gz vaso -t
  cd "$START_DIR"
done

echo "Waiting for VASO_coreg.txt to be created via manual alignment with ITK-snap. Press Ctrl+C to cancel"

while [ ! -f VASO_coreg.txt ]; do
  sleep 3600  # check every hour 
done

echo "Manual registration completed. Continuing with script..."

## Once VASO_coreg is created, then apply it to all VASO volumes. If the process was cancelled, create a new .sh script by pasting the lines below.  

while IFS= read -r line1 && IFS= read -r line2 <&3; do
cd $line1/vaso_split/ || continue
for n in $VOLUMES; do
antsApplyTransforms --interpolation BSpline[5] -d 3 -i vaso$n.nii.gz -r "$START_DIR"/$line2/MP2RAGE-UNI_defaced_MPRAGEised.nii.gz -t VASO_coreg.txt -o vaso$n-coreg.nii.gz
done
fslmerge -t VASO_LN_coreg.nii.gz *-coreg.nii.gz
cd "$START_DIR"
done < FOLDERS_no_PA.txt 3< FOLDERS-ANAT.txt

echo "Visually Investigate the Goodness of Alignment for BOLD and VASO "

# Merge all co-registered BOLD data into a single NIfTI (.nii) file, which will serve as the input for a MATLAB-based sampling script. While the code can be executed separately for each run followed by merging the resulting .mat files, we choose instead to combine all runs beforehand to generate a single .mat file as the final output. Similarly, merge all co-registered VASO data. Note that this step requires a lot of RAM around 500 GB.  

subjects=(001 002 003 004 005 006)

for subj in "${subjects[@]}"; do
    echo "Processing subject $subj..."
    
    # Find coreg BOLD files
    files=$(find $subj -type f -name "NORDIC_bold_coreg_resampled.nii.gz" | sort)
   
    outdir="${subj}/func"

    3dTcat -prefix "${outdir}/NORDIC_BOLD_coreg_merged_all.nii.gz" $files

    echo "Done with $subj"
done
  
# Apply the same for VASO 
for subj in "${subjects[@]}"; do
    echo "Processing subject $subj..."
    
    # Find coreg BOLD files
    files=$(find $subj -type f -name "VASO_LN_coreg.nii.gz" | sort)
   
    outdir="${subj}/func"

    3dTcat -prefix "${outdir}/NORDIC_VASO_coreg_merged_all.nii.gz" $files

    echo "Done with $subj"
done

# Following co-registration, obtain usual QC metrics including mean image and tSNR maps 

cat FOLDERS_no_PA.txt | while read line; do
  mkdir $line/QC_bold $line/QC_vaso
  fslmaths $line/vaso_split/VASO_LN_coreg.nii.gz -Tmean $line/QC_vaso/vaso_mean.nii.gz
  fslmaths $line/vaso_split/VASO_LN_coreg.nii.gz -Tstd $line/QC_vaso/vaso_tsd.nii.gz
  fslmaths $line/QC_vaso/vaso_mean.nii.gz -div $line/QC_vaso/vaso_tsd.nii.gz $line/QC_vaso/vaso_tsnr.nii.gz
  fslmaths $line/NORDIC_bold_coreg_resampled.nii.gz -Tmean $line/QC_bold/bold_mean.nii.gz
  fslmaths $line/NORDIC_bold_coreg_resampled.nii.gz -Tstd $line/QC_bold/bold_tsd.nii.gz
  fslmaths $line/QC_bold/bold_mean.nii.gz -div $line/QC_bold/bold_tsd.nii.gz $line/QC_bold/bold_tsnr.nii.gz
done
  

#--------------------------------------------------
# Section2: HC segmentation to subfields using HippUnfold package (here is run through singularity container) 
#--------------------------------------------------
echo "Starting with HC segmentation"

while IFS= read -r line1 && IFS= read -r line2 <&3; do 
 mkdir $line1/HUinput $line1/HUoutput_T1
 mkdir -p $line1/HUinput/sub-$line2/anat
 cp $line1/MP2RAGE-UNI_defaced_MPRAGEised.nii.gz $line1/HUinput/sub-$line2/anat/T1w.nii.gz  
done <  FOLDERS-ANAT.txt 3< Hippunfold_indices.txt

# khanlab_hippunfold_latest.sif is located in my home directory, change the path accordingly otherwise it will complain 
SIF_PATH="/home/kahmadi/khanlab_hippunfold_latest.sif"

# Check if the SIF file exists before proceeding
if [ ! -f "$SIF_PATH" ]; then
  echo "ERROR: Cannot find SIF file at: $SIF_PATH"
  echo "Please edit the script and update the correct path to khanlab_hippunfold_latest.sif"
  exit 1
fi

cat FOLDERS-ANAT.txt | while read line; do
  singularity run -e "$SIF_PATH" $line/HUinput/ $line/HUoutput_T1 participant -p --cores all --modality T1w
done 

#--------------------------------------------------
# Section3: Whole brain segmentation using FSL FAST to get cropped masks of CSF, WM (to be used with 'aCompCor') and GM (to be used later during cluster-level inference in GLM analysis with SPM)
#--------------------------------------------------

echo "preparation for FAST segmentation"

cat FOLDERS-ANAT.txt | while read line; do 
  mkdir $line/FSL_fast; cp $line/MP2RAGE-UNI_defaced_MPRAGEised.nii.gz $line/FSL_fast
  # Skull stripping (adjust -f if needed)
  bet $line/FSL_fast/MP2RAGE-UNI_defaced_MPRAGEised.nii.gz $line/FSL_fast/T1 -f 0.3 # double check the output, you may need to change -f value to 0.2 in some cases 
  
  fast -n 3 -t 1 -o $line/FSL_fast/fast_T1 -b bias -B $line/FSL_fast/biasedRemoved $line/FSL_fast/T1.nii.gz
done 

# get masks of each tissue type

cat FOLDERS-ANAT.txt | while read line; do
  fslmaths $line/FSL_fast/fast_T1_pveseg.nii.gz -uthr 1 $line/FSL_fast/csf_new.nii.gz
  fslmaths $line/FSL_fast/fast_T1_pveseg.nii.gz -thr 3 $line/FSL_fast/WM_new.nii.gz
  fslmaths $line/FSL_fast/fast_T1_pveseg.nii.gz -uthr 2 -thr 2 $line/FSL_fast/GM_new.nii.gz
  fslmaths $line/FSL_fast/WM_new.nii.gz -bin $line/FSL_fast/WM_new_bin.nii.gz
  fslmaths $line/FSL_fast/GM_new.nii.gz -bin $line/FSL_fast/GM_new_bin.nii.gz
done 

# crop them to macth functional slab 
awk -F/ '!seen[$1]++' FOLDERS_no_PA.txt > unique_no_PA.txt
paste FOLDERS-ANAT.txt unique_no_PA.txt | while IFS=$'\t' read -r line1 line2; do
  fslmaths $line1/FSL_fast/WM_new_bin.nii.gz -mas $line2/QC_vaso/vaso_mean.nii.gz $line1/FSL_fast/WM_new_bin_cropped.nii.gz
  fslmaths $line1/FSL_fast/csf_new.nii.gz -mas $line2/QC_vaso/vaso_mean.nii.gz $line1/FSL_fast/csf_new_cropped.nii.gz
  fslmaths $line1/FSL_fast/GM_new_bin.nii.gz -mas $line2/QC_vaso/vaso_mean.nii.gz $line1/FSL_fast/GM_new_bin_cropped.nii.gz
done

# The outputs of FAST segmenetation may contain a few csf or wm voxels inside HC i.e., incorrectly segmented voxels. To make sure that no HC voxel is included in those masks erode them before feeding them to aCompCor. Manual edits might be necessary at this stage.  

cat FOLDERS-ANAT.txt | while read line; do
  fslmaths $line/FSL_fast/WM_new_bin_cropped.nii.gz -kernel gauss 0.8 -ero WM_new_bin_cropped_eroded.nii.gz
  fslmaths $line/FSL_fast/csf_new_cropped.nii.gz -kernel gauss 0.8 -ero csf_new_cropped_eroded.nii.gz 
done 


#--------------------------------------------------
# Section4: Create a study-specific T1-template to be able to run 2nd-level analysis
#--------------------------------------------------

echo "creating T1-template"

mkdir Template_from_T1s

while IFS= read -r line1 && IFS= read -r line2 <&3; do 
 cp $line1/FSL_fast/T1.nii.gz Template_from_T1s/T1_sub_"$line2".nii.gz
done <  FOLDERS-ANAT.txt 3< Hippunfold_indices.txt

cd Template_from_T1s
antsMultivariateTemplateConstruction2.sh -d 3 -o SST_ -c 2 -j 4 -k 1 -r 1 -t SyN -m CC T1_sub0*.nii.gz


# In case you want to check whether the second-level results fall into HC subfields and restrict the cluster-level inferences to the GM mask of the template, use Hippunfold and FAST on the T1 template.

mkdir FSL_fast HUinput HUoutput_T1
mkdir -p HUinput/sub-0000/anat
cp SST_template0.nii.gz HUinput/sub-0000/anat/T1w.nii.gz 
singularity run -e "$SIF_PATH" HUinput/ HUoutput_T1 participant -p --cores all --modality T1w


## GET GM mask of the template T1 and run Hippunfold on template T1
cp SST_template0.nii.gz FSL_fast/

# Running FAST directly on SST_template0.nii.gz may fail because the background of template image may have non-zero values. To solve this, first create a masked brain using 3dAutomask and then run fast

cd FSL_fast/
3dAutomask -prefix masked.nii.gz  -apply_prefix SST_template0_maksed.nii.gz -dilate 2 SST_template0.nii.gz   

fast -n 3 -t 1 -o fast_T1 -b bias -B biasedRemoved SST_template0_maksed.nii.gz
fslmaths fast_T1_pveseg.nii.gz -uthr 2 -thr 2 GM_new.nii.gz
fslmath GM_new.nii.gz -bin GM_new_bin.nii.gz

cd "$START_DIR"
 echo "FINISHED"
