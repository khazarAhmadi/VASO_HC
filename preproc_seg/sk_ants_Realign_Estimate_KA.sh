#! /bin/bash
################################################################################
#
# PREPARATIONS
#
################################################################################

if [[ ${#ANTSPATH} -le 3 ]]; then
    setPath >&2
fi

ANTS=${ANTSPATH}/antsRegistration

if [[ ! -s ${ANTS} ]]; then
    echo "ANTS can't be found. Please (re)define $ANTSPATH in your environment."
    exit
fi

################################################################################
# Simple formatting

bold=$(tput bold)
normal=$(tput sgr0)

################################################################################

function Help() {
    cat <<HELP

Usage:

$(basename $0) ${bold}-f${normal} Fixed Image ${bold} ${bold}-x${normal} Fixed Mask ${bold} -a${normal} Path to 4D fMRI data ${bold}-b${normal} Path to opposite PE data ${bold}-n${normal} number of threads

--------------------------------------------------------------------------------
Input arguments:

    -f: Fixed Image (created if not specified. Specify as /path/to/data/fixed.nii.gz)

    -x: Fixed Mask (created if not specified. Specify as /path/to/data/fixedMask.nii.gz)

    -t: TR in seconds (e.g. 2.0)

    -a: Path to 4D fMRI data (required. Specify as /path/to/data/fMRI.nii.gz)

    -b: Path to opposite PE data (required. Specify as /path/to/data/oPE.nii.gz)

    -n: Number of threads (default = '16'. Specify more threads if available.)

    -s: staged processing

--------------------------------------------------------------------------------

Example:

$(basename $0) -t 2.0 -a /home/user/folder/file.nii.gz -b /home/user/folder/file_ope.nii.gz -n 24

This command creates the fixed file, fixed mask, estimates motion & distortion-field,
saves all required transformations for the Reslice script. -s = 1, stops after making the fixed images.

--------------------------------------------------------------------------------
Script was created by: Sriranga Kashyap (09-2020)
--------------------------------------------------------------------------------
Requires ANTs to be installed and $ANTSPATH defined in your environment.
--------------------------------------------------------------------------------
ANTs can be downloaded here: https://github.com/ANTsX/ANTs
References to cite:
    1) http://www.ncbi.nlm.nih.gov/pubmed/20851191
    2) http://www.frontiersin.org/Journal/10.3389/fninf.2013.00039/abstract
--------------------------------------------------------------------------------

HELP
    exit 1
}

################################################################################
function reportParameters() {
    cat <<REPORTPARAMETERS

--------------------------------------------------------------------------------
    ${bold} Processes Initialised ${normal}
--------------------------------------------------------------------------------
    ANTSPATH is $ANTSPATH

    Fixed Image         : ${bold} ${fixed_data_name}.nii.gz ${normal}
    Fixed Mask          : ${bold} ${fixed_mask_data_name}.nii.gz ${normal}
    TR                  : ${bold} $tr ${normal}
    4D fMRI data        : ${bold} ${func_data_name}.nii.gz ${normal}
    Opposite PE data    : ${bold} ${ope_data_name}.nii.gz ${normal}
    Number of threads   : ${bold} $n_threads ${normal}
--------------------------------------------------------------------------------

REPORTPARAMETERS
}

################################################################################

if [[ "$1" == "-h" || $# -eq 0 ]]; then
    Help >&2
fi

################################################################################
#
# DEFAULTS
#
################################################################################
interpolation_type=LanczosWindowedSinc
n_threads=24
iterations=6
################################################################################
#
# PARSE INPUT ARGUMENTS
#
################################################################################

while getopts "f:x:t:a:b:n:z:s:" OPT; do
    case $OPT in
    h) #help
        Help
        exit 0
        ;;
    f) # fixed data path
        fixed_data_path=$OPTARG
        ;;
    x) # fixed data mask
        fixed_data_mask=$OPTARG
        ;;
    t) # tr
        tr=$OPTARG
        ;;
    a) # 4D fMRI data path
        func_data_path=$OPTARG
        ;;
    b) # opposite phase-encoding data path
        ope_data_path=$OPTARG
        ;;
    n) # transformation type
        n_threads=$OPTARG
        ;;
    z) # template reference
        reftemp_data_path=$OPTARG
        ;;
    s) # staged processing
        staged_proc=$OPTARG
        ;;
    \?) # report error
        echo "$HELP" >&2
        exit 1
        ;;
    esac
done

################################################################################
#
# SET NUMBER OF THREADS
#
################################################################################

ORIGINALNUMBEROFTHREADS=${ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS}
ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$n_threads
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS

################################################################################
#
# START PROCESSING FILES
#
################################################################################

start_time=$(date +%s)

# GET INPUTS
if [ -z "$staged_proc" ]; then
    echo "--------------------------------------------------------------------------------"
    echo "++++ Not doing staged processing."
else
    echo "--------------------------------------------------------------------------------"
    echo "++++ Staged processing enabled."

fi

if [ -z "$fixed_data_path" ]; then
    echo "--------------------------------------------------------------------------------"
    echo "++++ Fixed Image not specified."
    # initialise file names
    fxd_full=$(basename ${func_data_path})
    fxd_prefix=${fxd_full%%.*}
    fixed_data=${fxd_prefix}_fixed.nii.gz
    fixed_data_name=${fxd_prefix}_fixed
else
    fixed_data=$(basename ${fixed_data_path})
    fixed_data_name=${fixed_data%%.*}
fi

if [ -z "$fixed_data_mask" ]; then
    echo " "
    echo "++++ Fixed Mask not specified."
    # initialise file names
    fxd_mask_full=$(basename ${func_data_path})
    fxd_mask_prefix=${fxd_mask_full%%.*}
    fixed_mask_data=${fxd_mask_prefix}_fixedMask.nii.gz
    fixed_mask_data_name=${fxd_mask_prefix}_fixedMask

else
    fixed_mask=$(basename ${fixed_data_mask})
    fixed_mask_data_name=${fixed_mask%%.*}
fi

func_data=$(basename ${func_data_path})
func_data_name=${func_data%%.*}

ope_data=$(basename ${ope_data_path})
ope_data_name=${ope_data%%.*}

if [ -z "$reftemp_data_path" ]; then
    echo " "
    echo "++++ Reference Template not specified."
    echo "--------------------------------------------------------------------------------"
else
    reftemp_data=$(basename ${reftemp_data_path})
    reftemp_data_name=${reftemp_data%%.*}
    echo " "
    echo "++++ Reference Template found: $reftemp_data_name"
    echo "--------------------------------------------------------------------------------"
fi

################################################################################
#
# REPORT INPUT PARAMETERS
#
################################################################################

reportParameters

################################################################################
# CREATE TEMPORARY DIRECTORIES

#data_prefix=sub${subject_id}_${modality}_run${run_number}
data_split=$(dirname $func_data_path)/${func_data_name}_split
data_moco=$(dirname $func_data_path)/${func_data_name}_moco
data_mats=$(dirname $func_data_path)/${func_data_name}_mats

if [ -d "$data_split" ]; then
    echo "-----> $data_split exists."
elif [ -d "$data_moco" ]; then
    echo "-----> $data_moco exists."
elif [ -d "$data_mats" ]; then
    echo "-----> $data_mats exists."
else
    mkdir $data_split
    mkdir $data_moco
    mkdir $data_mats
fi

#ope_prefix=sub${subject_id}_${modality}_ope_run${run_number}
ope_split=$(dirname $ope_data_path)/${ope_data_name}_split
ope_moco=$(dirname $ope_data_path)/${ope_data_name}_moco_out
ope_mats=$(dirname $ope_data_path)/${ope_data_name}_mats

if [ -d "$ope_split" ]; then
    echo "-----> $ope_split exists."
elif [ -d "$ope_moco" ]; then
    echo "-----> $ope_moco exists."
elif [ -d "$ope_mats" ]; then
    echo "-----> $ope_mats exists."
else
    mkdir $ope_split
    mkdir $ope_moco
    mkdir $ope_mats

    echo "-----> Temporary directories were created."
    echo " "
fi

################################################################################
# GET NUMBER OF VOLUMES FROM HEADER

n_vols_func=$($ANTSPATH/PrintHeader $(dirname $func_data_path)/${func_data_name}.nii.gz | grep Dimens | cut -d ',' -f 4 | cut -d ']' -f 1)

n_vols_ope=$($ANTSPATH/PrintHeader $(dirname $ope_data_path)/${ope_data_name}.nii.gz | grep Dimens | cut -d ',' -f 4 | cut -d ']' -f 1)

################################################################################
# DISASSEMBLE OPPOSITE PHASE-ENCODING DATA

$ANTSPATH/ImageMath \
    4 \
    ${ope_split}/${ope_data_name}_.nii.gz \
    TimeSeriesDisassemble \
    $(dirname $ope_data_path)/${ope_data_name}.nii.gz

echo "-----> Opposite phase data was disassembled into its ${bold} $n_vols_ope ${normal} constituent volumes."
echo " "

################################################################################
# START MOTION CORRECTION ON OPPOSITE PHASE-ENCODING DATA

basevol=1000 # ANTs indexing
fromvol=$(($basevol + 1))
nthvol=$(($basevol + $n_vols_ope - 1)) # Zero indexing

if [ -f "$(dirname $ope_data_path)/${ope_data_name}_MoCorr_meanTemplate.nii.gz" ]; then
    echo "-----> Motion corrected opposite phase-encoding template found."
else

    echo "-----> Starting motion correction on opposite phase-encoding data."
    echo " "

    start_time0=$(date +%s)

    for volume in $(eval echo "{$fromvol..$nthvol}"); do

        echo -ne "--------------------> Registering $volume"\\r

        FIXED=${ope_split}/${ope_data_name}_${basevol}.nii.gz
        MOVING=${ope_split}/${ope_data_name}_${volume}.nii.gz
        OUTPUT=${ope_moco}/${ope_data_name}_${volume}

        $ANTSPATH/antsRegistration \
            --verbose 0 \
            --float 1 \
            --dimensionality 3 \
            --use-histogram-matching 1 \
            --interpolation BSpline[4] \
            --collapse-output-transforms 1 \
            --output [ ${OUTPUT}_ , ${OUTPUT}_Warped.nii.gz , 1 ] \
            --winsorize-image-intensities [ 0.005 , 0.995 ] \
            --initial-moving-transform [ $FIXED , $MOVING , 1 ] \
            --transform Rigid[0.1] \
            --metric MI[ $FIXED , $MOVING , 1 , 64 , Regular , 0.25] \
            --convergence [ 500x250 , 1e-6 , 10 ] \
            --shrink-factors 2x1 \
            --smoothing-sigmas 1x0vox
        #            --transform Affine[0.1] \
        #            --metric MI[ $FIXED , $MOVING , 1 , 64 , Regular , 0.25] \
        #            --convergence [ 50 , 1e-6 , 10 ] \
        #            --shrink-factors 1 \
        #            --smoothing-sigmas 0vox

    done

    end_time0=$(date +%s)
    nettime0=$(expr $end_time0 - $start_time0)

    echo " "
    echo "-----> Motion estimation took in $(($nettime0 / 3600))h:$(($nettime0 % 3600 / 60))m:$(($nettime0 % 60))s."
    echo " "

    $ANTSPATH/ImageMath \
        4 \
        $(dirname $ope_data_path)/${ope_data_name}_MoCorr.nii.gz \
        TimeSeriesAssemble \
        $tr \
        0 \
        ${ope_split}/${ope_data_name}_${basevol}.nii.gz $ope_moco/*_Warped.nii.gz

    $ANTSPATH/AverageImages \
        3 \
        $(dirname $ope_data_path)/${ope_data_name}_MoCorr_meanTemplate.nii.gz \
        0 \
        ${ope_split}/${ope_data_name}_${basevol}.nii.gz \
        $ope_moco/${ope_data_name}_*_Warped.nii.gz &>/dev/null

    echo "-----> Created motion corrected mean template."
    echo " "

    rm -rf $ope_moco
    rm -rf $ope_mats
    rm -rf $ope_split
    rm $(dirname $ope_data_path)/${ope_data_name}_MoCorr.nii.gz

    echo "-----> Temporary files cleaned up."
    echo " "
fi

################################################################################
# DISASSEMBLE 4D FUNCTIONAL DATA
basevol=1000                            # ANTs indexing
nthvol=$(($basevol + $n_vols_func - 1)) # Zero indexing

if [ -f "${data_split}/${func_data_name}_${nthvol}.nii.gz" ]; then
    echo "-----> Disassembled timeseries data exists."
else
    $ANTSPATH/ImageMath \
        4 \
        ${data_split}/${func_data_name}_.nii.gz \
        TimeSeriesDisassemble \
        $(dirname $func_data_path)/${func_data_name}.nii.gz
fi

echo "-----> 4D data was disassembled into its ${bold} $n_vols_func ${normal} constituent volumes."
echo " "

################################################################################
# CHECK IF FIXED DATA WAS GIVEN, IF NOT MAKE ONE.

if [ -z "$fixed_data_path" ]; then

    ################################################################################
    # START MOTION CORRECTION ON FIRST 5 VOLS

    basevol=1000 # ANTs indexing
    fromvol=$(($basevol + 1))
    nthvol=$(($basevol + 4)) # Zero indexing

    echo "-----> Making fixed image."
    echo " "

    start_time0=$(date +%s)

    for volume in $(eval echo "{$fromvol..$nthvol}"); do

        echo -ne "--------------------> Registering $volume"\\r

        FIXED=${data_split}/${func_data_name}_${basevol}.nii.gz
        MOVING=${data_split}/${func_data_name}_${volume}.nii.gz
        OUTPUT=${data_moco}/${func_data_name}_${volume}

        $ANTSPATH/antsRegistration \
            --verbose 0 \
            --float 1 \
            --dimensionality 3 \
            --use-histogram-matching 1 \
            --interpolation $interpolation_type \
            --collapse-output-transforms 1 \
            --output [ ${OUTPUT}_ , ${OUTPUT}_Warped.nii.gz , 1 ] \
            --winsorize-image-intensities [ 0.005 , 0.995 ] \
            --initial-moving-transform [ $FIXED , $MOVING , 1 ] \
            --transform Rigid[0.1] \
            --metric MI[ $FIXED , $MOVING , 1 , 64 , Regular , 0.25] \
            --convergence [ 500x250 , 1e-6 , 10 ] \
            --shrink-factors 2x1 \
            --smoothing-sigmas 1x0vox \
            --transform Affine[0.1] \
            --metric MI[ $FIXED , $MOVING , 1 , 64 , Regular , 0.25] \
            --convergence [ 50 , 1e-6 , 10 ] \
            --shrink-factors 1 \
            --smoothing-sigmas 0vox

    done

    end_time0=$(date +%s)
    nettime0=$(expr $end_time0 - $start_time0)

    echo " "
    echo "-----> Motion estimation took in $(($nettime0 / 3600))h:$(($nettime0 % 3600 / 60))m:$(($nettime0 % 60))s."
    echo " "

    $ANTSPATH/ImageMath 4 \
        $(dirname $func_data_path)/${func_data_name}_MoCorr.nii.gz \
        TimeSeriesAssemble \
        $tr \
        0 \
        ${data_split}/${func_data_name}_${basevol}.nii.gz $data_moco/*_Warped.nii.gz

    $ANTSPATH/AverageImages \
        3 \
        $(dirname $func_data_path)/${func_data_name}_fixed.nii.gz \
        0 \
        ${data_split}/${func_data_name}_${basevol}.nii.gz \
        $data_moco/${func_data_name}_1001_Warped.nii.gz \
        $data_moco/${func_data_name}_1002_Warped.nii.gz \
        $data_moco/${func_data_name}_1003_Warped.nii.gz \
        $data_moco/${func_data_name}_1004_Warped.nii.gz &>/dev/null

    echo "-----> Created motion corrected mean template."
    echo " "

fi

################################################################################
# CHECK IF FIXED MASK WAS GIVEN, IF NOT MAKE ONE.

if [ -z "$fixed_data_mask" ]; then
    echo "-----> Making automatic fixed mask. Maybe sub-optimal, but its fine."
    echo ""
    ## Run N4 many times
    for i in {1..3}; do
        if [ "$i" -eq "1" ]; then
            ImageMath 3 \
                $(dirname $func_data_path)/${func_data_name}_fixed_Trunc.nii.gz \
                TruncateImageIntensity \
                $(dirname $func_data_path)/${func_data_name}_fixed.nii.gz \
                0.005 0.995

            N4BiasFieldCorrection \
                --verbose 0 \
                --image-dimensionality 3 \
                --shrink-factor 4 \
                --rescale-intensities 1 \
                --bspline-fitting [200] \
                --convergence [500x200x100x50,1e-9] \
                --input-image $(dirname $func_data_path)/${func_data_name}_fixed_Trunc.nii.gz \
                --output $(dirname $func_data_path)/${func_data_name}_fixed_N4.nii.gz
        else
            N4BiasFieldCorrection \
                --verbose 0 \
                --image-dimensionality 3 \
                --shrink-factor 4 \
                --rescale-intensities 1 \
                --bspline-fitting [200] \
                --convergence [50x50x50x50,1e-9] \
                --input-image $(dirname $func_data_path)/${func_data_name}_fixed_N4.nii.gz \
                --output $(dirname $func_data_path)/${func_data_name}_fixed_N4.nii.gz
        fi
    done

    ## DENOISE BEFORE MASKING
    DenoiseImage \
        --image-dimensionality 3 \
        --input-image $(dirname $func_data_path)/${func_data_name}_fixed_N4.nii.gz \
        --noise-model Rician \
        --output $(dirname $func_data_path)/${func_data_name}_fixed_N4_DEN.nii.gz

    $ANTSPATH/ImageMath \
        3 \
        $(dirname $func_data_path)/${func_data_name}_fixedMask.nii.gz \
        ThresholdAtMean \
        $(dirname $func_data_path)/${func_data_name}_fixed_N4_DEN.nii.gz

    ## MORPHOLOGICAL EROSION TO REDUCE EDGE EFFECTS
    $ANTSPATH/ImageMath \
        3 \
        $(dirname $func_data_path)/${func_data_name}_fixedMask.nii.gz \
        ME \
        $(dirname $func_data_path)/${func_data_name}_fixedMask.nii.gz 2

    ## GET THE LARGEST COMPONENT
    $ANTSPATH/ImageMath \
        3 \
        $(dirname $func_data_path)/${func_data_name}_fixedMask.nii.gz \
        GetLargestComponent \
        $(dirname $func_data_path)/${func_data_name}_fixedMask.nii.gz

    ## FILL HOLES
    $ANTSPATH/ImageMath \
        3 \
        $(dirname $func_data_path)/${func_data_name}_fixedMask.nii.gz \
        FillHoles \
        $(dirname $func_data_path)/${func_data_name}_fixedMask.nii.gz

    ## MORPHOLOGICAL DILATION
    $ANTSPATH/ImageMath \
        3 \
        $(dirname $func_data_path)/${func_data_name}_fixedMask.nii.gz \
        MD \
        $(dirname $func_data_path)/${func_data_name}_fixedMask.nii.gz 2

    ## MORPHOLOGICAL EROSION
    $ANTSPATH/ImageMath \
        3 \
        $(dirname $func_data_path)/${func_data_name}_fixedMask.nii.gz \
        ME \
        $(dirname $func_data_path)/${func_data_name}_fixedMask.nii.gz 1

    $ANTSPATH/ImageMath \
        3 \
        $(dirname $func_data_path)/${func_data_name}_fixedMask.nii.gz \
        MD \
        $(dirname $func_data_path)/${func_data_name}_fixedMask.nii.gz 2

    # Remove unwanted data
    rm $(dirname $func_data_path)/${func_data_name}_fixed_Trunc.nii.gz \
        $(dirname $func_data_path)/${func_data_name}_fixed_N4.nii.gz \
        $(dirname $func_data_path)/${func_data_name}_fixed_N4_DEN.nii.gz

fi

if [ -z "$staged_proc" ]; then

    ################################################################################
    # START MOTION COMPENSATION ON FUNCTIONAL DATA

    basevol=1000 # ANTs indexing

    nthvol=$(($basevol + $n_vols_func - 1)) # Zero indexing

    if [ -f "$(dirname $func_data_path)/${func_data_name}_MoCorr_meanTemplate.nii.gz" ]; then
        echo "-----> Motion corrected functional template found."
    else

        echo "-----> Starting motion correction on functional data."
        echo " "

        start_time0=$(date +%s)

        for volume in $(eval echo "{$basevol..$nthvol}"); do

            echo -ne "--------------------> Registering $volume"\\r

            if [ -z "$fixed_data_path" ]; then
                FIXED=$(dirname $func_data_path)/${func_data_name}_fixed.nii.gz
            fi

            FIXED=$(dirname $func_data_path)/${fixed_data_name}.nii.gz

            MOVING=${data_split}/${func_data_name}_${volume}.nii.gz

            if [ -z "$fixed_data_mask" ]; then
                MASK=$(dirname $func_data_path)/${func_data_name}_fixedMask.nii.gz
            fi

            MASK=$(dirname $func_data_path)/${fixed_mask_data_name}.nii.gz

            OUTPUT=${data_moco}/${func_data_name}_${volume}

            $ANTSPATH/antsRegistration \
                --verbose 0 \
                --float 1 \
                --use-estimate-learning-rate-once 1 \
                --dimensionality 3 \
                --random-seed 42 \
                --use-histogram-matching 1 \
                --interpolation $interpolation_type \
                --collapse-output-transforms 1 \
                --output [ ${OUTPUT}_ , ${OUTPUT}_Warped.nii.gz , 1 ] \
                --winsorize-image-intensities [ 0.005 , 0.995 ] \
                --masks [ $MASK , 1 ] \
                --initial-moving-transform [ $FIXED , $MOVING , 1 ] \
                --transform Rigid[0.1] \
                --metric MI[ $FIXED , $MOVING , 1 , 64 , Regular , 0.25] \
                --convergence [ 500x250 , 1e-6 , 10 ] \
                --shrink-factors 2x1 \
                --smoothing-sigmas 1x0vox \
                --transform Affine[0.1] \
                --metric MI[ $FIXED , $MOVING , 1 , 64 , Regular , 0.25 ] \
                --convergence [ 50 , 1e-6 , 10 ] \
                --shrink-factors 1 \
                --smoothing-sigmas 0vox

            $ANTSPATH/ConvertTransformFile 3 \
                ${OUTPUT}_0GenericAffine.mat \
                ${OUTPUT}_ants2itk.mat \
                --hm \
                --ras

            c3d_affine_tool \
                -ref $FIXED \
                -src $MOVING \
                ${OUTPUT}_ants2itk.mat \
                -ras2fsl \
                -o ${OUTPUT}_itk2fsl.mat

        done

        echo "-----> Compiling motion parameters file in mm and deg."
        echo " "

        for volume in $(eval echo "{$basevol..$nthvol}"); do

            matrix=${data_moco}/${func_data_name}_${volume}_itk2fsl.mat
            # Use 'avscale' to create file containing Translations (in mm) and Rotations (in deg)
            mm=$(${FSLDIR}/bin/avscale --allparams $matrix $FIXED | grep "Translations" | awk '{print $5 " " $6 " " $7}')
            mmx=$(echo $mm | cut -d " " -f 1)
            mmy=$(echo $mm | cut -d " " -f 2)
            mmz=$(echo $mm | cut -d " " -f 3)
            radians=$(${FSLDIR}/bin/avscale --allparams $matrix $FIXED | grep "Rotation Angles" | awk '{print $6 " " $7 " " $8}')
            radx=$(echo $radians | cut -d " " -f 1)
            degx=$(echo "$radx * (180 / 3.14159)" | sed 's/[eE]+\?/*10^/g' | bc -l)
            rady=$(echo $radians | cut -d " " -f 2)
            degy=$(echo "$rady * (180 / 3.14159)" | sed 's/[eE]+\?/*10^/g' | bc -l)
            radz=$(echo $radians | cut -d " " -f 3)
            degz=$(echo "$radz * (180 / 3.14159)" | sed 's/[eE]+\?/*10^/g' | bc -l)
            # The "%.6f" formatting specifier allows the numeric value to be as wide as it needs to be to accomodate the number
            # Then we mandate (include) a single space as a delimiter between values.
            #echo $(printf "%.6f" $mmx) $(printf "%.6f" $mmy) $(printf "%.6f" $mmz) $(printf "%.6f" $degx) $(printf "%.6f" $degy) $(printf "%.6f" $degz) >>$(dirname $func_data_path)/${func_data_name}_MoCorr.params
            echo $mmx $mmy $mmz $degx $degy $degz >>$(dirname $func_data_path)/${func_data_name}_MoCorr.params

            # Absolute RMS
            if [[ $volume -eq $nthvol ]]; then
                matrix1=${data_moco}/${func_data_name}_1000_itk2fsl.mat
                matrix2=${data_moco}/${func_data_name}_${volume}_itk2fsl.mat
                ${FSLDIR}/bin/rmsdiff $matrix1 $matrix2 $FIXED >>$(dirname $func_data_path)/${func_data_name}_MoCorr.rmsabs
            else
                matrix1=${data_moco}/${func_data_name}_1000_itk2fsl.mat
                matrix2=${data_moco}/${func_data_name}_$((${volume} + 1))_itk2fsl.mat
                ${FSLDIR}/bin/rmsdiff $matrix1 $matrix2 $FIXED >>$(dirname $func_data_path)/${func_data_name}_MoCorr.rmsabs
            fi
            # Relative RMS
            if [[ $volume -eq $basevol ]]; then
                matrix1=${data_moco}/${func_data_name}_${volume}_itk2fsl.mat
                matrix2=${data_moco}/${func_data_name}_${volume}_itk2fsl.mat
                ${FSLDIR}/bin/rmsdiff $matrix1 $matrix2 $FIXED >>$(dirname $func_data_path)/${func_data_name}_MoCorr.rmsrel
            else
                matrix1=${data_moco}/${func_data_name}_$((${volume} - 1))_itk2fsl.mat
                matrix2=${data_moco}/${func_data_name}_${volume}_itk2fsl.mat
                ${FSLDIR}/bin/rmsdiff $matrix1 $matrix2 $FIXED >>$(dirname $func_data_path)/${func_data_name}_MoCorr.rmsrel
            fi
        done

        echo "-----> Generating plots."
        echo " "

        ${FSLDIR}/bin/fsl_tsplot -i $(dirname $func_data_path)/${func_data_name}_MoCorr.params -t 'Translations (mm)' -u 1 --start=1 --finish=3 -a x,y,z -w 640 -h 144 -o $(dirname $func_data_path)/${func_data_name}_MoCorr_translations.png
        ${FSLDIR}/bin/fsl_tsplot -i $(dirname $func_data_path)/${func_data_name}_MoCorr.params -t 'Rotations (deg)' -u 1 --start=4 --finish=6 -a x,y,z -w 640 -h 144 -o $(dirname $func_data_path)/${func_data_name}_MoCorr_rotations.png
        ${FSLDIR}/bin/fsl_tsplot -i $(dirname $func_data_path)/${func_data_name}_MoCorr.rmsabs,$(dirname $func_data_path)/${func_data_name}_MoCorr.rmsrel, -t 'Mean Displacements (mm)' -u 1 -a absolute,relative -w 640 -h 144 -o $(dirname $func_data_path)/${func_data_name}_MoCorr_rms.png

        end_time0=$(date +%s)
        nettime0=$(expr $end_time0 - $start_time0)

        echo " "
        echo "-----> Motion estimation took in $(($nettime0 / 3600))h:$(($nettime0 % 3600 / 60))m:$(($nettime0 % 60))s."
        echo " "

        $ANTSPATH/ImageMath \
            4 \
            $(dirname $func_data_path)/${func_data_name}_MoCorr_nativeSpace.nii.gz \
            TimeSeriesAssemble \
            $tr \
            0 \
            $data_moco/*_Warped.nii.gz

        $ANTSPATH/AverageImages \
            3 \
            $(dirname $func_data_path)/${func_data_name}_MoCorr_meanTemplate.nii.gz \
            0 \
            $data_moco/*_Warped.nii.gz

        echo "-----> Created motion corrected mean template for distortion correction."
        echo " "

        cp -r ${data_moco}/*.mat ${data_mats}/
        rm -rf $data_moco
        # rm $(dirname $func_data_path)/${func_data_name}_MoCorr.nii.gz

        echo "-----> Temporary files cleaned up."
        echo " "
    fi

    ################################################################################
    # REGISTER THE OPPOSITE PE TEMPLATE TO THE FUNCTIONAL TEMPLATE - RIGID ONLY
    FIXED=$(dirname $func_data_path)/${func_data_name}_MoCorr_meanTemplate.nii.gz
    MOVING=$(dirname $ope_data_path)/${ope_data_name}_MoCorr_meanTemplate.nii.gz
    OUTPUT=$(dirname $ope_data_path)/${ope_data_name}_MoCorr_meanTemplate_reg2func

    $ANTSPATH/antsRegistration \
        --verbose 0 \
        --use-estimate-learning-rate-once 1 \
        --float 1 \
        --random-seed 42 \
        --dimensionality 3 \
        --use-histogram-matching 1 \
        --interpolation $interpolation_type \
        --collapse-output-transforms 1 \
        --output [ ${OUTPUT}_ , ${OUTPUT}.nii.gz , 1 ] \
        --winsorize-image-intensities [ 0.005 , 0.995 ] \
        --initial-moving-transform [ $FIXED , $MOVING , 1 ] \
        --transform Rigid[0.1] \
        --metric MI[ $FIXED , $MOVING , 1 , 64, Regular, 0.4 ] \
        --convergence [ 100 , 1e-6 , 10 ] \
        --shrink-factors 1 \
        --smoothing-sigmas 0vox

    echo "-----> Registered Opposite PE template to the functional template."
    echo " "

    ################################################################################
    # RUN ANTS TEMPLATE CONSTRUCTION TO FIND MIDDLE GROUND
    FUNC_TEMP=$(dirname $func_data_path)/${func_data_name}_MoCorr_meanTemplate
    OPE_TEMP=$(dirname $ope_data_path)/${ope_data_name}_MoCorr_meanTemplate_reg2func
    OUTPUT=$(dirname $func_data_path)/${func_data_name}_DistCorr_

    if [ -z "$reftemp_data_name" ]; then
        echo "-----> Starting distortion correction."
        echo " "
        #    $ANTSPATH/sk_antsDistCorrScript.sh \
        $ANTSPATH/antsMultivariateTemplateConstruction2.sh \
            -d 3 \
            -a 0 \
            -c 2 \
            -g 0.1 \
            -j $n_threads \
            -n 1 \
            -l 1 \
            -r 1 \
            -i $iterations \
            -k 1 \
            -f 4x2x1 \
            -s 2x1x0vox \
            -q 50x25x15 \
            -t SyN \
            -m CC \
            -o $OUTPUT \
            ${FUNC_TEMP}.nii.gz ${OPE_TEMP}.nii.gz &>/dev/null
    else
        echo "-----> Starting distortion correction. Also using reference template."
        echo " "
        #    $ANTSPATH/sk_antsDistCorrScript.sh \
        $ANTSPATH/antsMultivariateTemplateConstruction2.sh \
            -d 3 \
            -a 0 \
            -c 2 \
            -g 0.1 \
            -j $n_threads \
            -n 1 \
            -l 1 \
            -r 1 \
            -i $iterations \
            -k 1 \
            -f 4x2x1 \
            -s 2x1x0vox \
            -q 50x25x15 \
            -z $reftemp_data_path \
            -t SyN \
            -m CC \
            -o $OUTPUT \
            ${FUNC_TEMP}.nii.gz ${OPE_TEMP}.nii.gz $reftemp_data_path &>/dev/null
    fi

    echo "-----> Created distortion corrected template."
    echo " "

    $ANTSPATH/ImageMath 3 \
        $(dirname $func_data_path)/$(basename ${OUTPUT})template0.nii.gz \
        UnsharpMask \
        $(dirname $func_data_path)/$(basename ${OUTPUT})template0.nii.gz \
        0.7 2 0 0

    # Do clean up
    rm \
        $(dirname $func_data_path)/$(basename ${OUTPUT})$(basename ${OPE_TEMP})10GenericAffine.mat \
        $(dirname $func_data_path)/$(basename ${OUTPUT})$(basename ${OPE_TEMP})11InverseWarp.nii.gz \
        $(dirname $func_data_path)/$(basename ${OUTPUT})$(basename ${OPE_TEMP})11Warp.nii.gz \
        $(dirname $func_data_path)/$(basename ${OUTPUT})template0GenericAffine.mat \
        $(dirname $func_data_path)/$(basename ${OUTPUT})template0warp.nii.gz \
        $(dirname $func_data_path)/$(basename ${OUTPUT})templatewarplog.txt \
        $(dirname $func_data_path)/$(basename ${OUTPUT})template0$(basename ${OPE_TEMP})1WarpedToTemplate.nii.gz \
        $(dirname $func_data_path)/$(basename ${OUTPUT})template0$(basename ${FUNC_TEMP})0WarpedToTemplate.nii.gz \
        $(dirname $func_data_path)/$(basename ${OUTPUT})template0$(basename ${OPE_TEMP})Repaired.nii.gz \
        $(dirname $func_data_path)/$(basename ${OUTPUT})template0$(basename ${FUNC_TEMP})Repaired.nii.gz

    mv \
        $(dirname $func_data_path)/$(basename ${OUTPUT})$(basename ${FUNC_TEMP})01Warp.nii.gz \
        $(dirname $func_data_path)/$(basename ${OUTPUT})01Warp.nii.gz

    mv \
        $(dirname $func_data_path)/$(basename ${OUTPUT})$(basename ${FUNC_TEMP})01InverseWarp.nii.gz \
        $(dirname $func_data_path)/$(basename ${OUTPUT})01InverseWarp.nii.gz

    mv \
        $(dirname $func_data_path)/$(basename ${OUTPUT})$(basename ${FUNC_TEMP})00GenericAffine.mat \
        $(dirname $func_data_path)/$(basename ${OUTPUT})00GenericAffine.mat

    echo "-----> Temporary files cleaned up."
    echo " "

    tar -zcf \
        $(dirname $func_data_path)/${func_data_name}_MoCorr_DistCorr_transforms.tar.gz \
        ${data_mats}/*.mat \
        $(dirname $func_data_path)/$(basename ${OUTPUT})00GenericAffine.mat \
        $(dirname $func_data_path)/$(basename ${OUTPUT})01Warp.nii.gz \
        $(dirname $func_data_path)/$(basename ${OUTPUT})01InverseWarp.nii.gz

    echo "-----> Transformations saved."
    echo " "

    end_time0=$(date +%s)
    nettime0=$(expr $end_time0 - $start_time0)

    echo " "
    echo "-----> Completed combined estimation in $(($nettime0 / 3600))h:$(($nettime0 % 3600 / 60))m:$(($nettime0 % 60))s."
    echo " "
else
    echo "-----> Staged processing enabled, so only made the template images, mask and quick undistortion. Need to re-run this script without -s "
    echo " "

    # ################################################################################
    # # RUN ANTS TEMPLATE CONSTRUCTION TO FIND MIDDLE GROUND
    # FUNC_TEMP=$(dirname $func_data_path)/${func_data_name}_MoCorr_fixed
    # OPE_TEMP=$(dirname $ope_data_path)/${ope_data_name}_MoCorr_meanTemplate_reg2func
    # OUTPUT=$(dirname $func_data_path)/${func_data_name}_DistCorr_quick

    # if [ -z "$reftemp_data_name" ]; then
    #     echo "-----> Starting distortion correction."
    #     echo " "
    #     #    $ANTSPATH/sk_antsDistCorrScript.sh \
    #     $ANTSPATH/antsMultivariateTemplateConstruction2.sh \
    #     -d 3 \
    #     -a 0 \
    #     -c 2 \
    #     -g 0.25 \
    #     -j $n_threads \
    #     -n 1 \
    #     -l 1 \
    #     -r 0 \
    #     -i 3 \
    #     -k 1 \
    #     -f 4x2x1 \
    #     -s 4x2x0vox \
    #     -q 15x10x10 \
    #     -t SyN \
    #     -m CC[2] \
    #     -o $OUTPUT \
    #     ${FUNC_TEMP}.nii.gz ${OPE_TEMP}.nii.gz 2>&1 | tee $(dirname $func_data_path)/${func_data_name}_DistCorr.log
    # else
    #     echo "-----> Starting distortion correction. Also using reference template."
    #     echo " "
    #     #    $ANTSPATH/sk_antsDistCorrScript.sh \
    #     $ANTSPATH/antsMultivariateTemplateConstruction2.sh \
    #     -d 3 \
    #     -a 0 \
    #     -c 2 \
    #     -g 0.25 \
    #     -j $n_threads \
    #     -n 1 \
    #     -l 1 \
    #     -r 0 \
    #     -i 3 \
    #     -k 1 \
    #     -f 4x2x1 \
    #     -s 4x2x0vox \
    #     -q 15x10x10 \
    #     -z $reftemp_data_path \
    #     -t SyN \
    #     -m CC[2] \
    #     -o $OUTPUT \
    #     ${FUNC_TEMP}.nii.gz ${OPE_TEMP}.nii.gz $reftemp_data_path 2>&1 | tee $(dirname $func_data_path)/${func_data_name}_DistCorr.log
    # fi

    # echo "-----> Created distortion corrected template."
    # echo " "

    # $ANTSPATH/ImageMath 3 \
    # $(dirname $func_data_path)/$(basename ${OUTPUT})template0.nii.gz \
    # UnsharpMask \
    # $(dirname $func_data_path)/$(basename ${OUTPUT})template0.nii.gz \
    # 0.7 2 0 0

    # # Do clean up
    # rm \
    # $(dirname $func_data_path)/$(basename ${OUTPUT})$(basename ${OPE_TEMP})10GenericAffine.mat \
    # $(dirname $func_data_path)/$(basename ${OUTPUT})$(basename ${OPE_TEMP})11InverseWarp.nii.gz \
    # $(dirname $func_data_path)/$(basename ${OUTPUT})$(basename ${OPE_TEMP})11Warp.nii.gz \
    # $(dirname $func_data_path)/$(basename ${OUTPUT})$(basename ${FUNC_TEMP})01InverseWarp.nii.gz \
    # $(dirname $func_data_path)/$(basename ${OUTPUT})template0GenericAffine.mat \
    # $(dirname $func_data_path)/$(basename ${OUTPUT})template0warp.nii.gz \
    # $(dirname $func_data_path)/$(basename ${OUTPUT})templatewarplog.txt \
    # $(dirname $func_data_path)/$(basename ${OUTPUT})template0$(basename ${OPE_TEMP})1WarpedToTemplate.nii.gz \
    # $(dirname $func_data_path)/$(basename ${OUTPUT})template0$(basename ${FUNC_TEMP})0WarpedToTemplate.nii.gz

    # mv \
    # $(dirname $func_data_path)/$(basename ${OUTPUT})$(basename ${FUNC_TEMP})01Warp.nii.gz \
    # $(dirname $func_data_path)/$(basename ${OUTPUT})01Warp.nii.gz

    # mv \
    # $(dirname $func_data_path)/$(basename ${OUTPUT})$(basename ${FUNC_TEMP})00GenericAffine.mat \
    # $(dirname $func_data_path)/$(basename ${OUTPUT})00GenericAffine.mat

    # echo "-----> Temporary files cleaned up."
    # echo " "

    # tar -zcf \
    # $(dirname $func_data_path)/${func_data_name}_MoCorr_DistCorr_transforms.tar.gz \
    # ${data_mats}/*.mat \
    # $(dirname $func_data_path)/$(basename ${OUTPUT})00GenericAffine.mat \
    # $(dirname $func_data_path)/$(basename ${OUTPUT})01Warp.nii.gz

    # echo "-----> Transformations saved."
    # echo " "

    end_time0=$(date +%s)
    nettime0=$(expr $end_time0 - $start_time0)

    echo " "
    echo "-----> Completed combined estimation in $(($nettime0 / 3600))h:$(($nettime0 % 3600 / 60))m:$(($nettime0 % 60))s."
    echo " "
fi
