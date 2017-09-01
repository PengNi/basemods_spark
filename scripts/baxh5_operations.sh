#!/usr/bin/env bash
# FIXME don't forget to let users set nproc

# get input params-----------------------------------
SEYMOUR_HOME=$1
temp_output_folder=$2
baxh5_filepath=$3
reference_filepath=$4
referencesa_filepath=$5
blasr_cmph5file=$6
kernel_num=$7


# Setting up SMRTpipe environment
echo "Setting up ENV on $(uname -n) for methylation detection"
source $SEYMOUR_HOME/etc/setup.sh


# set variables--------------------------------------
baxh5filename=${baxh5_filepath##*/}
fileprefix=$(echo ${baxh5filename%.bax.h5})

baxh5_fofn_file="${fileprefix}.fofn"
if [ -e $temp_output_folder/$baxh5_fofn_file ]
then
rm $temp_output_folder/$baxh5_fofn_file
fi
echo $baxh5_filepath >> $temp_output_folder/$baxh5_fofn_file

filtered_summary_csvfile="${fileprefix}.filtered_summary.csv"
filterd_regions_fofnfile="${fileprefix}.filterd_regions.fofn"

blasrtmp_samfile="${fileprefix}.tmp.sam"
blasrtmp_filterd_samfile="${fileprefix}.tmp.filtered.sam"
# blasr_cmph5file="${fileprefix}.aligned_reads.cmp.h5"


# task for P_Filter.filter---------------------------
echo "Task filter started on $(date)"
filter_plsh5.py --debug --filter='MinReadScore=0.7500,MinSRL=50,MinRL=50' --trim='True' --outputDir=$temp_output_folder --outputSummary=$temp_output_folder/$filtered_summary_csvfile --outputFofn=$temp_output_folder/$filterd_regions_fofnfile $temp_output_folder/$baxh5_fofn_file || exit $?
echo -e "Task filter finished at $(date)\n\n"


# task for P_Mapping.align ---------------------------
echo "Task align started on $(date)"
# Task 1
blasr $temp_output_folder/$baxh5_fofn_file $reference_filepath -sam -out $temp_output_folder/$blasrtmp_samfile -regionTable $temp_output_folder/$filterd_regions_fofnfile -minSubreadLength 50 -minReadLength 50 -randomSeed 1 -concordant -useQuality -minMatch 12 -bestn 10 -minPctIdentity 70.0 -sa $referencesa_filepath -placeRepeatsRandomly -nproc $kernel_num
samFilter $temp_output_folder/$blasrtmp_samfile $reference_filepath $temp_output_folder/$blasrtmp_filterd_samfile -minAccuracy 75 -minLength 50 -seed 1 -hitPolicy randombest
samtoh5 $temp_output_folder/$blasrtmp_filterd_samfile $reference_filepath $temp_output_folder/$blasr_cmph5file -readType standard
rm $temp_output_folder/$blasrtmp_samfile
rm $temp_output_folder/$blasrtmp_filterd_samfile
echo "Task align: Task 1 completed at $(date)"
# Task 2
echo 'Alignment Complete' || exit $?
echo "Task align: Task 2 completed at $(date)"
## Task 3
## task 3 can be commented
#loadChemistry.py $temp_output_folder/$baxh5_fofn_file $temp_output_folder/$blasr_cmph5file || exit $?
#echo "Task align: Task 3 completed at $(date)"
# Task 4
loadPulses $temp_output_folder/$baxh5_fofn_file $temp_output_folder/$blasr_cmph5file -metrics DeletionQV,IPD,InsertionQV,PulseWidth,QualityValue,MergeQV,SubstitutionQV,DeletionTag -byread || exit $?
echo 'LoadPulses Complete' || exit $?
echo "Task align: Task 4 completed at $(date)"

echo -e "Task align finished on $(date)\n\n"



# exit-----------------------------------------------
rcode=$?
echo "baxh5 operations finished on $(date) with exit code ${rcode}."
exit ${rcode}