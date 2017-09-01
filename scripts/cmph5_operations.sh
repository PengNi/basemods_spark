#!/usr/bin/env bash
# FIXME don't forget to let users set nproc
# FIXME ipdsummary.py --identify add all methylation types

# Setting up SMRTpipe environment
echo "Setting up ENV on $(uname -n) for methylation detection"
SEYMOUR_HOME=/home/smrtanalysis/smrtana/install/smrtanalysis_2.3.0.140936
source $SEYMOUR_HOME/etc/setup.sh

# get input params-----------------------------------
temp_output_folder=$1
cmph5_filepath=$2
ref_chunk_info=$3
reference_filepath=$4

tmpcOnrEX_gff=$5
tmpcc5Wn6_csv=$6
kernel_num=$7

# set variables--------------------------------------
ipdSummaryParamsPath=$SEYMOUR_HOME/analysis/etc/algorithm_parameters/2014-09/kineticsTools

cmph5filename=${cmph5_filepath##*/}
fileprefix=$(echo ${cmph5filename%.cmp.h5})

cmph5_tmp_filename="${fileprefix}_TMP"
base_mod_contig_txt="${fileprefix}.base_mod_contig_ids.txt"
temp_kinetics_h5="${fileprefix}.temp_kinetics.h5"
#tmpcOnrEX_gff="${fileprefix}.tmpcOnrEX.gff"
#tmpcc5Wn6_csv="${fileprefix}.tmpcc5Wn6.csv"
#modifications_gff_gz="${fileprefix}.modifications.gff.gz"
#modifications_csv_gz="${fileprefix}.modifications.csv.gz"



# task for P_Mapping.sort ---------------------------
echo "Task sort started on $(date)"
# Task 1
cmph5tools.py -vv sort --deep --inPlace $cmph5_filepath || exit $?
echo "Task sort: Task 1 completed at $(date)"

echo -e "Task sort finished on $(date)\n\n"


## task for P_Mapping.repack --------------------------
#echo "Task repack started on $(date)"
## Task 1
#((which h5repack && (h5repack -f GZIP=1 $cmph5_filepath $temp_output_folder/$cmph5_tmp_filename && mv /$temp_output_folder/$cmph5_tmp_filename $cmph5_filepath)) || echo 'no h5repack found, continuing w/out') || exit $?
#echo "Task repack: Task 1 completed at $(date)"
#
#echo -e "Task repack finished on $(date)\n\n"


# task for P_ModificationDetection.writeContigList ----
echo "Task writeContigList started on $(date)"
# Task 1
rm -f $temp_output_folder/$base_mod_contig_txt || exit $?
echo "Task writeContigList: Task 1 completed at $(date)"
# Task 2
echo "ref000001" >> $temp_output_folder/$base_mod_contig_txt || exit $?
echo "Task writeContigList: Task 2 completed at $(date)"

echo -e "Task writeContigList finished on $(date)\n\n"


# task for P_ModificationDetection.computeModifications --
echo "Task computeModifications with nproc 39. Started on $(date)"
# Task 1
ipdSummary.py -v -W $temp_output_folder/$base_mod_contig_txt --methylFraction --identify m6A,m4C --paramsPath $ipdSummaryParamsPath --numWorkers $kernel_num --summary_h5 $temp_output_folder/$temp_kinetics_h5 --gff $temp_output_folder/$tmpcOnrEX_gff --csv $temp_output_folder/$tmpcc5Wn6_csv --reference $reference_filepath --refChunkInfo $ref_chunk_info $cmph5_filepath || exit $?
echo "Task computeModifications: Task 1 completed at $(date)"
## Task 2
#gzip --no-name -c $temp_output_folder/$tmpcOnrEX_gff > $temp_output_folder/$modifications_gff_gz || exit $?
#echo "Task computeModifications: Task 2 completed at $(date)"
## Task 3
#unlink $temp_output_folder/$tmpcOnrEX_gff || exit $?
#echo "Task computeModifications: Task 3 completed at $(date)"
## Task 4
#gzip --no-name -c $temp_output_folder/$tmpcc5Wn6_csv > $temp_output_folder/$modifications_csv_gz || exit $?
#echo "Task computeModifications: Task 4 completed at $(date)"
## Task 5
#unlink $temp_output_folder/$tmpcc5Wn6_csv || exit $?
#echo "Task computeModifications: Task 5 completed at $(date)"

echo -e "Task computeModifications finished on $(date)\n\n"



# exit-----------------------------------------------
rcode=$?
echo "cmph5 operations finished on $(date) with exit code ${rcode}."
exit ${rcode}










