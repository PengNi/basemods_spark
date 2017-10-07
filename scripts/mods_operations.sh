#!/usr/bin/env bash
# FIXME add the part of "makeMotifReport.sh" or not?


# get input params-----------------------------------
SEYMOUR_HOME=$1
temp_output_folder=$2
modsgff_filepath=$3
modscsv_filepath=$4
reference_filepath=$5
motifs_gff_gz=$6


# Setting up SMRTpipe environment
echo "Setting up ENV on $(uname -n) for methylation detection"
source $SEYMOUR_HOME/etc/setup.sh


# set variables--------------------------------------
xsl_filepath=$SEYMOUR_HOME/analysis/etc/xsl/brandedGraphReport.xsl

modsgff_filename=${modsgff_filepath##*/}
fileprefix=$(echo ${modsgff_filename%.modifications.gff})

motif_summary_csv_filename="${fileprefix}.motif_summary.csv"
motif_summary_xml_filename="${fileprefix}.motif_summary.xml"
motif_summary_html_filename="${fileprefix}.motif_summary.html"
motif_tmpPQUF7M_gff="${fileprefix}.tmpPQUF7M.gff"


# task for P_MotifFinder.findMotifs -----------------
echo "Task findMotifs started on $(date)"
# Task 1
java -Xmx8000m -jar ${SEYMOUR_HOME}/analysis/lib/java/motif-maker-0.2.one-jar.jar find --minScore 30 --gff $modsgff_filepath --fasta $reference_filepath --output $temp_output_folder/$motif_summary_csv_filename --xml $temp_output_folder/$motif_summary_xml_filename || exit $?
echo "Task 1 completed at $(date)"
# Task 2
saxonb9 -xsl:$xsl_filepath -s:$temp_output_folder/$motif_summary_xml_filename -o:$temp_output_folder/$motif_summary_html_filename || exit $?
echo "Task 2 completed at $(date)"

echo -e "Task findMotifs finished on $(date)\n\n"


# task for  P_MotifFinder.makeMotifGff ---------------------------------------
echo "Task makeMotifGff started on $(date)"
# Task 1
java -Xmx8000m -jar ${SEYMOUR_HOME}/analysis/lib/java/motif-maker-0.2.one-jar.jar reprocess -m $temp_output_folder/$motif_summary_csv_filename --gff $modsgff_filepath --fasta $reference_filepath --csv $modscsv_filepath --output $temp_output_folder/$motif_tmpPQUF7M_gff || exit $?
echo "Task 1 completed at $(date)"
# Task 2
gzip --no-name -c $temp_output_folder/$motif_tmpPQUF7M_gff > $temp_output_folder/$motifs_gff_gz || exit $?
echo "Task 2 completed at $(date)"
# Task 3
unlink $temp_output_folder/$motif_tmpPQUF7M_gff || exit $?
echo "Task 3 completed at $(date)"

echo -e "Task makeMotifGff finished on $(date)\n\n"



# exit-----------------------------------------------
rcode=$?
echo "modification operations finished on $(date) with exit code ${rcode}."
exit ${rcode}
