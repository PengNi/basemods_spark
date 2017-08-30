#!/usr/bin/python
# coding=utf-8
from pyspark import SparkContext, SparkConf, SparkFiles
from subprocess import Popen, PIPE
import h5py
import shlex
import os
# import shutil
import numpy

from bash5_splitting import split_holenumbers, makeOffsetsDataStructure, get_movieName
from fasta_info import getRefInfoFromFastaFiles
from movie_chemistry import getMoviesChemistry
from CmpH5Format import CmpH5Format

# FIXME: the file path, especially the ref file path should be set by user.
# FIXME: better in a configure file?
TEMP_OUTPUT_FOLDER = "/tmp/basmod_spark_data"
SCRIPTS_FOLDER = "/home/hadoop/workspace_py/basemods_spark/scripts"
reference_directory = '/home/hadoop/workspace_py/basemods_spark/data/lambda/sequence'
ref_filename = 'lambda.fasta'
ref_sa_filename = 'lambda.fasta.sa'
metadataxml_dir = '/home/hadoop/workspace_py/basemods_spark/data/lambda_v210'
metadataxml_name = 'm130802_062611_ethan_c100542982550000001823084811241306_s1_p0.metadata.xml'

H5GROUP = h5py._hl.group.Group
H5DATASET = h5py._hl.dataset.Dataset
refMaxLength = 3e12
COLUMNS = 60
PAD = 15


# ---------------------------------------------------------------------------------
# FIXME: 1. how to convert a dataset to an rdd without reading the whole dataset into memory?
# FIXME: 2. 现在的做法是对于每一个dataset都转成一个RDD，然后再用union+groupByKey合并。
# FIXME:    根据每个key的信息，直接将所有该key的dataset读成rdd里的一个元素，省掉union和group，
# FIXME:    是不是更好？（可能的情况：这种做法会更快，但会更耗内存？）
# convert baxh5 to RDD------------------------------------------------
def baxh5toRDD(sc, baxh5file, folds=1, numpartitions=3):
    f = h5py.File(baxh5file, "r")
    holenumbers = f['/PulseData/BaseCalls/ZMW/HoleNumber'].value
    if folds > len(holenumbers):
        folds = len(holenumbers)
    hole_splitspots = split_holenumbers(holenumbers, folds)
    hole2range = makeOffsetsDataStructure(f)  # todo: how to make it faster
    basecall_splitspots = get_basecall_range_of_each_holesblock(hole_splitspots,
                                                                holenumbers, hole2range)
    moviename = get_movieName(f)
    # print(hole_splitspots)

    rdds = []

    # datasets in PulseData/BaseCalls/ZMW
    # FIXME: use (for key in f['']) or (for key in f[''].keys())?
    for key in f['/PulseData/BaseCalls/ZMW']:
        rdds.append(convert_dataset_in_zmw_to_rdd(sc, numpartitions, hole_splitspots,
                                                  holenumbers, f,
                                                  '/PulseData/BaseCalls/ZMW/' + str(key),
                                                  moviename))

    # datasets in PulseData/BaseCalls/ZMWMetrics
    for key in f['/PulseData/BaseCalls/ZMWMetrics']:
        rdds.append(convert_dataset_in_zmw_to_rdd(sc, numpartitions, hole_splitspots,
                                                  holenumbers, f,
                                                  '/PulseData/BaseCalls/ZMWMetrics/' + str(key),
                                                  moviename))

    # datasets in PulseData/BaseCalls
    for key in f['/PulseData/BaseCalls']:
        if not (str(key) == 'ZMWMetrics' or str(key) == 'ZMW'):
            rdds.append(convert_dataset_in_basecalls_to_rdd(sc, numpartitions, hole_splitspots,
                                                            holenumbers, basecall_splitspots, f,
                                                            '/PulseData/BaseCalls/' + str(key),
                                                            moviename))

    # PulseData/Regions
    rdds.append(convert_regions_dataset_to_rdd(sc, numpartitions, hole_splitspots, holenumbers, f,
                                               '/PulseData/Regions', moviename))

    # baxh5 attrs
    baxh5attrs = sc.broadcast(get_baxh5_attrs(f))

    f.close()
    if len(rdds) == 1:
        wholeinfo_rdd = rdds[0]
    elif len(rdds) > 1:
        wholeinfo_rdd = sc.union(rdds)\
            .groupByKey() \
            .map(lambda (x, y): (x, list(y)))\
            .map(lambda (x, y): (x, (y, baxh5attrs.value)))
    else:
        wholeinfo_rdd = None
        print("baxh5tordd wrong!")

    # print(wholeinfo_rdd.first())
    # -----------------------------------------------------------------
    return wholeinfo_rdd


def convert_dataset_in_zmw_to_rdd(sc, numpartitions, hole_splitspots, holenumbers,
                                  f, datasetname, moviename):
    sdata = f[datasetname].value
    dtype = sdata.dtype
    sshape = sdata.shape

    scparas = list()
    if len(sshape) == 1:
        for holess in hole_splitspots:
            scparas.append(((moviename, (holenumbers[holess[0]], holenumbers[holess[1]-1])),
                            ((datasetname, dtype), sdata[holess[0]:holess[1]])))
    else:
        for holess in hole_splitspots:
            scparas.append(((moviename, (holenumbers[holess[0]], holenumbers[holess[1]-1])),
                            ((datasetname, dtype), sdata[holess[0]:holess[1], :])))

    return sc.parallelize(scparas, numpartitions)


def convert_dataset_in_basecalls_to_rdd(sc, numpartitions, hole_splitspots, holenumbers,
                                        basecall_splitspots, f, datasetname, moviename):
    sdata = f[datasetname].value
    dtype = sdata.dtype
    sshape = sdata.shape

    scparas = list()
    if len(sshape) == 1:
        for i in range(0, len(hole_splitspots)):
            scparas.append(((moviename, (holenumbers[hole_splitspots[i][0]],
                                         holenumbers[hole_splitspots[i][1]-1])),
                            ((datasetname, dtype),
                             sdata[basecall_splitspots[i][0]:basecall_splitspots[i][1]])))
    else:
        for i in range(0, len(hole_splitspots)):
            scparas.append(((moviename, (holenumbers[hole_splitspots[i][0]],
                                         holenumbers[hole_splitspots[i][1]-1])),
                            ((datasetname, dtype),
                             sdata[basecall_splitspots[i][0]:basecall_splitspots[i][1], :])))
    return sc.parallelize(scparas, numpartitions)


def convert_regions_dataset_to_rdd(sc, numpartitions, hole_splitspots, holenumbers,
                                   f, datasetname, moviename):
    regiondata = f[datasetname].value
    dtype = regiondata.dtype
    sshape = regiondata.shape

    scparas = list()
    locs_start = 0
    for holess in hole_splitspots[:-1]:
        holenumbers_tmp = set(holenumbers[holess[0]:holess[1]])
        for i in range(locs_start, sshape[0]):
            if regiondata[i, 0] not in holenumbers_tmp:
                scparas.append(((moviename, (holenumbers[holess[0]], holenumbers[holess[1]-1])),
                                ((datasetname, dtype), regiondata[locs_start:i, :])))
                locs_start = i
                break
    scparas.append(((moviename, (holenumbers[hole_splitspots[-1][0]],
                                 holenumbers[hole_splitspots[-1][1]-1])),
                    ((datasetname, dtype), regiondata[locs_start:, :])))
    return sc.parallelize(scparas, numpartitions)


def get_basecall_range_of_each_holesblock(hole_splitspots, holenumbers, holerange):
    basecall_splitspots = []
    for hole_splitspot in hole_splitspots:
        begin, end = holerange[holenumbers[hole_splitspot[0]]][0], \
                     holerange[holenumbers[hole_splitspot[1]-1]][1]
        basecall_splitspots.append((begin, end))
    return basecall_splitspots


# FIXME: 1. how to get a list of all the datasets and groups of an h5file at once?
# FIXME: h5obj.visit(printname) can only print the list
# FIXME: 2. how to check a item in h5obj is a group or a dataset?
# get baxh5 attrs------------------------------------------------------
def get_baxh5_attrs(baxh5obj):
    attrslist = []
    # ScanData
    attrslist.append((('/ScanData', H5GROUP), get_h5item_attrs(baxh5obj, '/ScanData')))
    for key in baxh5obj['/ScanData'].keys():
        pathtmp = '/ScanData/' + str(key)
        if isinstance(baxh5obj[pathtmp], H5GROUP):
            attrslist.append(((pathtmp, H5GROUP), get_h5item_attrs(baxh5obj, pathtmp)))
        else:
            attrslist.append(((pathtmp, H5DATASET), get_h5item_attrs(baxh5obj, pathtmp)))
    # PulseData
    attrslist.append((('/PulseData/BaseCalls', H5GROUP),
                      get_h5item_attrs(baxh5obj, '/PulseData/BaseCalls')))
    attrslist.append((('/PulseData/Regions', H5DATASET),
                      get_h5item_attrs(baxh5obj, '/PulseData/Regions')))
    return attrslist


def get_h5item_attrs(h5obj, itempath):
    attrslist = []
    itemattrs = h5obj[itempath].attrs
    for key, val in itemattrs.items():
        attrslist.append((key, val))
    return attrslist


# operations for each rdd element in baxh5RDD-------------------------
def basemods_pipeline_baxh5_operations(keyval):
    """
    keyval: an element of baxh5RDD
    ((moviename, holerange), ([((datasetname, dtype), dataset_value),...],
    [((name, type(group or dataset)), [(attriname, attrival), ]),...]))
    :param keyval:
    :return: aligned reads
    """
    fileinfo, filecontent = keyval

    name_prefix = fileinfo[0] + "." + str(fileinfo[1][0]) + "-" + str(fileinfo[1][1])
    baxh5file = name_prefix + ".bax.h5"
    reference_path = SparkFiles.get(ref_filename)
    referencesa_path = SparkFiles.get(ref_sa_filename)
    # metadataxml_path = SparkFiles.get(metadataxml_name)
    cmph5file = name_prefix + ".aligned_reads.cmp.h5"

    if not os.path.isdir(TEMP_OUTPUT_FOLDER):
        os.mkdir(TEMP_OUTPUT_FOLDER, 0777)
    else:
        os.chmod(TEMP_OUTPUT_FOLDER, 0o777)

    # # copy metadata.xml from spark tmp dir to TEMP_OUTPUT_FOLDER
    # os.chmod(metadataxml_path, 0o777)
    # shutil.copy2(metadataxml_path, TEMP_OUTPUT_FOLDER + "/" + metadataxml_name)

    # write baxh5obj to file
    baxh5_dir = TEMP_OUTPUT_FOLDER + '/baxh5'
    if not os.path.isdir(baxh5_dir):
        os.mkdir(baxh5_dir)
    baxh5path = baxh5_dir + '/' + baxh5file
    writebaxh5(filecontent, baxh5path)

    # baxh5 operations (filter, align(blasr, filter, samtoh5), loadchemistry, loadpulse)
    baxh5_operations = "{scripts_fold}/baxh5_operations.sh {temp_output_folder} {baxh5_filepath}" \
                       " {reference_filepath} {referencesa_filepath} {cmph5_filename}".\
        format(scripts_fold=SCRIPTS_FOLDER,
               temp_output_folder=TEMP_OUTPUT_FOLDER,
               baxh5_filepath=baxh5path,
               reference_filepath=reference_path,
               referencesa_filepath=referencesa_path,
               cmph5_filename=cmph5file)
    baxh5_process = Popen(shlex.split(baxh5_operations), stdout=PIPE, stderr=PIPE)
    baxh5_out, baxh5_error = baxh5_process.communicate()

    if "[Errno" in baxh5_error.strip() or "error" in baxh5_error.strip().lower():
        raise ValueError("baxh5 process failed to complete! (Error)\n stdout: {} \n stderr: {}".
                         format(baxh5_out, baxh5_error))

    if baxh5_process.returncode != 0:
        raise ValueError("baxh5 process failed to complete! (Non-zero return code)\n stdout: {} \n"
                         " stderr: {}".format(baxh5_out, baxh5_error))
    else:
        print("\nbaxh5 process logging:\n stdout:{} \n".format(baxh5_out))
        cmph5_filepath = '/'.join([TEMP_OUTPUT_FOLDER, cmph5file])
        return split_reads_in_cmph5(cmph5_filepath)


def writebaxh5(filecontent, filepath):
    """

    :param filecontent:
    :param filepath:
    :return:
    """
    f = h5py.File(filepath, "w")
    h5data, h5attr = filecontent
    for datasetinfo in h5data:
        datasetattr, datasetdata = datasetinfo
        datasetname = datasetattr[0]
        f.create_dataset(datasetname, data=datasetdata)
    for itemattrinfo in h5attr:
        item_name_type, itemattrs = itemattrinfo
        if item_name_type[1] == H5GROUP:
            itemtmp = f.require_group(item_name_type[0])
        else:
            if item_name_type[0] in f:
                itemtmp = f[item_name_type[0]]
            else:
                # FIXME: this is not recommended, all data of datasets
                # FIXME: should be written in h5data (the last for loop)
                itemtmp = f.create_dataset(item_name_type[0], (1,))
        for itemattr_name, itemattr_val in itemattrs:
            itemtmp.attrs[itemattr_name] = itemattr_val
    f.close()


def split_reads_in_cmph5(cmph5_filepath):
    """

    :param cmph5_filepath:
    :return: info about each read, [(key, value), ...]
    key: ((reffullname, refmd5), target_start, target_end)
    value: dict(), {'/../AlnIndex': a line of AlnIndex, 'AlnArray': numpy.array,
    'DeletionQV': numpy.array, ...}
    """
    cmph5 = h5py.File(cmph5_filepath, 'r')
    ch5Format = CmpH5Format()

    aI = cmph5[ch5Format.ALN_INDEX]
    aIlen = aI.shape[0]
    refId2refFullName = refGroupId2RefInfoIdentifier(cmph5, ch5Format)
    # FIXME: it is better to cal md5 in getRefInfoFromFastaFiles(), but don't know how
    refId2refMd5 = refGroupId2RefInfoIdentifier(cmph5, ch5Format, "MD5")
    movieId2MovieInfo = MovieId2MovieInfo(cmph5, ch5Format)

    reads_rdd = [None] * aIlen
    for i in range(0, aIlen):
        alnIndexTmp = aI[i, :]

        read_key = ((refId2refFullName[alnIndexTmp[ch5Format.REF_ID]],
                     refId2refMd5[alnIndexTmp[ch5Format.REF_ID]]),
                    alnIndexTmp[ch5Format.TARGET_START],
                    alnIndexTmp[ch5Format.TARGET_END])

        read_val = dict()
        # AlnIndex dataset
        read_val[ch5Format.ALN_INDEX] = alnIndexTmp
        alignInfoPath = cmph5[ch5Format.ALN_GROUP_PATH][numpy.where(
            cmph5[ch5Format.ALN_GROUP_ID].value == alnIndexTmp[ch5Format.ALN_ID])[0][0]]
        # align info datasets
        offset_begin = alnIndexTmp[ch5Format.OFFSET_BEGIN]
        offset_end = alnIndexTmp[ch5Format.OFFSET_END]
        for key in cmph5[alignInfoPath].keys():
            read_val[str(key)] = cmph5[alignInfoPath + '/' + str(key)][offset_begin:offset_end]
        # movie info
        movieid = alnIndexTmp[ch5Format.MOVIE_ID]
        read_val.update(movieId2MovieInfo[movieid])

        reads_rdd[i] = (read_key, read_val)
    cmph5.close()
    return reads_rdd


def refGroupId2RefInfoIdentifier(cmph5, ch5Format, refidentifier='FullName'):
    """

    :param cmph5:
    :param ch5Format:
    :param refidentifier: 'FullName' or 'MD5'
    :return:
    """
    refID = cmph5[ch5Format.REF_GROUP_ID]
    refInfoID = cmph5[ch5Format.REF_GROUP_INFO_ID]
    _refInfoID = cmph5['/'.join([ch5Format.REF_INFO, 'ID'])]
    refid = cmph5['/'.join([ch5Format.REF_INFO, refidentifier])]

    refid2refid = {}
    refshape = refID.shape
    for i in range(0, refshape[0]):
        refid2refid[refID[i]] = refid[numpy.where(_refInfoID.value == refInfoID[i])[0][0]]
    return refid2refid


def MovieId2MovieInfo(cmph5, ch5Format, infonames="Name,FrameRate"):
    """

    :param cmph5:
    :param ch5Format:
    :param infonames:names of datasets in MovieInfo, delimited by ","
    "ID,Name,Exp,FrameRate,Run,BindingKit,SequencingKit,SoftwareVersion"
    :return:
    """
    movieinfo = cmph5[ch5Format.MOVIE_INFO]
    movieid = movieinfo["ID"]

    res = {}
    for i in range(0, movieid.shape[0]):
        res[movieid[i]] = {}
        for datasetname in infonames.split(","):
            res[movieid[i]]["/".join([ch5Format.MOVIE_INFO, datasetname.strip()])] = \
                movieinfo[datasetname.strip()][0]
    return res


# group reads of cmp.h5 file------------------------------------------
def _queueChunksForReference(numHits, refLength):
    """
    from ipdsummary.py
    Compute the chunk extents and queue up the work for a single reference
    :param numHits: num of reads aligned to the ref
    :param refLength: ref length
    :return:
    """

    # Number of hits on current reference
    # refGroupId = refInfo.ID
    # numHits = (self.cmph5.RefGroupID == refGroupId).sum()

    # Don't process reference groups with 0 hits.  They may not exist?
    if numHits == 0:
        return

    # Maximum chunk size (set no larger than 1Mb for now)
    MAX_BLOCK_SIZE = 25000

    # Maximum number of hits per chunk
    MAX_HITS = 5000
    # nBases = min(refInfo.Length, refMaxLength)
    nBases = min(refLength, refMaxLength)

    # Adjust numHits if we are only doing part of the contig
    numHits = (numHits * nBases) / refLength

    nBlocks = max([numHits / MAX_HITS, nBases / (MAX_BLOCK_SIZE - 1) + 1])

    # Including nBases / (MAX_BLOCK_SIZE - 1) + 1 in nBlocks calculation:
    # E. coli genome: this should be ~ 10.
    # Human genome: ought to be largest & is meant to ensure that blockSize < MAX_BLOCK_SIZE.

    # Block layout
    blockSize = min(nBases, max(nBases / nBlocks + 1, 1000))
    blockStarts = numpy.arange(0, nBases, step=blockSize)
    blockEnds = blockStarts + blockSize
    blocks = zip(blockStarts, blockEnds)

    return blocks


# operations for each rdd element in cmph5ReadsRDD-------------------------
def basemods_pipeline_cmph5_operations(keyval, moviechemistry, refinfo):
    """

    :param keyval:
    key: ((reffullname, refmd5), (ref_splitting_range_start, ref_splitting_range_end, ref_splitting_range_folds))
    val: [read_value1, read_value2, ...] (read_value is the same as value in read_keyval)
    :param moviechemistry:
    :param refinfo:
    :return: (reffullname, (ref_start, {'csv': csvContent, 'gff': gffContent, }))
    """
    ((reffullname, refmd5), (ref_start, ref_end, ref_folds)) = keyval[0]
    reads_info = list(keyval[1])

    if len(reads_info) > 0:
        if not os.path.isdir(TEMP_OUTPUT_FOLDER):
            os.mkdir(TEMP_OUTPUT_FOLDER, 0777)
        else:
            os.chmod(TEMP_OUTPUT_FOLDER, 0o777)

        # setting paths and variables
        name_prefix = reffullname + '.' + str(ref_start) + '-' + str(ref_end)
        reference_path = SparkFiles.get(ref_filename)
        cmph5_filename = name_prefix + ".cmp.h5"
        cmph5path = '/'.join([TEMP_OUTPUT_FOLDER, cmph5_filename])

        tmpcOnrEX_gff = name_prefix + ".tmpcOnrEX.gff"
        tmpcc5Wn6_csv = name_prefix + ".tmpcc5Wn6.csv"

        refchunkinfo = ','.join(["1", str(ref_start), str(ref_end), str(ref_folds)])

        # writing cmph5 file
        writecmph5(cmph5path, reads_info, reffullname, refmd5, refinfo, moviechemistry)

        # cmph5 operations (sort, repack, computeModifications)
        cmph5_operations = "{scripts_fold}/cmph5_operations.sh {temp_output_folder} {cmph5_filepath}" \
                           " {ref_chunk_info} {reference_filepath} {gff_filename} {csv_filename}"\
            .format(scripts_fold=SCRIPTS_FOLDER,
                    temp_output_folder=TEMP_OUTPUT_FOLDER,
                    cmph5_filepath=cmph5path,
                    ref_chunk_info=refchunkinfo,
                    reference_filepath=reference_path,
                    gff_filename=tmpcOnrEX_gff,
                    csv_filename=tmpcc5Wn6_csv)
        cmph5_process = Popen(shlex.split(cmph5_operations), stdout=PIPE, stderr=PIPE)
        cmph5_out, cmph5_error = cmph5_process.communicate()

        if "[Errno" in cmph5_error.strip() or "error" in cmph5_error.strip().lower():
            raise ValueError("cmph5 process failed to complete! (Error)\n stdout: {} \n stderr: {}".
                             format(cmph5_out, cmph5_error))

        if cmph5_process.returncode != 0:
            raise ValueError("cmph5 process failed to complete! (Non-zero return code)\n stdout: {} \n"
                             " stderr: {}".format(cmph5_out, cmph5_error))
        else:
            print("\ncmph5 process logging:\n stdout:{} \n".format(cmph5_out))
            gff_filepath = '/'.join([TEMP_OUTPUT_FOLDER, tmpcOnrEX_gff])
            csv_filepath = '/'.join([TEMP_OUTPUT_FOLDER, tmpcc5Wn6_csv])

            gffContent, csvContent = "", ""
            with open(gff_filepath) as rf:
                for line in rf:
                    if not line.startswith("##"):
                        gffContent += line.strip() + "\n"
            with open(csv_filepath) as rf:
                next(rf)
                for line in rf:
                    csvContent += line.strip() + "\n"

            return reffullname, (ref_start, {'csv': csvContent, 'gff': gffContent, })
    else:
        raise ValueError("no reads to form a cmph5 file")


def writecmph5(filepath, reads_info, reffullname, refmd5, refinfo, moviechemistry):
    DATASET_ALNINDEX = "AlnInfo/AlnIndex"
    DATASET_MOVIEINFO = {"MovieInfo/Name", "MovieInfo/FrameRate", }
    # DATASET_MOVIECHEMISTRY = {"BindingKit", "SequencingKit", "SoftwareVersion", }
    DATASET_ALIGNSEQINFO = {"QualityValue", "IPD", "DeletionTag", "PulseWidth", "MergeQV",
                            "SubstitutionQV", "InsertionQV", "DeletionQV", "AlnArray", }

    ALIGNREF_PATH = "ref000001"
    ALIGNSEQ_GROUP = '/'.join([ALIGNREF_PATH, "rg1-0"])
    CMPH5_ALN_GROUP_ID = 1
    CMPH5_ALN_GROUP_PATH = '/' + ALIGNSEQ_GROUP
    REFGROUPPATH = '/' + ALIGNREF_PATH
    REFGROUPID = REFINFOID = 1

    f = h5py.File(filepath, "w")
    f_format = CmpH5Format()

    # add AlnGroup
    f.create_dataset(f_format.ALN_GROUP_ID, data=[CMPH5_ALN_GROUP_ID, ],
                     dtype="int32", maxshape=(None,), chunks=(256,))
    f.create_dataset(f_format.ALN_GROUP_PATH, dtype=h5py.special_dtype(vlen=str),
                     data=[CMPH5_ALN_GROUP_PATH, ], maxshape=(None,),
                     chunks=(256,))

    reads_len = len(reads_info)
    alignseq_len = [0] * reads_len
    alnindex_type = reads_info[0][DATASET_ALNINDEX].dtype
    alnindex_ncol = reads_info[0][DATASET_ALNINDEX].shape[0]

    # set alnindex dataset
    alnindex = numpy.ndarray((reads_len, alnindex_ncol), alnindex_type)
    for i in range(0, reads_len):
        alnindextmp = reads_info[i][DATASET_ALNINDEX]
        alnindex[i, :] = alnindextmp
        alignseq_len[i] = alnindextmp[f_format.OFFSET_END] - alnindextmp[f_format.OFFSET_BEGIN]

    alignseq_range = numpy.hstack(([0], numpy.cumsum(alignseq_len)))
    alignseq_range = numpy.array(alignseq_range, dtype=numpy.int32)
    # adjust values in alnindex dataset
    for i in range(0, reads_len):
        alnindex[i, f_format.ID] = i+1
        alnindex[i, f_format.ALN_ID] = CMPH5_ALN_GROUP_ID
        alnindex[i, f_format.MOLECULE_ID] = i
        alnindex[i, f_format.OFFSET_BEGIN] = alignseq_range[i]
        alnindex[i, f_format.OFFSET_END] = alignseq_range[i+1]

    datasetnames = list(reads_info[0].keys())
    # add align seq info
    # datasetnames.remove(DATASET_ALNINDEX)
    for datasetname in datasetnames:
        if datasetname in DATASET_ALIGNSEQINFO:
            datasettype = reads_info[0][datasetname].dtype
            datasetContent = numpy.ndarray((alignseq_range[-1],),
                                           datasettype)
            for i in range(0, reads_len):
                datasetContent[alignseq_range[i]:alignseq_range[i + 1]] = reads_info[i][datasetname]
            f.create_dataset('/'.join([ALIGNSEQ_GROUP, datasetname]), data=datasetContent)
            # else:
            #     raise ValueError("the {} is not proper".format(datasetname))

    # add movieinfo
    movieinfo_items = list(DATASET_MOVIEINFO.intersection(datasetnames))
    movieinfo_content = set()
    for i in range(0, reads_len):
        mcontemp = ()
        for j in range(0, len(movieinfo_items)):
            mcontemp += (reads_info[i][movieinfo_items[j]], )
        movieinfo_content.add(mcontemp)
    movienames = list()
    movieinfo_content = list(movieinfo_content)
    for i in range(0, len(movieinfo_items)):
        movie_item_values = []
        for mcont in movieinfo_content:
            movie_item_values.append(mcont[i])
        if movieinfo_items[i] == "MovieInfo/Name":
            f.create_dataset(movieinfo_items[i], dtype=h5py.special_dtype(vlen=str),
                             data=movie_item_values)
            movienames = movie_item_values
        elif movieinfo_items[i] == "MovieInfo/FrameRate":
            f.create_dataset(movieinfo_items[i],
                             data=numpy.array(movie_item_values, dtype=numpy.float32))
        else:
            raise ValueError("{} name is wrong for movieinfo".format(movieinfo_items[i]))
    if len(movienames) > 0:
        chem_datasetnames = list(moviechemistry[list(moviechemistry.keys())[0]].keys())
        for chem_name in chem_datasetnames:
            chem_values = []
            for moviename in movienames:
                chem_values.append(moviechemistry[moviename][chem_name])
            f.create_dataset('/'.join(["MovieInfo", chem_name]), dtype=h5py.special_dtype(vlen=str),
                             data=chem_values)
        # MovieInfo/ID
        f.create_dataset("MovieInfo/ID", dtype="int32",
                         data=range(1, len(movienames)+1))
        # adjust movie_id in alnindex
        moviename2movieid = {}
        for i in range(0, len(movienames)):
            moviename2movieid[movienames[i]] = i + 1
        for i in range(0, reads_len):
            alnindex[i, f_format.MOVIE_ID] = moviename2movieid[reads_info[i]["MovieInfo/Name"]]

    # add refinfo
    f.create_dataset("RefInfo/FullName", dtype=h5py.special_dtype(vlen=str),
                     data=[reffullname, ])
    f.create_dataset("RefInfo/ID", dtype="int32",
                     data=[REFINFOID, ])
    f.create_dataset("RefInfo/Length", dtype="int32",
                     data=[refinfo[reffullname].getSeqLength(), ])
    # f.create_dataset("RefInfo/MD5", dtype=h5py.special_dtype(vlen=str),
    #                  data=[refinfo[reffullname].getMd5(), ])
    f.create_dataset("RefInfo/MD5", dtype=h5py.special_dtype(vlen=str),
                     data=[refmd5, ])
    # add refgroup
    f.create_dataset("RefGroup/ID", dtype="int32",
                     data=[REFGROUPID, ])
    f.create_dataset("RefGroup/Path", dtype=h5py.special_dtype(vlen=str),
                     data=[REFGROUPPATH, ])
    f.create_dataset("RefGroup/RefInfoID", dtype="int32",
                     data=[REFINFOID, ])
    # adjust ref_id in alnindex
    alnindex[:, f_format.REF_ID] = REFGROUPID

    f.create_dataset(DATASET_ALNINDEX, data=alnindex)
    f.create_group('FileLog')
    # FIXME: needs to figure out the source of attr "Version", and
    # FIXME: how to transport from previous data to here
    f.attrs['Version'] = '2.0.0'

    f.flush()
    f.close()


def creat_redundant_reads(read_keyval, ref_splitting_info):
    """
    check if a read needs to be copyed for two ref chunks
    :param read_keyval:
    :param ref_splitting_info:
    :return: [(key, val), ]
    key: ((reffullname, refmd5), (ref_splitting_range_start, ref_splitting_range_end, ref_splitting_range_folds))
    val: the same as value in read_keyval
    """
    # FIXED: classify the reads by covering ref_splitting_range (start, end)
    # FIXED: or (start-pad, end+pad)?. For now, use (start, end)
    # FIXED: answer: use (start, end) is enough
    (rkey, rval) = read_keyval
    (read_refid, target_start, target_end) = rkey
    read_ref_splitting_info = ref_splitting_info[read_refid]

    res_keyval = []
    for ref_range in read_ref_splitting_info:
        refstart, refend = ref_range[0], ref_range[1]
        # refstart, refend = refstart - PAD, refend + PAD
        if not (target_end <= refstart or target_start >= refend):
            res_keyval.append(((read_refid, ref_range), rval))
    return res_keyval


def adjust_ref_splitting_info(ref_splitting_info, factor=2):
    for refid in ref_splitting_info.keys():
        ref_splitting_range = ref_splitting_info[refid]

        chunk_len = len(ref_splitting_range)
        folds = int(numpy.ceil(chunk_len / float(factor)))
        eachchunk = chunk_len / folds
        chunk_yu = chunk_len % folds

        chunk_steps = [0] * folds
        for i in range(0, folds):
            chunk_steps[i] = eachchunk
            if i < chunk_yu:
                chunk_steps[i] += 1
        chunk_steps_range = numpy.hstack(([0], numpy.cumsum(chunk_steps)))

        ref_splitting_adjustedrange = []
        for i in range(0, folds):
            sbegin, send = chunk_steps_range[i], chunk_steps_range[i+1]

            ref_splitting_adjustedrange.append((ref_splitting_range[sbegin][0],
                                                ref_splitting_range[send-1][1],
                                                send - sbegin))

        ref_splitting_info[refid] = ref_splitting_adjustedrange
    return ref_splitting_info


# operations for each rdd element in modificationRDD-----------------------
def basemods_pipeline_modification_operations(keyval, refinfo):
    reffullname = keyval[0]
    modsinfo = list(keyval[1])

    if not os.path.isdir(TEMP_OUTPUT_FOLDER):
        os.mkdir(TEMP_OUTPUT_FOLDER, 0777)
    else:
        os.chmod(TEMP_OUTPUT_FOLDER, 0o777)

    # setting paths and variables
    name_prefix = reffullname + ".modification"
    reference_path = SparkFiles.get(ref_filename)
    gfffilename = name_prefix + ".gff"
    csvfilename = name_prefix + ".csv"
    motifs_gff_gz_filename = reffullname + ".motifs.gff.gz"

    gfffilepath = '/'.join([TEMP_OUTPUT_FOLDER, gfffilename])
    csvfilepath = '/'.join([TEMP_OUTPUT_FOLDER, csvfilename])
    writemodificationinfo(modsinfo, reffullname, refinfo, gfffilepath, csvfilepath)

    # modification operations (findMotifs, makeMotifGff)
    mods_operations = "{scripts_fold}/mods_operations.sh {temp_output_folder} {gff_filepath}" \
                      " {csv_filepath} {reference_filepath} {motifsgffgz_filename}" \
        .format(scripts_fold=SCRIPTS_FOLDER,
                temp_output_folder=TEMP_OUTPUT_FOLDER,
                gff_filepath=gfffilepath,
                csv_filepath=csvfilepath,
                reference_filepath=reference_path,
                motifsgffgz_filename=motifs_gff_gz_filename)
    mods_process = Popen(shlex.split(mods_operations), stdout=PIPE, stderr=PIPE)
    mods_out, mods_error = mods_process.communicate()

    if "[Errno" in mods_error.strip() or "error" in mods_error.strip().lower():
        raise ValueError("mods process failed to complete! (Error)\n stdout: {} \n stderr: {}".
                         format(mods_out, mods_error))

    if mods_process.returncode != 0:
        raise ValueError("mods process failed to complete! (Non-zero return code)\n stdout: {} \n"
                         " stderr: {}".format(mods_out, mods_error))
    else:
        print("\nmods process logging:\n stdout:{} \n".format(mods_out))
        return motifs_gff_gz_filename


def writemodificationinfo(modsinfo, reffullname, refinfo, gfffilepath, csvfilepath):
    """

    :param modsinfo: [(ref_start, {'csv': csvinfo, 'gff': gffinfo, }), ]
    :param reffullname:
    :param refinfo:
    :param gfffilepath:
    :param csvfilepath:
    :return:
    """
    GFF_HEADER = "##gff-version 3\n" \
                 "##source ipdSummary.py v2.0\n" \
                 "##source-commandline /.../ipdSummary.py \n" \
                 "##sequence-region ref000001 1 {}\n".format(refinfo[reffullname].getSeqLength())
    CSV_HEADER = "refName,tpl,strand,base,score,tMean,tErr,modelPrediction," \
                 "ipdRatio,coverage,frac,fracLow,fracUp\n"

    modsinfo.sort(key=lambda x: x[0])

    with open(gfffilepath, mode='w') as wf:
        wf.write(GFF_HEADER)
        for mods_slice in modsinfo:
            wf.write(mods_slice[1]['gff'])
    with open(csvfilepath, mode='w') as wf:
        wf.write(CSV_HEADER)
        for mods_slice in modsinfo:
            wf.write(mods_slice[1]['csv'])
# -----------------------------------------------------------------------------------------


def basemods_pipe():
    conf = SparkConf().setAppName("Spark-based Pacbio BaseMod pipeline")
    sc = SparkContext(conf=conf)

    # files need to be shared in each node
    sc.addFile('/'.join([reference_directory, ref_filename]))
    sc.addFile('/'.join([reference_directory, ref_sa_filename]))
    # sc.addFile('/'.join([metadataxml_dir, metadataxml_name]))

    # baxh5 file operations
    baxh5_dir = '/home/hadoop/workspace_py/basemods_spark/data/lambda_v210/Analysis_Results'
    baxh5_folds = 1
    baxh5_numpartitions = 1
    # FIXME: will this (append, union) work when the file is large/the memory is not enough?
    baxh5rdds = []
    for filename in os.listdir(baxh5_dir):
        if filename.endswith(".bax.h5"):
            baxh5rdds.append(baxh5toRDD(sc, '/'.join([baxh5_dir, filename]), baxh5_folds, baxh5_numpartitions))
    all_baxh5rdds = sc.union(baxh5rdds)

    # cmph5 file operations
    # FIXME: 1. for the rdd contains every reads in cmph5, which is better: sort first and then group, or
    # FIXME:    group first and then sort?
    # FIXME: 2. keeping "aligned_reads_rdd" in memory saves time, but is a huge waste of memory.
    # FIXME:    need to figure out how not to use "aligned_reads_rdd" twice, so that can save memory.
    aligned_reads_rdd = all_baxh5rdds.flatMap(basemods_pipeline_baxh5_operations).persist()
    # FIXME: x[0][0] is ref's fullname, x[0] is (ref_fullname, ref_md5), check split_reads_in_cmph5()
    ref_indentifiers_count = aligned_reads_rdd.map(lambda (x, y): x[0]).countByValue()

    # reference info to be shared to each node
    # FIXME: 1. (IMPORTANT) for now, don't know how to cal md5 of a sequence,
    # FIXME:    so the md5 info can only be carried with each read.
    # FIXME: 2. keep sequence of the ref in refinfo or not? : not
    refinfos = sc.broadcast(
        getRefInfoFromFastaFiles(['/'.join([reference_directory, ref_filename]), ]))

    # get the ref chunks
    ref_splitting_info = {}
    for ref_id in ref_indentifiers_count.keys():
        # FIXME: ref_id is (ref_fullname, ref_md5), check split_reads_in_cmph5()
        # FIXME: ref_id[0] is ref's fullname
        ref_splitting_info[ref_id] = _queueChunksForReference(ref_indentifiers_count[ref_id],
                                                              refinfos.value[ref_id[0]].getSeqLength())
    # adjust ref_splitting_info
    dfactor = 2
    ref_splitting_info = sc.broadcast(adjust_ref_splitting_info(ref_splitting_info, dfactor))

    adjusted_reads_rdd = aligned_reads_rdd\
        .flatMap(lambda x: creat_redundant_reads(x, ref_splitting_info.value))
    cmph5rdd = adjusted_reads_rdd.groupByKey()

    # movie info to be shared to each node
    moviestriple = sc.broadcast(
        getMoviesChemistry(['/'.join([metadataxml_dir, metadataxml_name]), ]))
    # FIXME: how to cal MovieInfo.FrameRate?
    # FIXME: it was calculated in the "loadPulses (blasr/utils)" step, but the source code
    # FIXME: can't be found, so the FrameRate info can only be carried with every read
    # todo: cal FrameRate
    modification_rdd = cmph5rdd\
        .map(lambda (x, y): basemods_pipeline_cmph5_operations((x, y),
                                                               moviestriple.value,
                                                               refinfos.value))
    motif_rdd = modification_rdd.groupByKey()\
        .map(lambda (x, y): basemods_pipeline_modification_operations((x, y),
                                                                      refinfos.value))
    motif_rdd.count()

    aligned_reads_rdd.unpersist()


if __name__ == '__main__':
    print("spark start------------------------------------")
    basemods_pipe()
    print("spark end--------------------------------------\n"
          "-----------------------------------------------")
