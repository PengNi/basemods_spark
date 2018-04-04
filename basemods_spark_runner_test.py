#!/usr/bin/python
# coding=utf-8
from pyspark import SparkContext, SparkConf, SparkFiles
from pyspark.storagelevel import StorageLevel
from subprocess import Popen, PIPE
from itertools import groupby
from pbcore.io import BasH5Reader, CmpH5Reader
import h5py
import shlex
import os
import fnmatch
import numpy
import hashlib
import ConfigParser
import paramiko
import socket
import time
import random
import shutil
import xml.etree.ElementTree as ET

shell_script_baxh5 = 'baxh5_operations.sh'
shell_script_cmph5 = 'cmph5_operations.sh'
shell_script_mods = 'mods_operations.sh'
shell_script_sa = 'exec_sawriter.sh'
parameters_config = 'parameters_test.conf'

H5GROUP = h5py._hl.group.Group
H5DATASET = h5py._hl.dataset.Dataset

# for change file name
SPACE_ALTER = "_"

# contig attri names
SEQUENCE = 'sequence'
SEQUENCE_LENGTH = 'seqLength'
SEQUENCE_MD5 = 'md5'

refMaxLength = 3e12
COLUMNS = 60
PAD = 15

max_numpartitions = 10000
reads_shuffle_factor = 3

# sleep seconds when get data from master node/HDFS
MAX_SLEEP_SECONDS = 100

# for split reference to multi chunks
max_chunk_length = 25000
max_reads_per_chunk = 5000
limitation_readsnum = max_reads_per_chunk * 5

# dir for save the results in HDFS
HDFS_IPDINFO_DIR = '/fastipd'
HDFS_SAMIPDINFO_DIR = '/samipd'
HDFS_MODS_DIR = '/mods'


# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# read configfile to get parameters------------------------------------------
def getParametersFromFile(config_file):
    conf = ConfigParser.ConfigParser()
    conf.read(config_file)

    global TEMP_OUTPUT_FOLDER
    global SMRT_ANALYSIS_HOME
    global DATA_SAVE_MODE
    global REFERENCE_DIR
    global REF_FILENAME
    global REF_SA_FILENAME
    global CELL_DATA_DIR
    global CMPH5_DATA_DIR

    global HDFS_CMD
    global PROC_NUM
    global BAXH5_FOLDS
    # global REF_CHUNKS_FACTOR
    global READS_TRIM_STRATEGY
    global IPDMAXCOVERAGE
    global METHYLATION_TYPES

    global GET_IPD_FROM_BASH5
    global GET_IPD_FROM_CMPH5

    global MASTERNODE_IP
    global MASTERNODE_PORT
    global MASTERNODE_USERNAME
    global MASTERNODE_USERPASSWD

    global SPARK_EXECUTOR_MEMORY
    global SPARK_TASK_CPUS
    global SPARK_MEMORY_FRACTION
    global SPARK_MEMORY_STORAGEFRACTION

    TEMP_OUTPUT_FOLDER = conf.get("FilePath", "TEMP_OUTPUT_FOLDER")
    SMRT_ANALYSIS_HOME = conf.get("FilePath", "SMRT_ANALYSIS_HOME")
    DATA_SAVE_MODE = conf.get("FilePath", "DATA_SAVE_MODE")
    REFERENCE_DIR = conf.get("FilePath", "REFERENCE_DIR")
    REF_FILENAME = conf.get("FilePath", "REF_FILENAME")
    REF_SA_FILENAME = conf.get("FilePath", "REF_SA_FILENAME")
    CELL_DATA_DIR = conf.get("FilePath", "CELL_DATA_DIR")
    CMPH5_DATA_DIR = conf.get("FilePath", "CMPH5_DATA_DIR")

    HDFS_CMD = conf.get("PipelineArgs", "HDFS_CMD")
    PROC_NUM = conf.getint("PipelineArgs", "PROC_NUM")
    BAXH5_FOLDS = conf.getint("PipelineArgs", "BAXH5_FOLDS")
    # REF_CHUNKS_FACTOR = conf.getint("PipelineArgs", "REF_CHUNKS_FACTOR")
    READS_TRIM_STRATEGY = conf.get("PipelineArgs", "READS_TRIM_STRATEGY")
    IPDMAXCOVERAGE = conf.get("PipelineArgs", 'IPDMAXCOVERAGE')
    # METHYLATION_TYPES = conf.get("PipelineArgs", "METHYLATION_TYPES")
    METHYLATION_TYPES = trim_spaces(conf.get("PipelineArgs", "METHYLATION_TYPES"))

    GET_IPD_FROM_BASH5 = conf.get("PipelineArgs", "GET_IPD_FROM_BASH5")
    GET_IPD_FROM_CMPH5 = conf.get("PipelineArgs", "GET_IPD_FROM_CMPH5")

    MASTERNODE_IP = conf.get("MasterNodeInfo", "HOST")
    MASTERNODE_PORT = conf.getint("MasterNodeInfo", "HOSTPORT")
    MASTERNODE_USERNAME = conf.get("MasterNodeInfo", "USERNAME")
    MASTERNODE_USERPASSWD = conf.get("MasterNodeInfo", "USERPASSWD")

    SPARK_EXECUTOR_MEMORY = conf.get('SparkConfiguration', 'spark_executor_memory')
    SPARK_TASK_CPUS = conf.get('SparkConfiguration', 'spark_task_cpus')
    SPARK_MEMORY_FRACTION = conf.get('SparkConfiguration', 'spark_memory_fraction')
    SPARK_MEMORY_STORAGEFRACTION = conf.get('SparkConfiguration', 'spark_memory_storageFraction')
    return


# for trim spaces
def trim_spaces(ori_str):
    return str(ori_str).strip().replace(", ", ",")


# fasta_info.py-------------------------
def getRefInfoFromFastaFiles(filepaths):
    """

    :param filepaths: list of filepath
    :return:
    """
    contiginfo_dict = {}
    for fp in filepaths:
        contiginfo_dict.update(_getRefInfoFromFastaFile(fp))
    return contiginfo_dict


def _getRefInfoFromFastaFile(filepath):
    contiginfos = FastaInfo(filepath).getContigs()

    contiginfo_dict = {}
    for contiginfo in contiginfos:
        # contiginfo_dict[contiginfo.getContigName()] = contiginfo
        contiginfo_dict[contiginfo.getContigName()] = {}
        contiginfo_dict[contiginfo.getContigName()][SEQUENCE] = contiginfo.getSequence()
        contiginfo_dict[contiginfo.getContigName()][SEQUENCE_LENGTH] = contiginfo.getSeqLength()
        contiginfo_dict[contiginfo.getContigName()][SEQUENCE_MD5] = contiginfo.getMd5()
    del contiginfos
    return contiginfo_dict


class FastaInfo:
    def __init__(self, filepath):
        self._contigs = []  # list of ContigInfo()

        with open(filepath, mode='r') as rf:
            contigTmp = ContigInfo()
            firstline = next(rf)
            if not str(firstline).startswith(">"):
                print("read fasta file wrong!")
                raise ValueError
            else:
                contigTmp.setContigName(firstline.strip()[1:])

            sequencetmp = ""
            for line in rf:
                if line.startswith(">"):
                    contigTmp.setSequence(sequencetmp)
                    contigTmp.setSeqLength(len(sequencetmp))
                    contigTmp.setMd5(hashlib.md5(sequencetmp).hexdigest())
                    self._contigs.append(contigTmp)

                    sequencetmp = ""
                    contigTmp = ContigInfo()
                    contigTmp.setContigName(line.strip()[1:])
                else:
                    sequencetmp += line.strip()

            # the last contig
            contigTmp.setSequence(sequencetmp)
            contigTmp.setSeqLength(len(sequencetmp))
            contigTmp.setMd5(hashlib.md5(sequencetmp).hexdigest())
            self._contigs.append(contigTmp)

    def getContigs(self):
        return self._contigs


class ContigInfo:
    def __init__(self):
        self._contigName = ""
        self._sequence = ""
        self._seqLength = 0
        self._md5 = ""

    def getContigName(self):
        return self._contigName

    def setContigName(self, contigname):
        self._contigName = contigname

    def getSequence(self):
        return self._sequence

    def setSequence(self, sequence):
        self._sequence = sequence

    def getSeqLength(self):
        return self._seqLength

    def setSeqLength(self, seqlength):
        self._seqLength = seqlength

    def getMd5(self):
        return self._md5

    def setMd5(self, md5):
        self._md5 = md5


# CmpH5Format.py-----------------------
class CmpH5Format:
    # def __init__(self, cmpH5):
    def __init__(self):
        # if ('Version' in cmpH5.attrs):
        #     self.VERSION = cmpH5.attrs['Version']

        self.ALN_INFO = 'AlnInfo'
        self.REF_INFO = 'RefInfo'
        self.MOVIE_INFO = 'MovieInfo'
        self.REF_GROUP = 'RefGroup'
        self.ALN_GROUP = 'AlnGroup'
        self.ALN_INDEX_NAME = 'AlnIndex'
        self.FILE_LOG = 'FileLog'
        self.BARCODE_INFO = 'BarcodeInfo'

        self.ALN_INDEX = '/'.join([self.ALN_INFO, self.ALN_INDEX_NAME])
        self.REF_GROUP_ID = '/'.join([self.REF_GROUP, 'ID'])
        self.REF_GROUP_PATH = '/'.join([self.REF_GROUP, 'Path'])
        self.REF_GROUP_INFO_ID = '/'.join([self.REF_GROUP, 'RefInfoID'])

        self.REF_OFFSET_TABLE = '/'.join([self.REF_GROUP, 'OffsetTable'])
        self.ALN_GROUP_ID = '/'.join([self.ALN_GROUP, 'ID'])
        self.ALN_GROUP_PATH = '/'.join([self.ALN_GROUP, 'Path'])

        # Movie Info
        self.MOVIE_INFO_ID = '/'.join([self.MOVIE_INFO, 'ID'])
        self.MOVIE_INFO_NAME = '/'.join([self.MOVIE_INFO, 'Name'])
        self.MOVIE_INFO_EXP = '/'.join([self.MOVIE_INFO, 'Exp'])
        self.MOVIE_INFO_FRAME_RATE = '/'.join([self.MOVIE_INFO, 'FrameRate'])
        self.MOVIE_INFO_RUN = '/'.join([self.MOVIE_INFO, 'Run'])
        self.MOVIE_INFO_BINDING_KIT = '/'.join([self.MOVIE_INFO, 'BindingKit'])
        self.MOVIE_INFO_SEQUENCING_KIT = '/'.join([self.MOVIE_INFO, 'SequencingKit'])
        self.MOVIE_INFO_SOFTWARE_VERSION = '/'.join([self.MOVIE_INFO, 'SoftwareVersion'])

        (self.ID, self.ALN_ID, self.MOVIE_ID, self.REF_ID, self.TARGET_START,
         self.TARGET_END, self.RC_REF_STRAND, self.HOLE_NUMBER, self.SET_NUMBER,
         self.STROBE_NUMBER, self.MOLECULE_ID, self.READ_START, self.READ_END,
         self.MAP_QV, self.N_MATCHES, self.N_MISMATCHES, self.N_INSERTIONS,
         self.N_DELETIONS, self.OFFSET_BEGIN, self.OFFSET_END, self.N_BACK,
         self.N_OVERLAP) = range(0, 22)

        # self.extraTables = ['/'.join([self.ALN_INFO, x]) for x in
        #                     cmpH5[self.ALN_INFO].keys()
        #                     if not x == self.ALN_INDEX_NAME]
        # sorting
        self.INDEX_ATTR = "Index"
        self.INDEX_ELTS = ['REF_ID', 'TARGET_START', 'TARGET_END']


# movie_chemistry.py-----------------
def getMoviesChemistry(filepaths):
    """

    :param filepaths: list of filepath
    :return:
    """
    moviestriple = {}
    for filepath in filepaths:
        moviestriple.update(MovieChemistry(filepath).getMovieTriple())
    return moviestriple


class ChemistryLookupError(Exception):
    pass


class MovieChemistry:
    def __init__(self, filepath):
        self._movietriple = {}
        if str(filepath).endswith(".xml"):
            moviename = str(filepath).split("/")[-1].split(".")[0]
            triple = self.tripleFromMetadataXML(filepath)

            tripledict = dict()
            # tripledict['Name'] = moviename
            tripledict['BindingKit'] = triple[0]
            tripledict['SequencingKit'] = triple[1]
            tripledict['SoftwareVersion'] = triple[2]
            self._movietriple[moviename] = tripledict
        else:
            pass

    def getMovieTriple(self):
        return self._movietriple

    def tripleFromMetadataXML(self, metadataXmlPath):
        """
        from pbcore.io.BasH5IO
        Scrape the triple from the metadata.xml, or exception if the file
        or the relevant contents are not found
        """
        nsd = {None: "http://pacificbiosciences.com/PAP/Metadata.xsd",
               "pb": "http://pacificbiosciences.com/PAP/Metadata.xsd"}
        try:
            tree = ET.parse(metadataXmlPath)
            root = tree.getroot()
            bindingKit = root.find("pb:BindingKit/pb:PartNumber", namespaces=nsd).text
            sequencingKit = root.find("pb:SequencingKit/pb:PartNumber", namespaces=nsd).text
            # The instrument version is truncated to the first 3 dot components
            instrumentControlVersion = root.find("pb:InstCtrlVer", namespaces=nsd).text
            verComponents = instrumentControlVersion.split(".")[0:2]
            instrumentControlVersion = ".".join(verComponents)
            return (bindingKit, sequencingKit, instrumentControlVersion)
        except Exception as e:
            raise ChemistryLookupError, \
                ("Could not find, or extract chemistry information from, %s" % (metadataXmlPath,))


# bash5_splitting.py-------------------
def split_holenumbers(holenumbers, folds=1):
    holenumbers_len = len(holenumbers)
    if folds > holenumbers_len:
        folds = holenumbers_len
    onefoldlen_base = holenumbers_len / folds
    fold_yu = holenumbers_len % folds

    hole_splitspots = [onefoldlen_base] * folds
    hole_splitspots[0] += fold_yu
    endOffset = numpy.cumsum(hole_splitspots)
    beginOffset = numpy.hstack(([0], endOffset[0:-1]))
    offsets = zip(beginOffset, endOffset)
    return offsets


def makeOffsetsDataStructure(baxh5obj):
    """
    from pbcore.io.BasH5IO.py
    :param baxh5obj:
    :return:
    """
    numEvent = baxh5obj["/PulseData/BaseCalls/ZMW/NumEvent"].value
    holeNumber = baxh5obj["/PulseData/BaseCalls/ZMW/HoleNumber"].value
    endOffset = numpy.cumsum(numEvent)
    beginOffset = numpy.hstack(([0], endOffset[0:-1]))
    offsets = zip(beginOffset, endOffset)
    return dict(zip(holeNumber, offsets))


def get_movieName(baxh5file):
    """
    from pbcore.io.BasH5IO.py.movieName
    :param baxh5file:
    :return:
    """
    movieNameAttr = baxh5file["/ScanData/RunInfo"].attrs["MovieName"]

    # In old bas.h5 files, attributes of ScanData/RunInfo are stored as
    # strings in arrays of length one.
    if (isinstance(movieNameAttr, (numpy.ndarray, list)) and
                len(movieNameAttr) == 1):
        movieNameString = movieNameAttr[0]
    else:
        movieNameString = movieNameAttr

    if not isinstance(movieNameString, basestring):
        raise TypeError("Unsupported movieName {m} of type {t}."
                        .format(m=movieNameString,
                                t=type(movieNameString)))
    return movieNameString


# scp.py---------------------------------------------------
def ssh_scp_put(ip, port, user, password, local_file, remote_file):
    """

    :param ip:
    :param port: int
    :param user:
    :param password:
    :param local_file:
    :param remote_file:
    :return:
    """
    flag = 0
    try:
        paramiko.util.log_to_file('paramiko.log')
        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        ssh.connect(ip, port, user, password)
        sftp = paramiko.SFTPClient.from_transport(ssh.get_transport())
        sftp = ssh.open_sftp()
        sftp.put(local_file, remote_file)
        flag = 1
        sftp.close()
        ssh.close()
    except OSError:
        print("wrong put connection {} {} {} {} {}".format(ip, port, user, local_file, remote_file))
    finally:
        return flag


def ssh_scp_get(ip, port, user, password, remote_file, local_file):
    """

    :param ip:
    :param port: int
    :param user:
    :param password:
    :param remote_file:
    :param local_file:
    :return:
    """
    flag = 0
    try:
        paramiko.util.log_to_file('paramiko.log')
        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        ssh.connect(ip, port, user, password)
        sftp = paramiko.SFTPClient.from_transport(ssh.get_transport())
        sftp = ssh.open_sftp()
        sftp.get(remote_file, local_file)
        flag = 1
        sftp.close()
        ssh.close()
    except OSError:
        print("wrong get connection {} {} {} {} {}".format(ip, port, user, remote_file, local_file))
    finally:
        return flag


def worker_put_master(mip, mport, muser, mpassword, wfile, mfile, max_sleep_seconds=1):
    try:
        attemp_times = 100
        for i in range(0, attemp_times):
            expected_sleep_seconds = random.randint(0, max_sleep_seconds) * (i + 1)
            actual_sleep_seconds = expected_sleep_seconds \
                if expected_sleep_seconds < max_sleep_seconds else max_sleep_seconds
            time.sleep(actual_sleep_seconds)
            issuccess = ssh_scp_put(mip, mport, muser, mpassword,
                                    wfile, mfile)
            if issuccess > 0:
                print("{} {} success".format(mip, mfile))
                # todo write a success flag file
                break
    except Exception:
        print('wrong connection remote_filepath: {}'.format(mfile))
    finally:
        print('done transmitting data from local node to {}'.format(mfile))
        return mfile
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------------
# calculate .fasta.sa file --------------------------------------------
def exec_sawriter(sa_script_path, ref_fasta_filepath):
    sa_operations = "{sa_script} {seymour_home} {ref_fasta_file}".\
        format(sa_script=sa_script_path,
               seymour_home=SMRT_ANALYSIS_HOME,
               ref_fasta_file=ref_fasta_filepath)
    sa_process = Popen(sa_operations, stdout=PIPE, stderr=PIPE, shell=True)
    sa_out, sa_error = sa_process.communicate(sa_operations)

    if "[Errno" in sa_error.strip() or "error" in sa_error.strip().lower():
        raise ValueError("sa process failed to complete!\n"
                         "try again or use blasr.sawriter to generate the fasta.sa file of your"
                         "REFERENCE(.fasta) file manually\n"
                         "(Error)\n stdout: {} \n stderr: {}".
                         format(sa_out, sa_error))

    if sa_process.returncode != 0:
        raise ValueError("sa process failed to complete!\n"
                         "try again or use blasr.sawriter to generate the fasta.sa file of your"
                         "REFERENCE(.fasta) file manually\n"
                         "(Non-zero return code)\n stdout: {} \n"
                         " stderr: {}".format(sa_out, sa_error))
    else:
        print("\nsa process logging:\n stdout:{} \n".format(sa_out))
        return os.path.basename(ref_fasta_filepath) + ".sa"


def rename(oriname):
    return oriname.replace(' ', SPACE_ALTER)\
        .replace('/', SPACE_ALTER)\
        .replace(':', SPACE_ALTER)\
        .replace(';', SPACE_ALTER)


# name ref contig file
def name_reference_contig_file(ref_name, contigname):
    ref_prefix, ref_ext = os.path.splitext(ref_name)
    contigname = rename(contigname)
    return ref_prefix + '.' + contigname + ref_ext


# for write ipd info --------------------------------------------------
def asciiFromQvs(a):
    return (numpy.clip(a, 0, 93).astype(numpy.uint8) + 33).tostring()


def zmwReads(inBasH5, readType='subreads'):
    """
    Extract all reads of the appropriate read type
    """
    for zmw in inBasH5:
        if readType == "ccs":
            r = zmw.ccsRead
            if r:
                yield r
        elif readType == "unrolled":
            yield zmw.read()
        else:
            for r in zmw.subreads:
                yield r


def write_ipd_of_bash5(inBasH5File, outIpdInfoFile):
    DELIMITER1 = "@"
    DELIMITER2 = "+quality value"
    DELIMITER3 = "+IPD value"

    start = time.time()
    inBasH5 = BasH5Reader(inBasH5File)
    issuccess = 0
    try:
        outIpdInfo = open(outIpdInfoFile, 'w')
        for zmwRead in zmwReads(inBasH5):
            readinfotmp = "\n".join([DELIMITER1 + zmwRead.readName,
                                     zmwRead.basecalls(),
                                     DELIMITER2,
                                     asciiFromQvs(zmwRead.QualityValue()),
                                     DELIMITER3,
                                     " ".join(map(str, zmwRead.IPD()))])
            outIpdInfo.write(readinfotmp + '\n')
        outIpdInfo.flush()
        outIpdInfo.close()
        print("the IPD info of {} has been extracted and written.\ncost {} seconds"
              .format(inBasH5File, time.time() - start))
        issuccess = 1
    except IOError:
        print("IOError while writing {}".format(outIpdInfoFile))
    finally:
        return issuccess


def get_ipdvalue_of_h5_to_master(inH5File, master_data_dir, h5Type="BAXH5",
                                 max_sleep_seconds=1):
    """

    :param inH5File:
    :param master_data_dir:
    :param h5Type: "BAXH5" or "CMPH5"
    :param max_sleep_seconds:
    :return:
    """
    if h5Type == "BAXH5":
        name_prefix = os.path.basename(inH5File).split(".bax.h5")[0]
        outIpdInfoFileName = name_prefix + ".fastipd"
        dirpath = os.path.dirname(inH5File)
        outIpdInfoFile = '/'.join([dirpath, outIpdInfoFileName])
        wissuccess = write_ipd_of_bash5(inH5File, outIpdInfoFile)
    elif h5Type == "CMPH5":
        name_prefix = os.path.basename(inH5File).split(".cmp.h5")[0]
        outIpdInfoFileName = name_prefix + ".samipd"
        dirpath = os.path.dirname(inH5File)
        outIpdInfoFile = '/'.join([dirpath, outIpdInfoFileName])
        wissuccess = write_ipd_of_cmph5(inH5File, outIpdInfoFile)
    else:
        print("arg h5Type is not set rightly")
        return inH5File
    if wissuccess:
        master_ip = MASTERNODE_IP
        master_port = MASTERNODE_PORT
        master_user = MASTERNODE_USERNAME
        master_passwd = MASTERNODE_USERPASSWD

        remote_filepath = '/'.join([master_data_dir, outIpdInfoFileName])
        try:
            attemp_times = 100
            for i in range(0, attemp_times):
                expected_sleep_seconds = random.randint(0, max_sleep_seconds) * (i + 1)
                actual_sleep_seconds = expected_sleep_seconds \
                    if expected_sleep_seconds < max_sleep_seconds else max_sleep_seconds
                time.sleep(actual_sleep_seconds)
                issuccess = ssh_scp_put(master_ip, master_port, master_user, master_passwd,
                                        outIpdInfoFile, remote_filepath)
                if issuccess > 0:
                    print("{} {} success".format(master_ip, remote_filepath))
                    # todo write a success flag file
                    break
        except Exception:
            print('wrong connection remote_filepath: {}'.format(remote_filepath))
        finally:
            os.remove(outIpdInfoFile)  # rm temp ipdinfo file
            print('done transmitting data from local node to {}'.format(remote_filepath))
            return inH5File
    else:
        print("failed to write {}".format(outIpdInfoFile))
        return inH5File


def get_ipdvalue_of_h5_to_hdfs(inH5File, hdfs_data_dir, h5Type="BAXH5", max_sleep_seconds=1):
    """

    :param inH5File:
    :param hdfs_data_dir:
    :param h5Type: "BAXH5" or "CMPH5"
    :param max_sleep_seconds:
    :return:
    """
    if h5Type == "BAXH5":
        name_prefix = os.path.basename(inH5File).split(".bax.h5")[0]
        dirpath = os.path.dirname(inH5File)
        outIpdInfoFileName = name_prefix + ".fastipd"
        outIpdInfoFile = '/'.join([dirpath, outIpdInfoFileName])
        wissuccess = write_ipd_of_bash5(inH5File, outIpdInfoFile)
    elif h5Type == "CMPH5":
        name_prefix = os.path.basename(inH5File).split(".cmp.h5")[0]
        outIpdInfoFileName = name_prefix + ".samipd"
        dirpath = os.path.dirname(inH5File)
        outIpdInfoFile = '/'.join([dirpath, outIpdInfoFileName])
        wissuccess = write_ipd_of_cmph5(inH5File, outIpdInfoFile)
    else:
        print("arg h5Type is not set rightly")
        return inH5File
    if wissuccess:
        cmd_output, cmd_errors = run_cmd_safe([HDFS_CMD, 'dfs', '-copyFromLocal', '-f',
                                               outIpdInfoFile,
                                               '/'.join([hdfs_data_dir, outIpdInfoFileName])],
                                              max_sleep_seconds)
        os.remove(outIpdInfoFile)  # rm temp ipdinfo file
    else:
        print("failed to write {}".format(outIpdInfoFile))
    return inH5File


def get_ipdvalue_of_h5_to_sharedfolder(inH5File, shared_dir, h5Type="BAXH5"):
    """

    :param inH5File:
    :param shared_dir:
    :param h5Type: "BAXH5" or "CMPH5"
    :return:
    """
    if h5Type == "BAXH5":
        name_prefix = os.path.basename(inH5File).split(".bax.h5")[0]
        outIpdInfoFileName = name_prefix + ".fastipd"

        dirpath = os.path.dirname(inH5File)
        if dirpath != shared_dir:
            dirpath = shared_dir

        outIpdInfoFile = '/'.join([dirpath, outIpdInfoFileName])
        wissuccess = write_ipd_of_bash5(inH5File, outIpdInfoFile)
    elif h5Type == "CMPH5":
        name_prefix = os.path.basename(inH5File).split(".cmp.h5")[0]
        outIpdInfoFileName = name_prefix + ".samipd"

        dirpath = os.path.dirname(inH5File)
        if dirpath != shared_dir:
            dirpath = shared_dir

        outIpdInfoFile = '/'.join([dirpath, outIpdInfoFileName])
        wissuccess = write_ipd_of_cmph5(inH5File, outIpdInfoFile)
    else:
        print("arg h5Type is not set rightly")
        return inH5File


# get baxh5file from master node---------------------------------------
def get_baxh5file_from_masternode(remote_filepath, local_temp_dir, max_sleep_seconds=1):
    master_ip = MASTERNODE_IP
    master_port = MASTERNODE_PORT
    master_user = MASTERNODE_USERNAME
    master_passwd = MASTERNODE_USERPASSWD

    if not os.path.isdir(local_temp_dir):
        try:
            os.mkdir(local_temp_dir, 0777)
        except:
            print('local temp directory {} exists.'.format(local_temp_dir))

    filename = os.path.basename(remote_filepath)
    local_filepath = '/'.join([local_temp_dir, filename])
    try:
        attemp_times = 100
        for i in range(0, attemp_times):
            expected_sleep_seconds = random.randint(0, max_sleep_seconds) * (i + 1)
            actual_sleep_seconds = expected_sleep_seconds \
                if expected_sleep_seconds < max_sleep_seconds else max_sleep_seconds
            time.sleep(actual_sleep_seconds)
            issuccess = ssh_scp_get(master_ip, master_port, master_user, master_passwd,
                                    remote_filepath, local_filepath)
            if issuccess > 0:
                print("{} {} success".format(master_ip, remote_filepath))
                # todo write a success flag file
                break
        # ----transmit ipd value to master node
        # FIXME: need to delete this when it is useless
        if str(GET_IPD_FROM_BASH5).lower() == 'yes':
            get_ipdvalue_of_h5_to_master(local_filepath, CELL_DATA_DIR, "BAXH5",
                                         max_sleep_seconds)
    except Exception:
        print('wrong connection local_filepath: {}'.format(local_filepath))
    finally:
        print('done transmitting data from master node to {}'.format(local_filepath))
        return local_filepath


def get_baxh5file_from_hdfs(hdfs_filepath, local_temp_dir, max_sleep_seconds=1):
    if not os.path.isdir(local_temp_dir):
        try:
            os.mkdir(local_temp_dir, 0777)
        except:
            print('local temp directory {} exists.'.format(local_temp_dir))
    filename = os.path.basename(hdfs_filepath)
    local_filepath = '/'.join([local_temp_dir, filename])
    cmd_output, cmd_errors = run_hdfs_get_cmd(hdfs_filepath, local_filepath,
                                              max_sleep_seconds)
    if "[Errno" in cmd_errors.strip() or "error" in cmd_errors.strip().lower():
        raise RuntimeError("hdfs get process failed to complete! (Error)\n stdout: {} \n stderr: {}".
                           format(cmd_output, cmd_errors))
    else:
        # print("\nhdfs get process logging:\n stdout:{} \n".format(cmd_output))
        # ----transmit ipd value to HDFS
        # FIXME: need to delete this when it is useless
        if str(GET_IPD_FROM_BASH5).lower() == 'yes':
            # hdfs_dir = os.path.dirname(hdfs_filepath)
            hdfs_dir = CELL_DATA_DIR + HDFS_IPDINFO_DIR
            get_ipdvalue_of_h5_to_hdfs(local_filepath, hdfs_dir, "BAXH5")
        return local_filepath


def get_baxh5file_from_sharedfolder(shareddir_filepath):
    try:
        # ----transmit ipd value to master node
        # FIXME: need to delete this when it is useless
        if str(GET_IPD_FROM_BASH5).lower() == 'yes':
            get_ipdvalue_of_h5_to_sharedfolder(shareddir_filepath, CELL_DATA_DIR,
                                               "BAXH5")
    except Exception:
        print('wrong connection shared_folder: {}'.format(shareddir_filepath))
    finally:
        print('done transmitting data from shared folder to {}'.format(shareddir_filepath))
        return shareddir_filepath


# convert baxh5 to list------------------------------------------------
def get_chunks_of_baxh5file(baxh5file, folds=1):
    f = h5py.File(baxh5file, "r")
    holenumbers = f['/PulseData/BaseCalls/ZMW/HoleNumber'].value
    if folds > len(holenumbers):
        folds = len(holenumbers)
    elif folds < 1:
        folds = 1
    hole_splitspots = split_holenumbers(holenumbers, folds)
    hole2range = makeOffsetsDataStructure(f)  # todo: how to make it faster
    basecall_splitspots = get_basecall_range_of_each_holesblock(hole_splitspots,
                                                                holenumbers, hole2range)
    moviename = get_movieName(f)

    chunk_data_info = []

    # datasets in PulseData/BaseCalls/ZMW
    # FIXME: use (for key in f['']) or (for key in f[''].keys())?
    for key in f['/PulseData/BaseCalls/ZMW']:
        chunk_data_info.extend(get_chunks_in_zmw_info(hole_splitspots, holenumbers,
                                                      '/PulseData/BaseCalls/ZMW/' + str(key),
                                                      moviename))

    # datasets in PulseData/BaseCalls/ZMWMetrics
    for key in f['/PulseData/BaseCalls/ZMWMetrics']:
        chunk_data_info.extend(get_chunks_in_zmw_info(hole_splitspots, holenumbers,
                                                      '/PulseData/BaseCalls/ZMWMetrics/' + str(key),
                                                      moviename))

    # datasets in PulseData/BaseCalls
    for key in f['/PulseData/BaseCalls']:
        if not (str(key) == 'ZMWMetrics' or str(key) == 'ZMW'):
            chunk_data_info.extend(get_chunks_in_basecalls_info(hole_splitspots, holenumbers,
                                                                basecall_splitspots,
                                                                '/PulseData/BaseCalls/' + str(key),
                                                                moviename))

    # PulseData/Regions
    chunk_data_info.extend(get_chunks_in_region_info(hole_splitspots, holenumbers, f,
                                                     '/PulseData/Regions', moviename))
    # baxh5attrs = get_baxh5_attrs(f)
    f.close()

    chunk_data_group = group_folds_of_one_baxh5file(chunk_data_info)
    del chunk_data_info

    print('done splitting {} to {} chunk(s)'.format(baxh5file, folds))
    return map(lambda x: add_each_fold_the_filepath(x, baxh5file),
               chunk_data_group)


def add_each_fold_the_filepath(x, baxh5filepath):
    return x[0], (baxh5filepath, x[1])


def group_folds_of_one_baxh5file(holerange_data):
    group_holerange_data = []
    for key, group in groupby(sorted(holerange_data), lambda x: x[0]):
        data_set = []
        for h5data in group:
            data_set.append(h5data[1])
        group_holerange_data.append((key, data_set))
    return group_holerange_data


def get_chunks_in_zmw_info(hole_splitspots,
                           holenumbers, datasetname, moviename):
    chunks_info = []
    for holess in hole_splitspots:
        chunks_info.append(((moviename, (holenumbers[holess[0]], holenumbers[holess[1] - 1])),
                            (datasetname, (holess[0], holess[1]))))
    return chunks_info


def get_chunks_in_basecalls_info(hole_splitspots, holenumbers,
                                 basecall_splitspots, datasetname, moviename):
    chunks_info = []
    for i in range(0, len(hole_splitspots)):
        chunks_info.append(((moviename, (holenumbers[hole_splitspots[i][0]],
                                         holenumbers[hole_splitspots[i][1] - 1])),
                            (datasetname, (basecall_splitspots[i][0],
                                           basecall_splitspots[i][1]))))
    return chunks_info


def get_chunks_in_region_info(hole_splitspots,
                              holenumbers, f, datasetname, moviename):
    regiondata = f[datasetname]
    sshape = regiondata.shape

    chunks_info = []
    locs_start = 0
    for holess in hole_splitspots[:-1]:
        holenumbers_tmp = set(holenumbers[holess[0]:holess[1]])
        for i in range(locs_start, sshape[0]):
            if regiondata[i, 0] not in holenumbers_tmp:
                chunks_info.append(((moviename, (holenumbers[holess[0]], holenumbers[holess[1] - 1])),
                                   (datasetname, (locs_start, i))))
                locs_start = i
                break
    chunks_info.append(((moviename, (holenumbers[hole_splitspots[-1][0]],
                                     holenumbers[hole_splitspots[-1][1] - 1])),
                       (datasetname, (locs_start, sshape[0]))))
    return chunks_info


def get_basecall_range_of_each_holesblock(hole_splitspots, holenumbers, holerange):
    basecall_splitspots = []
    for hole_splitspot in hole_splitspots:
        begin, end = holerange[holenumbers[hole_splitspot[0]]][0], \
                     holerange[holenumbers[hole_splitspot[1] - 1]][1]
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
    ((moviename, holerange), (filepath,
    [(datasetname, (dataset_begin_spot, dataset_end_spot)),...]))
    :param keyval:
    :return: aligned reads
    """
    fileinfo, filecontent = keyval

    if DATA_SAVE_MODE == 'MASTER' or DATA_SAVE_MODE == 'HDFS':
        reference_path = SparkFiles.get(REF_FILENAME)
        referencesa_path = SparkFiles.get(USED_REF_SA_FILENAME)
    elif DATA_SAVE_MODE == 'SHARED_FOLDER':
        reference_path = '/'.join([REFERENCE_DIR, REF_FILENAME])
        referencesa_path = '/'.join([REFERENCE_DIR, USED_REF_SA_FILENAME])
    else:
        reference_path = referencesa_path = ''
    baxh5_shell_file_path = SparkFiles.get(shell_script_baxh5)

    if not os.path.isdir(TEMP_OUTPUT_FOLDER):
        try:
            os.mkdir(TEMP_OUTPUT_FOLDER, 0777)
        except:
            print('temp directory {} exists.'.format(TEMP_OUTPUT_FOLDER))

    if BAXH5_FOLDS == 1:
        baxh5path = filecontent[0]
        name_prefix = os.path.basename(baxh5path).split('.bax.h5')[0]
    else:
        name_prefix = (fileinfo[0] + "." + str(fileinfo[1][0]) + "-" + str(fileinfo[1][1])).replace(' ', SPACE_ALTER)
        baxh5file = name_prefix + ".bax.h5"
        baxh5_dir = TEMP_OUTPUT_FOLDER
        baxh5path = baxh5_dir + '/' + baxh5file
        writebaxh5(filecontent, baxh5path)
    cmph5file = name_prefix + ".aligned_reads.cmp.h5"

    # baxh5 operations (filter, align(blasr, filter, samtoh5), loadchemistry, loadpulse)
    baxh5_operations = "{baxh5_operations_sh} {seymour_home} {temp_output_folder} {baxh5_filepath}" \
                       " {reference_filepath} {referencesa_filepath} {cmph5_filename} {proc_num}". \
        format(baxh5_operations_sh=baxh5_shell_file_path,
               seymour_home=SMRT_ANALYSIS_HOME,
               temp_output_folder=TEMP_OUTPUT_FOLDER,
               baxh5_filepath=baxh5path,
               reference_filepath=reference_path,
               referencesa_filepath=referencesa_path,
               cmph5_filename=cmph5file,
               proc_num=PROC_NUM)
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
        # write ipdvalue of cmph5 file ---
        if str(GET_IPD_FROM_CMPH5).lower() == 'yes':
            if DATA_SAVE_MODE == 'MASTER':
                get_ipdvalue_of_h5_to_master(cmph5_filepath, CELL_DATA_DIR, "CMPH5",
                                             MAX_SLEEP_SECONDS)
            elif DATA_SAVE_MODE == 'HDFS':
                # hdfs_dir = os.path.dirname(hdfs_filepath)
                hdfs_dir = CELL_DATA_DIR + HDFS_SAMIPDINFO_DIR
                get_ipdvalue_of_h5_to_hdfs(cmph5_filepath, hdfs_dir, "CMPH5")
            elif DATA_SAVE_MODE == 'SHARED_FOLDER':
                get_ipdvalue_of_h5_to_sharedfolder(cmph5_filepath, CELL_DATA_DIR,
                                                   "CMPH5")
            else:
                pass
        # --------------------------------
        # aligned_reads = split_reads_in_cmph5(cmph5_filepath)
        copy_cmph5_to_shared_folder(cmph5_filepath, CMPH5_DATA_DIR)
        # rm temp files --------
        if DATA_SAVE_MODE == 'MASTER' or DATA_SAVE_MODE == 'HDFS':
            os.remove(baxh5path)
        # remove_files_in_a_folder(TEMP_OUTPUT_FOLDER, name_prefix)
        # ----------------------
        # return aligned_reads
        return [1] * 100


def copy_cmph5_to_shared_folder(cmph5file, shared_dir):
    if not os.path.isdir(shared_dir):
        try:
            os.mkdir(shared_dir, 0777)
        except:
            print('shared_dir {} exists.'.format(shared_dir))
    run_cmd_safe(['cp', cmph5file, shared_dir], max_sleep_seconds=200)


def write_ipd_of_cmph5(inCmpH5File, outIpdInfoFile):
    start = time.time()
    inCmpH5 = CmpH5Reader(inCmpH5File)
    issuccess = 0
    try:
        outIpdInfo = open(outIpdInfoFile, 'w')
        outIpdInfo.write("#readID\trefID\trefStart\trefEnd\trefStrand\tMapQV\trefSeq\treadSeq\treadIPD\n")
        for alignment in inCmpH5:
            outIpdInfo.write('\t'.join([alignment.readName, alignment.referenceName,
                                        str(alignment.tStart),
                                        str(alignment.tEnd), str(alignment.RCRefStrand),
                                        str(alignment.mapQV),
                                        alignment.reference(), alignment.read(),
                                        " ".join(map(str, alignment.IPD()))]) + "\n")
        outIpdInfo.flush()
        outIpdInfo.close()
        print("the IPD info of {} has been extracted and written.\ncost {} seconds"
              .format(inCmpH5File, time.time() - start))
        issuccess = 1
    except IOError:
        print("IOError while writing {}".format(outIpdInfoFile))
    finally:
        return issuccess


def writebaxh5(filecontent, filepath):
    """
    a lighter way
    :param filecontent:
    :param filepath:
    :return:
    """
    wf = h5py.File(filepath, "w")
    if isinstance(filecontent, tuple):
        ori_h5_path, h5data = filecontent
    else:
        raise ValueError("the format of h5file info is wrong!")

    rf = h5py.File(ori_h5_path, "r")
    for datasetinfo in h5data:
        datasetname, (lbegin, lend) = datasetinfo
        sdata = rf[datasetname]
        sshape = sdata.shape
        if len(sshape) == 1:
            wf.create_dataset(datasetname, data=sdata[lbegin:lend])
        else:
            wf.create_dataset(datasetname, data=sdata[lbegin:lend, :])

    h5attr = get_baxh5_attrs(rf)
    for itemattrinfo in h5attr:
        item_name_type, itemattrs = itemattrinfo
        if item_name_type[1] == H5GROUP:
            itemtmp = wf.require_group(item_name_type[0])
        else:
            if item_name_type[0] in wf:
                itemtmp = wf[item_name_type[0]]
            else:
                # FIXME: this is not recommended, all data of datasets
                # FIXME: should be written in h5data (the last for loop)
                itemtmp = wf.create_dataset(item_name_type[0], (1,))
        for itemattr_name, itemattr_val in itemattrs:
            itemtmp.attrs[itemattr_name] = itemattr_val
    rf.close()
    wf.close()


def split_reads_in_cmph5(cmph5_filepath):
    """

    :param cmph5_filepath:
    :return: info about each read, [(key, value), ...]
    key: (reffullname, target_start, target_end)
    value: dict(), {'/../AlnIndex': a line of AlnIndex, 'AlnArray': numpy.array,
    'DeletionQV': numpy.array, ...}
    """
    cmph5 = h5py.File(cmph5_filepath, 'r')
    ch5Format = CmpH5Format()

    aI = cmph5[ch5Format.ALN_INDEX]
    aIlen = aI.shape[0]
    refId2refFullName = refGroupId2RefInfoIdentifier(cmph5, ch5Format)
    # refId2refMd5 = refGroupId2RefInfoIdentifier(cmph5, ch5Format, "MD5")
    movieId2MovieInfo = MovieId2MovieInfo(cmph5, ch5Format)

    reads_rdd = [None] * aIlen
    for i in range(0, aIlen):
        alnIndexTmp = aI[i, :]

        read_key = (refId2refFullName[alnIndexTmp[ch5Format.REF_ID]],
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
    MAX_BLOCK_SIZE = max_chunk_length

    # Maximum number of hits per chunk
    MAX_HITS = max_reads_per_chunk
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
    key: (reffullname, (ref_splitting_range_start, ref_splitting_range_end, ref_splitting_range_folds))
    val: [read_value1, read_value2, ...] (read_value is the same as value in read_keyval)
    :param moviechemistry:
    :param refinfo:
    :return: (reffullname, (ref_start, {'csv': csvContent, 'gff': gffContent, }))
    """
    (reffullname, (ref_start, ref_end, ref_folds)) = keyval[0]
    reads_info = [read_keyval[1] for read_keyval in list(keyval[1])]
    name_prefix = rename(reffullname) + '.' + str(ref_start) + '-' + str(ref_end)

    if len(reads_info) > 0:
        if not os.path.isdir(TEMP_OUTPUT_FOLDER):
            try:
                os.mkdir(TEMP_OUTPUT_FOLDER, 0777)
            except:
                print('temp directory {} exists.'.format(TEMP_OUTPUT_FOLDER))

        # setting paths and variables
        contig_filename = name_reference_contig_file(REF_FILENAME, reffullname)
        if DATA_SAVE_MODE == 'MASTER' or DATA_SAVE_MODE == 'HDFS':
            reference_path = SparkFiles.get(contig_filename)
        elif DATA_SAVE_MODE == 'SHARED_FOLDER':
            reference_path = os.path.join(REFERENCE_DIR, contig_filename)
        else:
            reference_path = ''

        cmph5_shell_file_path = SparkFiles.get(shell_script_cmph5)
        cmph5_filename = name_prefix + ".cmp.h5"
        cmph5path = '/'.join([TEMP_OUTPUT_FOLDER, cmph5_filename])

        modification_gff = name_prefix + ".modifications.gff"
        modification_csv = name_prefix + ".modifications.csv"

        refchunkinfo = ','.join(["1", str(ref_start), str(ref_end), str(ref_folds)])
        ipdMaxCoverage = IPDMAXCOVERAGE

        # writing cmph5 file
        writecmph5(cmph5path, reads_info, reffullname, refinfo, moviechemistry)

        # cmph5 operations (sort, repack, computeModifications)
        cmph5_operations = "{cmph5_operations_sh} {seymour_home} {temp_output_folder} {cmph5_filepath}" \
                           " {ref_chunk_info} {reference_filepath} {gff_filename} {csv_filename}" \
                           " {max_coverage} {proc_num} {methylation_type}" \
            .format(cmph5_operations_sh=cmph5_shell_file_path,
                    seymour_home=SMRT_ANALYSIS_HOME,
                    temp_output_folder=TEMP_OUTPUT_FOLDER,
                    cmph5_filepath=cmph5path,
                    ref_chunk_info=refchunkinfo,
                    reference_filepath=reference_path,
                    gff_filename=modification_gff,
                    csv_filename=modification_csv,
                    max_coverage=ipdMaxCoverage,
                    proc_num=SPARK_TASK_CPUS,
                    methylation_type=METHYLATION_TYPES)

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
            gff_filepath = '/'.join([TEMP_OUTPUT_FOLDER, modification_gff])
            csv_filepath = '/'.join([TEMP_OUTPUT_FOLDER, modification_csv])

            gffContent, csvContent = "", ""
            with open(gff_filepath) as rf:
                for line in rf:
                    if not line.startswith("##"):
                        gffContent += line.strip() + '\n'
                        break
                gffContent += "\n".join(line.strip() for line in rf) + '\n'
            with open(csv_filepath) as rf:
                next(rf)
                csvContent += "\n".join(line.strip() for line in rf) + '\n'
            # rm temp files --------
            remove_files_in_a_folder(TEMP_OUTPUT_FOLDER, name_prefix)
            # ----------------------
            return reffullname, (ref_start, {'csv': csvContent, 'gff': gffContent, })
    else:
        raise ValueError("no reads to form a cmph5 file")


def writecmph5(filepath, reads_info, reffullname, refinfo, moviechemistry):
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
        alnindex[i, f_format.ID] = i + 1
        alnindex[i, f_format.ALN_ID] = CMPH5_ALN_GROUP_ID
        alnindex[i, f_format.MOLECULE_ID] = i
        alnindex[i, f_format.OFFSET_BEGIN] = alignseq_range[i]
        alnindex[i, f_format.OFFSET_END] = alignseq_range[i + 1]

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
            mcontemp += (reads_info[i][movieinfo_items[j]],)
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
                         data=range(1, len(movienames) + 1))
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
                     data=[refinfo[reffullname][SEQUENCE_LENGTH], ])
    f.create_dataset("RefInfo/MD5", dtype=h5py.special_dtype(vlen=str),
                     data=[refinfo[reffullname][SEQUENCE_MD5], ])
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


def create_redundant_reads(read_keyval, ref_splitting_info,
                           danger_atomchunks, atomchunk2enlargedchunk):
    """
    check if a read needs to be copyed for two ref chunks,
    reads->atomchunks->enlargedchunks
    :param read_keyval:
    :param ref_splitting_info:
    :param danger_atomchunks: dict: {atomchunk1: ratio*100, ...}
    :param atomchunk2enlargedchunk:
    :return: [(key, val), ]
    key: (reffullname, (ref_splitting_range_start, ref_splitting_range_end, ref_splitting_range_folds))
    val: the same as value in read_keyval
    """
    # FIXED: classify the reads by covering ref_splitting_range (start, end)
    # FIXED: or (start-pad, end+pad)?. For now, use (start, end)
    # FIXED: answer: use (start, end) is enough
    (read_refid, target_start, target_end) = read_keyval[0]
    read_ref_splitting_info = ref_splitting_info[read_refid]

    res_keyval = []
    enlarged_chunks = set()
    for ref_range in read_ref_splitting_info:
        refstart, refend = ref_range[0], ref_range[1]
        # refstart, refend = refstart - PAD, refend + PAD
        if not (target_end <= refstart or target_start >= refend):
            ref_chunk = (read_refid, ref_range)
            if ref_chunk in danger_atomchunks.keys():
                if keepornot(danger_atomchunks[ref_chunk]):
                    enlarged_chunks.add(atomchunk2enlargedchunk[ref_chunk])
            else:
                enlarged_chunks.add(atomchunk2enlargedchunk[ref_chunk])
    for echunk in enlarged_chunks:
        res_keyval.append((echunk, read_keyval))
    return res_keyval


def keepornot(keep_ratio):
    randomNum = random.randint(1, 100)
    return 0 if randomNum > keep_ratio else 1


def get_atomchunks(read_key, ref_splitting_info):
    """
    check if a read needs to be copyed for two ref chunks
    :param read_key:
    :param ref_splitting_info:
    :return:
    """
    # FIXED: classify the reads by covering ref_splitting_range (start, end)
    # FIXED: or (start-pad, end+pad)?. For now, use (start, end)
    # FIXED: answer: use (start, end) is enough
    (read_refid, target_start, target_end) = read_key
    read_ref_splitting_info = ref_splitting_info[read_refid]

    chunksOfThisRead = []
    for ref_range in read_ref_splitting_info:
        refstart, refend = ref_range[0], ref_range[1]
        # refstart, refend = refstart - PAD, refend + PAD
        if not (target_end <= refstart or target_start >= refend):
            chunksOfThisRead.append((read_refid, ref_range))
    return chunksOfThisRead


def adjust_ref_splitting_info(ref_splitting_info, factor=2):
    factor = int(factor)
    if factor == 0:
        factor = 1
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
            sbegin, send = chunk_steps_range[i], chunk_steps_range[i + 1]

            ref_splitting_adjustedrange.append((ref_splitting_range[sbegin][0],
                                                ref_splitting_range[send - 1][1],
                                                send - sbegin))

        ref_splitting_info[refid] = ref_splitting_adjustedrange
    return ref_splitting_info


# def remove_partofreads_from_repeats(refchunk2reads):
#     trim_strategy = 'random'  # 'mapqv' or 'random'
#
#     strategies = ['random', 'mapqv']
#     if READS_TRIM_STRATEGY in strategies:
#         trim_strategy = READS_TRIM_STRATEGY
#
#     reads = list(refchunk2reads[1])
#     if len(reads) > limitation_readsnum:
#         if trim_strategy == 'mapqv':
#             ALN_INDEX = 'AlnInfo/AlnIndex'
#             MAPQV = 13
#             mapqvs = [read[1][ALN_INDEX][MAPQV] for read in reads]
#             sorted_idxs = [a[0] for a in sorted(enumerate(mapqvs), key=lambda x: x[1], reverse=True)]
#
#             keeped_reads = [reads[a] for a in sorted_idxs[0:limitation_readsnum]]
#             return refchunk2reads[0], keeped_reads
#         elif trim_strategy == 'random':
#             return refchunk2reads[0], random.sample(reads, limitation_readsnum)
#         else:
#             return refchunk2reads[0], random.sample(reads, limitation_readsnum)
#     else:
#         return refchunk2reads[0], reads


# FIXME: need to re-code this function, make it more pythonic
def resplit_refchunks(atomchunk2readsnum, max_chunknum):
    """

    :param atomchunk2readsnum: [(refinfo<tuple>, readsnum<int>), (), ]
    :param max_chunknum:
    :return:
    """
    atomlen = len(atomchunk2readsnum)
    if atomlen == 0:
        return None
    rechunks_info = []
    readnum_sum_tmp = atomchunk2readsnum[0][1]
    perchunk_tmp = [atomchunk2readsnum[0][0]]
    for i in range(1, atomlen):
        if atomchunk2readsnum[i][0][0] != atomchunk2readsnum[i-1][0][0]:
            if len(perchunk_tmp) > 0:
                rechunks_info.append(perchunk_tmp)
            perchunk_tmp = [atomchunk2readsnum[i][0]]
            readnum_sum_tmp = atomchunk2readsnum[i][1]
        else:
            if len(perchunk_tmp) >= max_chunknum:
                rechunks_info.append(perchunk_tmp)
                readnum_sum_tmp = atomchunk2readsnum[i][1]
                perchunk_tmp = [atomchunk2readsnum[i][0]]
                continue
            readnum_sum_tmp += atomchunk2readsnum[i][1]
            if readnum_sum_tmp >= limitation_readsnum:
                if readnum_sum_tmp - atomchunk2readsnum[i][1] > 0:
                    rechunks_info.append(perchunk_tmp)
                if atomchunk2readsnum[i][1] >= limitation_readsnum:
                    rechunks_info.append([atomchunk2readsnum[i][0]])
                    readnum_sum_tmp = 0
                    perchunk_tmp = []
                else:
                    readnum_sum_tmp = atomchunk2readsnum[i][1]
                    perchunk_tmp.append(atomchunk2readsnum[i][0])
            else:
                perchunk_tmp.append(atomchunk2readsnum[i][0])
    # check when for-loop is done
    if len(perchunk_tmp) > 0:
        rechunks_info.append(perchunk_tmp)
    atomchunk2enlargedchunk = {}
    for rechunks in rechunks_info:
        en_reffullname = rechunks[0][0]
        rstart = min(map(lambda x: x[1][0], rechunks))
        rend = max(map(lambda x: x[1][1], rechunks))
        folds = numpy.sum(map(lambda x: x[1][2], rechunks))
        enlargedchunk = (en_reffullname, (rstart, rend, folds))
        for rechunk in rechunks:
            atomchunk2enlargedchunk[rechunk] = enlargedchunk
    return atomchunk2enlargedchunk


# def remove_redundant_reads(enlargedchunk2readsinfo):
#     """
#
#     :param enlargedchunk2readsinfo: (enchunkinfo, [(chunkinfo, [read, read,]), ()])
#     :return: (enchunkinfo, [read, read,...])
#     """
#     chunkinfo = enlargedchunk2readsinfo[0]
#     readsgroup = sorted(list(enlargedchunk2readsinfo[1]), key=lambda x: x[0][1][0])
#     totalreads = []
#     totalreads.extend(readsgroup[0][1])
#     for i in range(1, len(readsgroup)):
#         (rstart, rend, rfold) = readsgroup[i][0][1]
#         for read in readsgroup[i][1]:
#             if read[0][1] >= rstart:
#                 totalreads.append(read)
#     del readsgroup
#     return chunkinfo, totalreads


# operations for each rdd element in modificationRDD-----------------------
def basemods_pipeline_modification_operations(keyval, refinfo):
    reffullname = keyval[0]
    modsinfo = list(keyval[1])

    if not os.path.isdir(TEMP_OUTPUT_FOLDER):
        try:
            os.mkdir(TEMP_OUTPUT_FOLDER, 0777)
        except:
            print('temp directory {} exists.'.format(TEMP_OUTPUT_FOLDER))

    # setting paths and variables
    # reference_path = SparkFiles.get(REF_FILENAME)
    contig_filename = name_reference_contig_file(REF_FILENAME, reffullname)
    reference_path = SparkFiles.get(contig_filename)
    mods_shell_file_path = SparkFiles.get(shell_script_mods)
    name_prefix = rename(reffullname)
    gfffilename = name_prefix + ".modifications.gff"
    csvfilename = name_prefix + ".modifications.csv"
    motifs_gff_gz_filename = name_prefix + ".motifs.gff.gz"

    gfffilepath = '/'.join([TEMP_OUTPUT_FOLDER, gfffilename])
    csvfilepath = '/'.join([TEMP_OUTPUT_FOLDER, csvfilename])
    writemodificationinfo(modsinfo, reffullname, refinfo, gfffilepath, csvfilepath)

    # modification operations (findMotifs, makeMotifGff)
    mods_operations = "{mods_operations_sh} {seymour_home} {temp_output_folder} {gff_filepath}" \
                      " {csv_filepath} {reference_filepath} {motifsgffgz_filename}" \
        .format(mods_operations_sh=mods_shell_file_path,
                seymour_home=SMRT_ANALYSIS_HOME,
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
                 "##sequence-region ref000001 1 {}\n".format(refinfo[reffullname][SEQUENCE_LENGTH])
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


def get_ip_of_node():
    node_ip = [l for l in ([ip for ip in socket.gethostbyname_ex(socket.gethostname())[2]
                            if not ip.startswith("127.")][:1], [[(s.connect(('8.8.8.8', 53)),
                                                                  s.getsockname()[0], s.close())
                                                                 for s in [socket.socket(socket.AF_INET,
                                                                                         socket.SOCK_DGRAM)]][
                                                                    0][1]])
               if l][0][0]
    return node_ip


def writemods_of_each_chromosome(keyval, refinfo, max_sleep_seconds=1):
    reffullname = keyval[0]
    modsinfo = list(keyval[1])

    if not os.path.isdir(TEMP_OUTPUT_FOLDER):
        try:
            os.mkdir(TEMP_OUTPUT_FOLDER, 0777)
        except:
            print('temp directory {} exists.'.format(TEMP_OUTPUT_FOLDER))

    name_prefix = rename(reffullname)
    gfffilename = name_prefix + ".modifications.gff"
    csvfilename = name_prefix + ".modifications.csv"

    gfffilepath = '/'.join([TEMP_OUTPUT_FOLDER, gfffilename])
    csvfilepath = '/'.join([TEMP_OUTPUT_FOLDER, csvfilename])
    writemodificationinfo(modsinfo, reffullname, refinfo, gfffilepath, csvfilepath)

    if DATA_SAVE_MODE == 'MASTER':
        master_ip = MASTERNODE_IP
        master_port = MASTERNODE_PORT
        master_user = MASTERNODE_USERNAME
        master_passwd = MASTERNODE_USERPASSWD

        mgfffilepath = '/'.join([CELL_DATA_DIR, gfffilename])
        mcsvfilepath = '/'.join([CELL_DATA_DIR, csvfilename])
        worker_put_master(master_ip, master_port, master_user,
                          master_passwd, gfffilepath, mgfffilepath,
                          max_sleep_seconds)
        worker_put_master(master_ip, master_port, master_user,
                          master_passwd, csvfilepath, mcsvfilepath,
                          max_sleep_seconds)
    elif DATA_SAVE_MODE == 'HDFS':
        hdfs_modsresult_dir = CELL_DATA_DIR + HDFS_MODS_DIR
        mgfffilepath = '/'.join([hdfs_modsresult_dir, gfffilename])
        mcsvfilepath = '/'.join([hdfs_modsresult_dir, csvfilename])
        cmd_output, cmd_errors = run_cmd_safe([HDFS_CMD, 'dfs', '-copyFromLocal', '-f',
                                               gfffilepath,
                                               mgfffilepath])
        cmd_output, cmd_errors = run_cmd_safe([HDFS_CMD, 'dfs', '-copyFromLocal', '-f',
                                               csvfilepath,
                                               mcsvfilepath])
    elif DATA_SAVE_MODE == 'SHARED_FOLDER':
        mgfffilepath = '/'.join([CELL_DATA_DIR, gfffilename])
        mcsvfilepath = '/'.join([CELL_DATA_DIR, csvfilename])
        run_cmd(['cp', gfffilepath, mgfffilepath])
        run_cmd(['cp', csvfilepath, mcsvfilepath])
    else:
        mgfffilepath = ''
        mcsvfilepath = ''
    return mgfffilepath, mcsvfilepath


# write each contig in a ref_fasta to a single file---------------
def write_ref_contigs(ref_dir, ref_name, ref_contigs):
    """

    :param ref_dir:
    :param ref_name:
    :param ref_contigs:
    :return:
    """
    contig_filepaths = []
    for contigkey in ref_contigs.keys():
        contigname = name_reference_contig_file(ref_name, contigkey)
        contigpath = '/'.join([ref_dir, contigname])
        with open(contigpath, 'w') as wf:
            wf.write('>' + contigkey + '\n')
            sequence = ref_contigs[contigkey][SEQUENCE]
            seqLength = ref_contigs[contigkey][SEQUENCE_LENGTH]
            xsteps = [s for s in range(0, seqLength, COLUMNS)]
            for i in range(0, len(xsteps)-1):
                wf.write(sequence[xsteps[i]:xsteps[i+1]] + '\n')
            wf.write(sequence[xsteps[-1]:seqLength] + '\n')
            del sequence
        contig_filepaths.append(contigpath)
    return contig_filepaths


# to clear temp folder in each worker node---------------
def rm_temp_folder(temp_folder):
    issuccess = 0
    if os.path.isdir(temp_folder):
        try:
            shutil.rmtree(temp_folder)
            issuccess = 1
            print("temp folder of {} has been deleted.".format(get_ip_of_node()))
        except OSError:
            print("something is wrong when deleting temp folder of {}, but don't worry.".format(get_ip_of_node()))
    return issuccess


# run_cmd--------------------------------------------------
# def run_cmd(args_list):
#     print('Running system command: {0}'.format(' '.join(args_list)))
#     proc = Popen(args_list, stdout=PIPE, stderr=PIPE)
#     (output, errors) = proc.communicate()
#     if proc.returncode:
#         raise RuntimeError(
#             'Error running command: %s. Return code: %d, Error: %s' % (
#                 ' '.join(args_list), proc.returncode, errors))
#     return output, errors

def run_cmd(args_list):
    proc = Popen(args_list, stdout=PIPE, stderr=PIPE)
    proc.communicate()
    return proc.returncode


def run_cmd_safe(args_list, max_sleep_seconds=1):
    print('Running system command: {0}'.format(' '.join(args_list)))
    attemp_times = 100
    for i in range(0, attemp_times):
        expected_sleep_seconds = random.randint(0, max_sleep_seconds) * (i + 1)
        actual_sleep_seconds = expected_sleep_seconds \
            if expected_sleep_seconds < max_sleep_seconds else max_sleep_seconds
        time.sleep(actual_sleep_seconds)
        proc = Popen(args_list, stdout=PIPE, stderr=PIPE)
        (output, errors) = proc.communicate()
        if not proc.returncode:
            print('Running system command: {0}, succeed!'.format(' '.join(args_list)))
            return output, errors
    print('Running system command: {0}, failed!'.format(' '.join(args_list)))
    return '', 'error'


def run_hdfs_get_cmd(hdfs_file, local_file, max_sleep_seconds=1):
    if os.path.isfile(local_file):
        os.remove(local_file)
    return run_cmd_safe([HDFS_CMD, 'dfs', '-get',
                         hdfs_file,
                         local_file],
                        max_sleep_seconds)


def mkdir_in_hdfs(thedir):
    hdfs_dir = thedir
    isFileOrDir = run_cmd([HDFS_CMD, 'dfs', '-test', '-e', hdfs_dir])  # 0 is true, 1 is false
    isDir = run_cmd([HDFS_CMD, 'dfs', '-test', '-d', hdfs_dir])  # 0 is true, 1 is false
    if isFileOrDir == 1:
        cmd_output, cmd_errors = run_cmd_safe([HDFS_CMD, 'dfs', '-mkdir', '-p',
                                               hdfs_dir])
    elif isDir == 1:
        raise ValueError("'{}' in your HDFS is a file, not a directory.\n"
                         "however, we need it to be a directory.\n"
                         "please fix it".format(hdfs_dir))


# rm files with specific prefix-pattern in a folder ----------
def remove_files_in_a_folder(file_directory, file_regex):
    for f in os.listdir(file_directory):
        if f.startswith(file_regex):
            os.remove('/'.join([file_directory, f]))
# ------------------------------------------------------------------------------------------


def basemods_pipe():
    pipe_start = time.time()

    # get global variables------------------------------------------------------------------
    # for reference .sa file name
    global USED_REF_SA_FILENAME
    abs_dir = os.path.dirname(os.path.realpath(__file__))
    getParametersFromFile('/'.join([abs_dir, parameters_config]))

    # start SparkContext--------------------------------------------------------------------
    SparkContext.setSystemProperty('spark.executor.memory', SPARK_EXECUTOR_MEMORY)
    SparkContext.setSystemProperty('spark.task.cpus', SPARK_TASK_CPUS)
    SparkContext.setSystemProperty('spark.memory.fraction', SPARK_MEMORY_FRACTION)
    SparkContext.setSystemProperty('spark.memory.storageFraction', SPARK_MEMORY_STORAGEFRACTION)
    conf = SparkConf().setAppName("Spark-based Pacbio BaseMod pipeline")
    sc = SparkContext(conf=conf)

    # create temp folder in master node-----------------------------------------------------
    if not os.path.isdir(TEMP_OUTPUT_FOLDER):
        try:
            os.mkdir(TEMP_OUTPUT_FOLDER, 0777)
        except:
            print('temp directory {} exists.'.format(TEMP_OUTPUT_FOLDER))

    # files need to be shared to each node--------------------------------------------------
    # reference----
    ref_dir = REFERENCE_DIR
    if DATA_SAVE_MODE == 'HDFS':
        local_ref_file = '/'.join([TEMP_OUTPUT_FOLDER, REF_FILENAME])
        cmd_output, cmd_errors = run_hdfs_get_cmd('/'.join([REFERENCE_DIR, REF_FILENAME]),
                                                  local_ref_file)
        sc.addFile(local_ref_file)
        ref_dir = TEMP_OUTPUT_FOLDER
    elif DATA_SAVE_MODE == 'MASTER':
        sc.addFile('/'.join([REFERENCE_DIR, REF_FILENAME]))
        ref_dir = REFERENCE_DIR
    elif DATA_SAVE_MODE == 'SHARED_FOLDER':
        ref_dir = REFERENCE_DIR
    else:
        print('please set DATA_SAVE_MODE in parameters.conf')
        return
    # reference.sa----
    if REF_SA_FILENAME == 'None' or (not REF_SA_FILENAME.endswith('.sa')):
        print('exec sawriter (blasr) to generate a .sa file for your reference.')
        sa_script_path = '/'.join([abs_dir, 'scripts', shell_script_sa])
        ref_sa_filename = exec_sawriter(sa_script_path, '/'.join([ref_dir, REF_FILENAME]))
        print('sawriter finished. {} is generated'.format(ref_sa_filename))
        if DATA_SAVE_MODE == 'MASTER' or DATA_SAVE_MODE == 'HDFS':
            sc.addFile('/'.join([ref_dir, ref_sa_filename]))
        USED_REF_SA_FILENAME = ref_sa_filename
    else:
        if DATA_SAVE_MODE == 'HDFS':
            local_refsa_file = '/'.join([TEMP_OUTPUT_FOLDER, REF_SA_FILENAME])
            cmd_output, cmd_errors = run_hdfs_get_cmd('/'.join([REFERENCE_DIR, REF_SA_FILENAME]),
                                                      local_refsa_file)
            sc.addFile(local_refsa_file)
        elif DATA_SAVE_MODE == 'MASTER':
            sc.addFile('/'.join([REFERENCE_DIR, REF_SA_FILENAME]))
        else:
            pass
        USED_REF_SA_FILENAME = REF_SA_FILENAME
    # ref contigs----
    refcontigs = getRefInfoFromFastaFiles(['/'.join([ref_dir, REF_FILENAME]), ])
    contig_filepaths = write_ref_contigs(ref_dir, REF_FILENAME, refcontigs)
    if DATA_SAVE_MODE == 'HDFS' or DATA_SAVE_MODE == 'MASTER':
        for contig_filepath in contig_filepaths:
            sc.addFile(contig_filepath)
    else:
        pass
    # del sequence
    for contig in refcontigs.keys():
        refcontigs[contig][SEQUENCE] = ''
    # scipts----
    sc.addFile('/'.join([abs_dir, 'scripts', shell_script_baxh5]))
    sc.addFile('/'.join([abs_dir, 'scripts', shell_script_cmph5]))
    sc.addFile('/'.join([abs_dir, 'scripts', shell_script_mods]))

    # STEP 1 baxh5->cmph5 (filter, blasr)---------------------------------------------------
    # get all files in cell_data_directory
    pacbio_data_dir = CELL_DATA_DIR
    baxh5_filenames = []
    metaxml_filenames = []
    if DATA_SAVE_MODE == 'MASTER' or DATA_SAVE_MODE == 'SHARED_FOLDER':
        for root, dirnames, filenames in os.walk(pacbio_data_dir):
            for filename in fnmatch.filter(filenames, '*.bax.h5'):
                baxh5_filenames.append(os.path.join(root, filename))
            for filename in fnmatch.filter(filenames, '*.metadata.xml'):
                metaxml_filenames.append(os.path.join(root, filename))
    elif DATA_SAVE_MODE == 'HDFS':
        cmd_output, cmd_errors = run_cmd_safe([HDFS_CMD, 'dfs', '-ls', '-R',
                                               pacbio_data_dir])
        fileitems = [word[-1] for word in [words.split(' ') for words in cmd_output.split('\n')]]
        # for baxh5_filenames and metadata.xml----
        for fileitem in fileitems:
            if fileitem.endswith('.bax.h5'):
                baxh5_filenames.append(fileitem)
            if fileitem.endswith('.metadata.xml'):
                filename = os.path.basename(fileitem)
                local_metadataxml_file = '/'.join([TEMP_OUTPUT_FOLDER, filename])
                cmd_output, cmd_errors = run_hdfs_get_cmd(fileitem, local_metadataxml_file)
                metaxml_filenames.append(local_metadataxml_file)

        # create folder in HDFS for saving mods result----
        mkdir_in_hdfs(CELL_DATA_DIR + HDFS_MODS_DIR)
        # create fastipd folder in HDFS if needed----
        if str(GET_IPD_FROM_BASH5).lower() == 'yes':
            mkdir_in_hdfs(pacbio_data_dir + HDFS_IPDINFO_DIR)
        # create samipd folder in HDFS if needed----
        if str(GET_IPD_FROM_CMPH5).lower() == 'yes':
            mkdir_in_hdfs(pacbio_data_dir + HDFS_SAMIPDINFO_DIR)
    else:
        pass

    baxh5_folds = BAXH5_FOLDS
    local_temp_dir = TEMP_OUTPUT_FOLDER
    # FIXME: how to do it smarter?---------------
    shuffle_factor = 1000
    numpartitions = len(baxh5_filenames) * shuffle_factor
    re_numpartitions = numpartitions if numpartitions < max_numpartitions else max_numpartitions
    baxh5nameRDD = sc.parallelize(baxh5_filenames, len(baxh5_filenames))\
        .repartition(re_numpartitions)\
        .coalesce(len(baxh5_filenames))

    if DATA_SAVE_MODE == 'MASTER':
        aligned_reads_rdd = baxh5nameRDD. \
            map(lambda x: get_baxh5file_from_masternode(x, local_temp_dir, MAX_SLEEP_SECONDS)). \
            flatMap(lambda x: get_chunks_of_baxh5file(x, baxh5_folds)). \
            flatMap(basemods_pipeline_baxh5_operations). \
            persist(StorageLevel.DISK_ONLY)  # use DISK_ONLY temporarily
            #persist(StorageLevel.MEMORY_AND_DISK)
    elif DATA_SAVE_MODE == 'HDFS':
        aligned_reads_rdd = baxh5nameRDD. \
            map(lambda x: get_baxh5file_from_hdfs(x, local_temp_dir)). \
            flatMap(lambda x: get_chunks_of_baxh5file(x, baxh5_folds)). \
            flatMap(basemods_pipeline_baxh5_operations). \
            persist(StorageLevel.DISK_ONLY)  # use DISK_ONLY temporarily
            #persist(StorageLevel.MEMORY_AND_DISK)
    elif DATA_SAVE_MODE == 'SHARED_FOLDER':
        aligned_reads_rdd = baxh5nameRDD. \
            map(lambda x: get_baxh5file_from_sharedfolder(x)). \
            flatMap(lambda x: get_chunks_of_baxh5file(x, baxh5_folds)). \
            flatMap(basemods_pipeline_baxh5_operations). \
            persist(StorageLevel.DISK_ONLY)  # use DISK_ONLY temporarily
            # persist(StorageLevel.MEMORY_AND_DISK)
    else:
        aligned_reads_rdd = None
        return

    # STEP 2 cmph5->mods.gff/csv (ipdSummary.py)--------------------------------------------
    # x[0] is ref_contig's fullname, check split_reads_in_cmph5()
    # type(ref_identifiers_count):dict
    # ref_identifiers_count = aligned_reads_rdd.map(lambda (x, y): x[0]).countByValue()
    aligned_reads_rdd.count()

    # clear persisted rdd and broadcast variables----------------------------------------------
    aligned_reads_rdd.unpersist()

    # # rm temp folder of each worker node and master node---------------------------------------
    # # can't guarantee rm every worker's temp folder
    # worker_num = 40
    # rdd_ele_num = worker_num * 10
    # rm_num = sc.range(rdd_ele_num, numSlices=rdd_ele_num)\
    #     .map(lambda x: TEMP_OUTPUT_FOLDER)\
    #     .map(rm_temp_folder)\
    #     .reduce(lambda x, y: x + y)
    # print("temp folders of {} worker node(s) have been deleted.".format(rm_num))
    # missuccess = rm_temp_folder(TEMP_OUTPUT_FOLDER)
    # if missuccess > 0:
    #     print("temp folder of master node has been deleted.")

    # exit-------------------------------------------------------------------------------------
    SparkContext.stop(sc)
    print('total time cost: {} seconds'.format(time.time() - pipe_start))


if __name__ == '__main__':
    print("spark start------------------------------------")
    basemods_pipe()
    print("spark end--------------------------------------\n"
          "-----------------------------------------------")
