#!/usr/bin/python
# coding=utf-8
from pyspark import SparkContext, SparkConf
from subprocess import Popen, PIPE
import os
import random
import time
import socket

worker_num = 2
HDFS_CMD = "/usr/local/hadoop/bin/hdfs"
TEMP_OUTPUT_FOLDER = "/tmp"
smrtanalysistar_path = "hdfs://127.0.0.1:9000/tools/smrtanalysis.tar.gz"
scripts_path = "hdfs://127.0.0.1:9000/tools/scripts"


def run_cmd(args_list):
    print('Running system command: {0}'.format(' '.join(args_list)))
    proc = Popen(args_list, stdout=PIPE, stderr=PIPE)
    (output, errors) = proc.communicate()
    if proc.returncode:
        output, errors = "wrong", "wrong"
        print("something is wrong")
    return output, errors


def get_ip_of_node():
    node_ip = [l for l in ([ip for ip in socket.gethostbyname_ex(socket.gethostname())[2]
                            if not ip.startswith("127.")][:1], [[(s.connect(('8.8.8.8', 53)),
                                                                  s.getsockname()[0], s.close())
                                                                 for s in [socket.socket(socket.AF_INET,
                                                                                         socket.SOCK_DGRAM)]][
                                                                    0][1]])
               if l][0][0]
    return node_ip


# used for mapPartitions
def get_files_from_hdfs(something):
    time.sleep(random.randint(0, 120))
    scripts_localpath = '/'.join([TEMP_OUTPUT_FOLDER, os.path.basename(scripts_path)])
    smrtanalysistar_localpath = '/'.join([TEMP_OUTPUT_FOLDER,
                                          os.path.basename(smrtanalysistar_path)])
    worker_ip = get_ip_of_node()
    resstr = ""

    # get scripts from hdfs and chmod
    if os.path.isdir(scripts_localpath):
        resstr += "scripts folder has already existed.\n"
    else:
        o1, e1 = run_cmd([HDFS_CMD, 'dfs', '-get', scripts_path, scripts_localpath])
        if o1 == "wrong" and e1 == "wrong":
            resstr += "MAYBE scripts folder has already existed.\n"
        else:
            run_cmd(['chmod', '-R', '+x', scripts_localpath])
            resstr += "the scripts folder has been got from HDFS successfully!\n"

    # get smrtanalysis.tar.gz from hdfs and uncompress it
    if (os.path.isfile(smrtanalysistar_localpath) or
            os.path.isfile(smrtanalysistar_localpath + '._COPYING_')):
        resstr += "smrtanalysis.tar.gz has already existed.\n"
        if os.path.isdir(smrtanalysistar_localpath.split('.tar.gz')[0]):
            resstr += "smrtanalysis folder has already existed.\n"
        else:
            pass
            # run_cmd(['tar', '-xhzvf', smrtanalysistar_localpath, '-C', TEMP_OUTPUT_FOLDER])
            # resstr += "smrtanalysis folder has beed created from tar.gz file.\n"
    else:
        o1, e1 = run_cmd([HDFS_CMD, 'dfs', '-get', smrtanalysistar_path, smrtanalysistar_localpath])
        if o1 == "wrong" and e1 == "wrong":
            resstr += "MAYBE smrtanalysis.tar.gz has already existed.\n"
        else:
            run_cmd(['tar', '-xhzvf', smrtanalysistar_localpath, '-C', TEMP_OUTPUT_FOLDER])
            resstr += "smrtanalysis folder has beed created from tar.gz file " \
                      "after getting tar.gz file from HDFS.\n"
    return [str(worker_ip) + ":\n" + resstr]


if __name__ == '__main__':
    pipe_start = time.time()

    SparkContext.setSystemProperty('spark.executor.instances', str(worker_num))
    conf = SparkConf().setAppName("distribute SMRT-Analysis and scripts")
    sc = SparkContext(conf=conf)

    shuffle_factor = 1000
    elesnum = worker_num * shuffle_factor
    issuccess = sc.range(elesnum, numSlices=elesnum)\
        .coalesce(worker_num)\
        .mapPartitions(get_files_from_hdfs)\
        .collect()

    for iss in issuccess:
        print(iss + '\n\n')

    sc.stop()
    print('total time cost: {} seconds'.format(time.time() - pipe_start))
