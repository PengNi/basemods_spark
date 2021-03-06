#!/usr/bin/python
# coding=utf-8
from pyspark import SparkContext, SparkConf
import os
import time
import socket
import shutil

worker_num = 40
TEMP_OUTPUT_FOLDERorFILE = '/tmp/basemods_spark_data'


def get_ip_of_node():
    node_ip = [l for l in ([ip for ip in socket.gethostbyname_ex(socket.gethostname())[2]
                            if not ip.startswith("127.")][:1], [[(s.connect(('8.8.8.8', 53)),
                                                                  s.getsockname()[0], s.close())
                                                                 for s in [socket.socket(socket.AF_INET,
                                                                                         socket.SOCK_DGRAM)]][
                                                                    0][1]])
               if l][0][0]
    return node_ip


# to clear temp folder or file in each worker node---------------
def rm_temp_folderorfile(temp_folderorfile):
    issuccess = 0
    if os.path.isdir(temp_folderorfile):
        try:
            shutil.rmtree(temp_folderorfile)
            issuccess = 1
            print("temp folder of {} has been deleted.".format(get_ip_of_node()))
        except OSError:
            print("something is wrong when deleting temp folder of {}, but don't worry.".format(get_ip_of_node()))
    elif os.path.isfile(temp_folderorfile):
        try:
            os.remove(temp_folderorfile)
            issuccess = 1
            print("temp file of {} has been deleted.".format(get_ip_of_node()))
        except OSError:
            print("something is wrong when deleting temp file of {}, but don't worry.".format(get_ip_of_node()))
    return issuccess


if __name__ == '__main__':
    pipe_start = time.time()

    SparkContext.setSystemProperty('spark.executor.instances', str(worker_num))
    conf = SparkConf().setAppName("clear tmp data of SMRT-Analysis spark pipeline")
    sc = SparkContext(conf=conf)

    # rm temp folder of each worker node and master node---------------------------------------
    # can't guarantee rm every worker's temp folder
    rdd_ele_num = worker_num * 10
    rm_num = sc.range(rdd_ele_num, numSlices=rdd_ele_num) \
        .map(lambda x: TEMP_OUTPUT_FOLDERorFILE) \
        .map(rm_temp_folderorfile) \
        .reduce(lambda x, y: x + y)
    print("temp folders of {} worker node(s) have been deleted.".format(rm_num))
    missuccess = rm_temp_folderorfile(TEMP_OUTPUT_FOLDERorFILE)
    if missuccess > 0:
        print("temp folder of master node has been deleted.")

    sc.stop()
    print('total time cost: {} seconds'.format(time.time() - pipe_start))
