#!/usr/bin/python
# coding=utf-8
from pyspark import SparkContext, SparkConf
from subprocess import Popen, PIPE
import time
import socket

worker_num = 40
TEMP_OUTPUT_FOLDER = '/tmp'


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


# detect a folder's size---------------
def diskusage(something):
    o1, e1 = run_cmd(['du', '-sh', TEMP_OUTPUT_FOLDER])
    return [str(get_ip_of_node()) + ':\n' + o1 + '\n' + e1]


if __name__ == '__main__':
    pipe_start = time.time()

    SparkContext.setSystemProperty('spark.executor.instances', str(worker_num))
    conf = SparkConf().setAppName("detect diskusage of SMRT-Analysis spark pipeline")
    sc = SparkContext(conf=conf)

    # rm temp folder of each worker node and master node---------------------------------------
    # can't guarantee rm every worker's temp folder
    rdd_ele_num = worker_num * 10
    ele_disk_usage = sc.range(rdd_ele_num, numSlices=rdd_ele_num) \
        .coalesce(worker_num)\
        .mapPartitions(diskusage) \
        .collect()
    for ele in ele_disk_usage:
        print(ele)

    sc.stop()
    print('total time cost: {} seconds'.format(time.time() - pipe_start))
