## reading HDF5 File in Spark

### Solution 1. Use the master node's memory

#### description
> 1. read each dataset in an H5 file to master node's memory.
> 2. use *sc.parallelize()* to distribute each dataset to worker node.

#### limitations
> 1. It takes huge memory of master node. It will not work when there are too many/large H5 files.

#### tests
> 1. see *basemods_spark_runner_memory.py*


### Solution 2. use a storage system/Network File System

#### description
> build a NFS in the cluster, so that each worker node can read data from the same directory and we only need one copy of the data.

#### limitations
> For now we don't know how to build a NFS, and may be this solution is not suitable for all situations.

#### tests
> Hint: if you do some tests, please log the time, the result of your tests and how you did it.


### Solution 2. use a storage system/Network File System

#### description
> build a NFS in the cluster, so that each worker node can read data from the same directory and we only need one copy of the data.

#### limitations
> For now we don't know how to build a NFS, and may be this solution is not suitable for all situations.

#### tests
> Hint: if you do some tests, please log the time, the result of your tests and how you did it.

### Solution 3. Copy all the data to each worker node

#### description
> Copy all the data to the same directory in each worker node.

#### limitations
> It will make the job a lot easier, but it seems that it is not smart enough. It costs a great deal of space of hard disks.

#### tests
> 


### Solution 4. Transforming H5 file format 

#### description
> Transform H5 file to a hadoop-friendly file format (TextFile, SquenceFile, HadoopInputFileFormat) first, then save the transformed files to HDFS/AWS3 etc. The transformation __can/need__ be done in a single machine or a spark cluster. 

#### limitations
> The transformation may cost a lot of time. Need to find an efficient way.

#### tests
> Hint: if you do some tests, please log the time, the result of your tests and how you did it.

