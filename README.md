## basemods on spark

### Overview
the spark version of basemods pipeline in SMRT-Analysis (v2.3.0)

### Set up the environment

The OSs must be **Linux**.

1. #### Hadoop/Spark

    Setting up an [Hadoop](http://hadoop.apache.org/)/[Spark](https://spark.apache.org/) cluster.

2. #### SMRT-Analysis

    For now, we are using SMRT-Analysis v2.3.0

    2.1 Download [smrtanalysis.tar.gz](https://1drv.ms/u/s!AgfGWBktzWTwgjc7p4vgxt15FPQE).
    
    2.2 Copy _smrtanalysis.tar.gz_ to **each worker node** of your Spark cluster. Then decompress it to a desired directory.
    ```sh
    # suppose you want to decompress smrtanalysis.tar.gz to /home/hadoop, 
    # using the following command: 
    tar -xhzvf smrtanalysis.tar.gz -C /home/hadoop
    ```
    **Note**:
    
    (1). The decompressed location of smrtanalysis.tar.gz must be the same on all worker nodes. Don't forget to set the variable *SMRT\_ANALYSIS\_HOME* in parameters.conf. (Suppose you have decompressed _smrtanalysis.tar.gz_ to _/home/hadoop_ on all worker nodes, then you have to set *SMRT\_ANALYSIS\_HOME=/home/hadoop/smrtanalysis* in parameters.conf)
    
    (2). To preserve symbolic links in the _tar.gz_ file, "_-h_" must be used when using _tar_ command.


3. #### Python 2.x and required Python libraries

    If the OSs of nodes (both master and workers) in your cluster don't have python 2.x installed, you should install it. Install  *numpy*, *h5py*, *paramiko*, *pbcore* in your python environment. Install package *py4j*, *pyspark* in your python environment if you need to.
    
    Note that **Python 2.7.13** (or higher) is strongly recommended (not necessary) because of the bug described in [issue #5](https://github.com/PengNi/basemods_spark/issues/5). 


### How to use basemods_spark

1. #### make the scripts executable

    If the scripts in the code of basemods_spark you download don't have execute permissions, you should make them executable.
    
    ```sh
    chmod +x basemods_spark/scripts/exec_sawriter.sh
    
    chmod +x basemods_spark/scripts/baxh5_operations.sh
    
    chmod +x basemods_spark/scripts/cmph5_operations.sh
    
    chmod +x basemods_spark/scripts/mods_operations.sh
    ```

2. #### copy your data
    Copy your data to the master node of your *Hadoop/Spark* cluster.


3. #### parameters in configure file
    Set the parameters in configure file 'parameters.conf'.


4. #### start Spark and use spark-submit to run the pipeline

    (1) start HDFS if you need to

    ```sh
    $HADOOP_HOME/sbin/start-dfs.sh
    ```

    (2) start Spark
    ```sh
    $SPARK_HOME/sbin/start-all.sh
    ```

    (3) submit your job
    ```sh
    $SPARK_HOME/bin/spark-submit basemods_spark_runner.py
    ```
