## basemods on spark

### Overview
the spark version of basemods pipeline in SMRT-Analysis

### Set up the environment
1. #### Hadoop/Spark

    Setting up an [Hadoop](http://hadoop.apache.org/)/[Spark](https://spark.apache.org/) cluster.

2. #### SMRT-Analysis

    2.1 Install [SMRT-Analysis](https://github.com/PacificBiosciences/SMRT-Analysis) (v2.3.0) in **each worker node** of your *Hadoop/Spark* cluster. (or use NFS?)
    
    2.2 In **each worker node**, replace *ipdsummary.py* and *KineticWorker.py* in SMRT-Analysis with the modified scripts in the directory "basemods_spark/scripts".
     + copy *basemods_spark/scripts/ipdSummary.py* to *$SMRT_HOME/install/smrtanalysis_2.3.0.140936/analysis/bin*

      ```sh
      chmod +x basemods_spark/scripts/ipdSummary.py
      
      sudo mv $SMRT_HOME/install/smrtanalysis_2.3.0.140936/analysis/bin/ipdSummary.py $SMRT_HOME/install/smrtanalysis_2.3.0.140936/analysis/bin/ipdSummary.py.bak
      
      sudo cp basemods_spark/scripts/ipdSummary.py $SMRT_HOME/install/smrtanalysis_2.3.0.140936/analysis/bin/
      ```

    + copy *basemods_spark/scripts/KineticWorker.py* to *$SMRT_HOME/install/smrtanalysis_2.3.0.140936/analysis/lib/python2.7/kineticsTools*
   
    ```sh
    sudo mv $SMRT_HOME/install/smrtanalysis_2.3.0.140936/analysis/lib/python2.7/kineticsTools/KineticWorker.py $SMRT_HOME/install/smrtanalysis_2.3.0.140936/analysis/lib/python2.7/kineticsTools/KineticWorker.py.bak
    
    sudo cp basemods_spark/scripts/KineticWorker.py $SMRT_HOME/install/smrtanalysis_2.3.0.140936/analysis/lib/python2.7/kineticsTools/
    ```

3. #### Python 2.x and required Python libraries

    If the OSs of nodes (both master and workers) in your cluster don't have python 2.x installed, you should install it. Install  *numpy*, *h5py*, *paramiko*, *pbcore* in your python environment. Install package *py4j*, *pyspark* in your python environment if you need to.
    
    Note that **Python 2.7.13** (or higher) is strongly recommended (not necessary) because the bug described in [issue #5](https://github.com/PengNi/basemods_spark/issues/5). 


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
