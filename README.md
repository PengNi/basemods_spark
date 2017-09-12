## basemods on spark

### Overview
the spark version of basemods pipeline in SMRT-Analysis

### Set up the environment
1. #### Hadoop/Spark

    Setting up an [Hadoop](http://hadoop.apache.org/)/[Spark](https://spark.apache.org/) cluster.

2. #### SMRT-Analysis

    2.1 Install [SMRT-Analysis](https://github.com/PacificBiosciences/SMRT-Analysis) in each worker node of your *Hadoop*/*Spark* cluster.
    
    2.2 In each worker node, replace *ipdsummary.py* and *KineticWorker.py* in SMRT-Analysis with the modified scripts in our code.
     + copy *basemods_spark/scripts/ipdSummary.py* to *$SMRT_HOME/install/smrtanalysis_2.3.0.140936/analysis/bin*
       
      ```sh
      chmod +x basemods_spark/scripts/ipdSummary.py
      
      sudo mv $SMRT_HOME/install/smrtanalysis_2.3.0.140936/analysis/bin/ipdSummary.py $SMRT_HOME/install/smrtanalysis_2.3.0.140936/analysis/bin/ipdSummary.py.bak
      
      sudo cp basemods_spark/scripts/ipdSummary.py $SMRT_HOME/install/smrtanalysis_2.3.0.140936/analysis/bin
      ```
    + copy *basemods_spark/scripts/KineticWorker.py* to *$SMRT_HOME/install/smrtanalysis_2.3.0.140936/analysis/lib/python2.7/kineticsTools*
    
    ```sh
    sudo mv $SMRT_HOME/install/smrtanalysis_2.3.0.140936/analysis/lib/python2.7/kineticsTools/KineticWorker.py $SMRT_HOME/install/smrtanalysis_2.3.0.140936/analysis/lib/python2.7/kineticsTools/KineticWorker.py.bak
    
    sudo cp basemods_spark/scripts/KineticWorker.py $SMRT_HOME/install/smrtanalysis_2.3.0.140936/analysis/lib/python2.7/kineticsTools/
    ```

3. #### Python 2.x and required Python libraries

    If your OS doesn't have python 2.x installed, you should install it. Install package *py4j*, *pyspark* to your python environment if you need to.


### How to use basemods_spark

1. #### make scripts in our code executable

If the scripts in our code don't have execute permissions, you should make them executable
    
    ```sh
    chmod +x basemods_spark/scripts/baxh5_operations.sh
    
    chmod +x basemods_spark/scripts/cmph5_operations.sh
    
    chmod +x basemods_spark/scripts/mods_operations.sh
    ```

2. #### use spark-submit

