# FastDPeak
Density Peak is not applicable for large scale data, due to the two quantities, i.e, density ρ and δ, are both obtained by brute force algorithm with complexity O(n2). Then, a simple but fast Density Peak Clustering, namely FastDPeak, is proposed, which runs in about O(nlog(n)) expected time in the intrinsic dimensionality. It replaces density with kNN-density, which is computed by fast kNN algorithm such as cover tree, yielding the huge improvement for density computations. Based on kNN-density, local density peaks and non-local density peaks are identiﬁed, and a fast algorithm, whichusestwodifferentstrategiestocomputeδ for them, is also proposed with complexity O(n). Experimental results show that FastDPeak is effective and outperforms other variants of DPeak.

***********************************************************************************
 The FastDPeak program was compiled under Windows using c++ with CodeBlocks 10.05.
*********************************************************************************** 

***********************************************************************************
 Files
***********************************************************************************
Unzip the "FastDPeak.zip" file, which will create folder mainly containing:
- a project file named "FastDPeak.cbp".
- a dataset folder named "data".
- a c++ file named “main.cpp” is the main function file.
- other c++ files
***********************************************************************************
 Environment configuration
***********************************************************************************

Step1:
- Download CodeBlocks in http://www.codeblocks.org/
- Download TDM-GCC-32 in http://tdm-gcc.tdragon.net/download/   

Step2:
- Open CodeBlock: choose “setting” ->”compiler and debugger”->”ToolChain executables”,and set the parameters like the “Figure 1” 

Step3:
***********************************************************************************
 Dataset Formation
***********************************************************************************
The dataset should be given in a text file with the following formation:

- each line represents a point with d numbers, where d is dimension:

For instance, the first 20 lines of the sample dataset "agg.txt" are shown as below:
1.2    1.6    
3      5      
2      4.6    
10     2      
2.1    4.1    
3.5    5.1    
6.6    1      
3.6    4.2    
3      10.3   
6.8    7.2    
1.1    1.5    
3.3    2.5    
8.6    9.2    
1.8    2.2    
11.2   8.4    
7.9    8.8    
9.2    9.8    
9.4    9      
10.2   2      
2.4    3.6    

In this example, there are 2 numbers in each line: the first line represents the first point with the coordinates (1.2   1.6), and its id is 1. Similarly, the rest nine lines specify the coordinates of the point with id = 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,
16,17,18,19,20 respectively.  

***********************************************************************************
 Data Download
***********************************************************************************

The data used in the paper can be downloaded from https://pan.baidu.com/s/1Yyawj0iwM2Mu8zrYioWb-w  Password: 3yje, as Figure 2 shows.

-“new_KDD_data.txt” is “KDD99”.
-“new_bio_train.txt” is “KDD04”.
-“new_BigCross500K.txt” is “BigCross”.
-“synthesis_1.txt” is “SYN1”.
-“synthesis_2.txt” is “SYN2”.
-“synthesis_3.txt” is “SYN3”.

***********************************************************************************
 An example of quick start
***********************************************************************************
Step1:
Open project “FastDPeak.cbp” in Codeblocks.

Step2:
Open “main.cpp”,as figure 4 shows.

In line 720 : float *raw_data = read_data(data_file_name," ",&dim, &data_size);
-The first parameter data_file_name represents the origin dataset.
-The second parameter dim represents dim of the origin dataset.
-The third parameter data_size represents size of the origin dataset.

In line 722: 
float *data_with_index = Generate_data_with_index(raw_data,data_size,dim,new_size);
-The first parameter raw_data represents the origin dataset.
-The second parameter data_size represents size of the origin dataset.
-The third parameter dim represents dim of the origin dataset.
-The fourth parameter new_size represents size of new dataset which we would like to classify.

In line 730-731: 
Fast_Density_Peak(K,data_with_index,new_size,batch_num,dim,local_peak_threshold
,log_file,results,cl,cl_results,dis_matrix,node_p_ptr);
-The first parameter K represents the K of KNN.
-The second parameter data_with_index represents new dataset which we would like to classify.data_size represents size of the origin dataset.
-The third parameter new_size represents size of new dataset which we would like to classify.
-The fourth parameter batch_num represents the number of batch in creating Covertree.
-The fifth parameter dim represents dim of the origin dataset.
-The sixth parameter local_peak_threshold represents threshold of local peak points.
-The sixth parameter log_file represents log file which record records of experiment.
-The seventh parameter results represents the tree of new dataset by FastDPeak.
-The eighth parameter cl represents the number of clusters.
-The ninth parameter cl_results represents the designate file in which FastDPeak will write the clustering result.
-The tenth parameter dis_matrix represents the matrix of distance.
-The eleventh parameter node_p_ptr represents the density and index of every point.

Step3:
-Press the “Build and run” button in CodeBlocks under release mode.

***********************************************************************************
 Output Format
***********************************************************************************

The clustering result is saved in “data\cl_results.txt”, and the output formation is:

1.2    1.6    2
3      5      1
2      4.6    1
10     2      1
2.1    4.1    1
3.5    5.1    1
6.6    1      1
3.6    4.2    1
3      10.3   1
6.8    7.2    2
1.1    1.5    1
3.3    2.5    2
8.6    9.2    1
1.8    2.2    1
11.2   8.4    2
7.9    8.8    2
9.2    9.8    1
9.4    9      1
10.2   2      2
2.4    3.6    2

The first two columns store coordinates of all point, and the last column represents the cluster ID.

For example, the first line (1.2  1.6  2) means that the first point is classified into cluster 2. If cluster ID is -1, it means this point is a noise.

***********************************************************************************
 Experiments
***********************************************************************************
Some experiments we did in paper are listed in experiment_records.doc

