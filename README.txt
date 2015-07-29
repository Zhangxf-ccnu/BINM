README file for Matlab source code supporting the paper "Identifying binary protein-protein interactions from affinity purication mass spectrometry data".


Contents of this archive
------------------------
This archive contains several Matlab scripts used to identify binary protein-protein interactions using the algorithm BINM presented in the above paper. 

(1) BINM.m: Matlab script for the core of BINM which learns parameter W_{dir} from a given AP-MS PPI network with weighted adjacent matrix W_{obs}.
 
(2) BINM_main.m: Matlab script for the main function of BINM which reads data from text files first and then identifies binary protein interactions using BINM.m. It also writes the confidence scores of direct interactions into a file "output_file_name" provided by the user or the default file "result.txt"

(3) BINM_demo.m: A simple Matlab script to test BINM. When a data set is chose, it can be run in a straightforward manner within a Matlab window.

This archive also contains a folder named as "data" which includes the four AP-MS datasets used in this study.

Please do not hesitate to contact Prof. Dao-Qing Dai at stsddq@mail.sysu.edu.cn (or Dr. Xiao-Fei Zhang at zhangxf@mail.ccnu.edu.cn) to seek any clarifications regarding any contents or operation of the archive.
