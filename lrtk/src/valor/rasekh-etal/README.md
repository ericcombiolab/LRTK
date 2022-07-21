********************************************************************************************

                VALOR : Discover inversions using long range information

********************************************************************************************

This algorithm is designed for discovering large inversions using long range information such as pooled clone sequencing and 10X linked reads. The input is a set of BAM files of mapped paired end reads. 

NOTE: In the BAM or bed files, it is expected to see the chromosome label as: chr# for example chr1, chr2, ..., chrX, chrY. Make sure the reference is the same or change the run.sh code.

NOTE: BAM files should be sorted by read name. (samtools sort -n) for bamtobed to work correctly.

NOTE: The algorithm does not rely on the reference. So hg18, hg19 or other are accepted. ONLY the chromosome name should be as : chr#



First you should set the variables. 

********************************************************************************************

                                edit src/Config.java

********************************************************************************************

Set the variables in the src/Config.java file. These are needed by the Java executables.

* You should set the READ_LENGTH, FRAG_SIZE, CLONE_MEAN, and CLONE_STD_DEV from the data.

* INV_MAX_SIZE and INV_MIN_SIZE are user specific.

* The rest will be set automatically.



-- Set the read information

parameter: READ_LENGTH (length of each read)

parameter: FRAG_SIZE (max length of a normal segment)



-- Set the physical statistics of clones. Changing these variables adaptively can improve the performance.

parameter: CLONE_MEAN (expected average of clone size: 150K for BAC and 40K for FOSMID)

parameter: CLONE_STD_DEV (expected standard deviation of the clones)



-- Set the inversion information. It is suggested to run the algorithm on narrowed ranges of inversion length.

parameter: INV_MIN_SIZE (minimum size of an inversion, should be at least 2*NORMAL_SIZE)

parameter: INV_MAX_SIZE = (maximum size of an inversion)




********************************************************************************************

                                   edit run.sh

********************************************************************************************

Edit the run.sh file:



-- The path to bam files

parameter: bamfile (the directory path of the pooled bam files containing the mapped reads)

parameter: count (the maximum number of pools)

parameter: suffix (the suffix of files)

Take notice that the files will be read as $bamfile1$suffix, $bamfile2$suffix, ..., $bamfile$count$suffix

NOTE: bamfiles should be sorted by read name. (samtools sort -n) for bamtobed to work correctly.



-- The paired end distance

parameter: FRAG_MAX (maximum length of a normal segment)

parameter: FRAG_MIN (the minimum length of a normal segment)

You should set the max and min to (mean + 3*std) and (mean - 3*std).



-- The output directory

parameter: outputDir (output/working directory)

The program will output the results here. It will make a temp folder for temporary files it needs which will be deleted at the end.

The inferred clones will be output as allRegions.tsv here along with the +/+ and -/- mapping reads.

Also it makes directories Clusters and Pools which the final inversion clusters and the split clones will be output.



-- The threshold to filter out the inferred clones by coverage

parameter: coverageThreshold (the coverage threshold for inferred clones by reads)

This step is optional. It is shown that many inferred clones are not reliable are are due to extension of discrete reads on duplicated areas. The default is 40%. But it is better to get the coverage for all inferred clones and plot the correlation between the size and coverage and come up with a meaningful number for your case.



********************************************************************************************

                                  run run.sh

********************************************************************************************

Then run the run.sh file:

sh run.sh

The output clusters will be in the output directory.

1. Before running it is advised to remove duplicated reads. This will save much time. Use rmdup but be aware that there are some bugs in the older versions.

2. After getting the inferred clones, it is advised to filter out the clones with lower than average coverage. Use coverageBed. However if you wish to skip this part you can comment it out.

3. After the execution and getting the inversions, check for overlaps with known gaps and remove the ones overlapping too much. Use bedtools intersect -r -f 0.01



NOTE: In the bam file or bed files, it is expected to see the chromosome label as: chr# for example chr1, chr2, ..., chrX, chrY. Make sure the reference is the same or change the run.sh code.





