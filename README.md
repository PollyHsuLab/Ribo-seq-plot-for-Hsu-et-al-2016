### Ribo-seq 3 nt periodicity plot 
This code was used for plotting the 3-nt periodicity from Ribo-seq data together with the RNA-seq result.  

#### What is required:
1. A gtf file (make sure the 1st column (chromosome name) is numeric, e.g. 10, not Chr10)  
2. A bam file for the corresponding RNA-seq reads  
3. A txt file generated from RiboTaper output p_sites_all file that containing p-site counts, p-site position, chromosome info, strand info.   

#### Note:
1. This is designed for Arabidopsis gtf file so it automatically extract the chromosome number from the 3rd character of the gene name. It has to be modified for other organism. Change the 'chr <- as.numeric(substr(YFG,3,3))' in the 'PLOT' function accordingly.

**Citation**: Super-resolution ribosome profiling reveals unannotated translation events in Arabidopsis. Proc Natl Acad Sci USA. 2016 Nov 8;113(45):E7126-E7135. doi: 10.1073/pnas.1614788113. Epub 2016 Oct 21
