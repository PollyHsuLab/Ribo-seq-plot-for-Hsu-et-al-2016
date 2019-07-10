### Ribo-seq 3 nt periodicity plot 
This code was used for plotting the 3-nt periodicity from Ribo-seq data together with the RNA-seq result.  

For example: A regular gene
![alt text](https://github.com/PollyHsuLab/Ribo-seq-plot-for-Hsu-et-al-2016/blob/master/Rplot01.jpg)

For example: A gene with uORF
![alt text](https://github.com/PollyHsuLab/Ribo-seq-plot-for-Hsu-et-al-2016/blob/master/Rplot02.jpg)

#### What is required:
1. A gtf file (make sure the 1st column (chromosome name) is numeric, e.g. 10, not Chr10)  
2. A bam file for the corresponding RNA-seq reads  
3. A txt file generated from RiboTaper output P_sites_all file that containing p-site counts, p-site position, chromosome info, strand info.   

#### Note:
1. This is designed for Arabidopsis gtf file so it automatically extract the chromosome number from the 3rd character of the gene name. It has to be modified for other organism. Change the `chr <- as.numeric(substr(YFG,3,3))` in the 'PLOT' function accordingly.  
2. you can plot uORFs For example: `PLOT(YFG ="AT2G18160",uORF="AT2G18162",CDSonly=F)`    
   However, if you find new uORFs and want to plot it, you will need to have the coordinate and uORF names in the gtf.  
3. gtf only, no gff  
4. You can plot different isoform: `PLOT(YFG ="AT2G46980",isoform=2)`  

**Citation**: Super-resolution ribosome profiling reveals unannotated translation events in Arabidopsis. Proc Natl Acad Sci USA. 2016 Nov 8;113(45):E7126-E7135. doi: 10.1073/pnas.1614788113. Epub 2016 Oct 21

### RiboPlotR release
Github: https://github.com/hsinyenwu/RiboPlotR
Paper: https://www.biorxiv.org/content/10.1101/694646v1.article-metrics
