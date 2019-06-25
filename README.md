# nanoKit_py3  
This kit design to convert **nanopore sequnecing data** fast5 to fastq read id/read name convertion. Also, I built a N50 calculator designed for assembly.   
:rocket:   
----  
## fast5toreadid.py   
  This tool is designed to extract million of fast5 files from nanopore sequencing to basecalled read id (read name)  
-i --input directory  
-s --save directory  
### Quick start  
    fast5toreadid.py -i fast5_test/ -s fast5toreadid_test_result/  
  
## single_fast5toreadid.py  
-i --input single fast5  
### Quick start  
    single_fast5toreadid.py -i fast5_test/0/GXB01298_20190601_FAK27254_GA10000_sequencing_run_juFosSP1_24_85131_read_4_ch_7_strand.fast5  
## FosmidSequencingRemoveVector.py  
-i --input fasta reads to trimming and filtering
-r --reference reference sequence to be trimmed from nanopore reads
-o --output output directory to store the trimmed sequence
### Quick start  
    FosmidSequencingRemoveVector.py -r fosmid_sample/pSMART_BAC.fa -i fosmid_sample/sample03.fa -o fosmid_sample/result/trimmed_reads.fa > fosmid_sample/result/not_vector.fa  
## Author  
* **Qiutao DING** - *Initial work* 
