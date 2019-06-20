# nanoKit_py3
this kit design to convert nanopore sequnecing data fast5 to fastq read id/read name convertion. Also, I built a N50 calculator designed for assembly.
##fast5toreadid.py
This tool is designed to extract million of fast5 files from nanopore sequencing to basecalled read id (read name)
'''
-i --input directory
-s --save directory
'''
Test tool using:

'''
fast5toreadid.py -i fast5_test/ -s fast5toreadid_test_result/
'''

##single_fast5toreadid.py
'''
-i --input single fast5
'''
This tool is designed to extract single fast5 file to basecalled read id (read name)

'''
fast5toreadid.py -i fast5_test/0/GXB01298_20190601_FAK27254_GA10000_sequencing_run_juFosSP1_24_85131_read_4_ch_7_strand.fast5
'''
##Author
* **Qiutao DING** - *Initial work* 
