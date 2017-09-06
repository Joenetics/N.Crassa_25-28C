# N.Crassa_25-28C
Masters Project examining the role of VVD protein in a wt and VVD-KO across two temperatures (25 &amp; 28C).
There are many files on here.

1)Some excel sheets for simple counting of genes involved in some basic comparisons
2)Both RStudio code I was given (OriginalRCode.r) and what i turned it into (AllTheRCode.r)
3)Some Python files;
  a)UpstreamFinder.py = finds upstream regions of genes based on upstream elements file (NC_Upstream.fasta))and macthing them to query accession codes
  b)UpstreamUniques.py = same principle, but only counts genes that ARE differentially expressed in VVD but NOT in the WT
  c)DUmmyMaker.py = Little program to generate random 1kb sequences with defined motifs randomly within; this is to test motif software.
  d)HTSeqMerger.py = little program to join HTSeq files into format suitable for DESeq2 package on RStudio. 
4)a)htseq_count = folder with all the htseq count files.
  b)MergedHTSeq.txt = All the HTSeq merged for different RStudio package; done with Python program I created 'HTSEqMerger.py'
