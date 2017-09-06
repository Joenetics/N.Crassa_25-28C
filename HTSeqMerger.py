#Open a number of HTSeqFiles, then merge their counts per gene; this is for RStudio packages that require merged CSVs.
#This is for 40 wells. I'm not sure how to create variable names in a loop in Python... plus, no regularity in naming files. so, i did it manually.


HTSeqFile1 = open("/Users/langyy/Desktop/Differential_expression/htseq_count/MV01_rep01_htseq_count",'r')
HTSeqFile2 = open("/Users/langyy/Desktop/Differential_expression/htseq_count/MV01_rep02_htseq_count",'r')
HTSeqFile3 = open("/Users/langyy/Desktop/Differential_expression/htseq_count/MV02_rep01_htseq_count",'r')
HTSeqFile4 = open("/Users/langyy/Desktop/Differential_expression/htseq_count/MV02_rep02_htseq_count",'r')
HTSeqFile5 = open("/Users/langyy/Desktop/Differential_expression/htseq_count/MV03_rep01_htseq_count",'r')
HTSeqFile6 = open("/Users/langyy/Desktop/Differential_expression/htseq_count/MV03_rep02_htseq_count",'r')
HTSeqFile7 = open("/Users/langyy/Desktop/Differential_expression/htseq_count/MV04_rep01_htseq_count",'r')
HTSeqFile8 = open("/Users/langyy/Desktop/Differential_expression/htseq_count/MV04_rep02_htseq_count",'r')
HTSeqFile9 = open("/Users/langyy/Desktop/Differential_expression/htseq_count/MV05_rep01_htseq_count",'r')
HTSeqFile10 = open("/Users/langyy/Desktop/Differential_expression/htseq_count/MV05_rep02_htseq_count",'r')
HTSeqFile11 = open("/Users/langyy/Desktop/Differential_expression/htseq_count/MV06_rep01_htseq_count",'r')
HTSeqFile12 = open("/Users/langyy/Desktop/Differential_expression/htseq_count/MV06_rep02_htseq_count",'r')
HTSeqFile13 = open("/Users/langyy/Desktop/Differential_expression/htseq_count/MV07_rep01_htseq_count",'r')
HTSeqFile14 = open("/Users/langyy/Desktop/Differential_expression/htseq_count/MV07_rep02_htseq_count",'r')
HTSeqFile15 = open("/Users/langyy/Desktop/Differential_expression/htseq_count/MV08_rep01_htseq_count",'r')
HTSeqFile16 = open("/Users/langyy/Desktop/Differential_expression/htseq_count/MV08_rep02_htseq_count",'r')
HTSeqFile17 = open("/Users/langyy/Desktop/Differential_expression/htseq_count/MV09_rep01_htseq_count",'r')
HTSeqFile18 = open("/Users/langyy/Desktop/Differential_expression/htseq_count/MV09_rep02_htseq_count",'r')
HTSeqFile19 = open("/Users/langyy/Desktop/Differential_expression/htseq_count/MV10_rep01_htseq_count",'r')
HTSeqFile20 = open("/Users/langyy/Desktop/Differential_expression/htseq_count/MV10_rep02_htseq_count",'r')
HTSeqFile21 = open("/Users/langyy/Desktop/Differential_expression/htseq_count/MV11_rep01_htseq_count",'r')
HTSeqFile22 = open("/Users/langyy/Desktop/Differential_expression/htseq_count/MV11_rep02_htseq_count",'r')
HTSeqFile23 = open("/Users/langyy/Desktop/Differential_expression/htseq_count/MV12_rep01_htseq_count",'r')
HTSeqFile24 = open("/Users/langyy/Desktop/Differential_expression/htseq_count/MV12_rep02_htseq_count",'r')
HTSeqFile25 = open("/Users/langyy/Desktop/Differential_expression/htseq_count/MV13_rep01_htseq_count",'r')
HTSeqFile26 = open("/Users/langyy/Desktop/Differential_expression/htseq_count/MV13_rep02_htseq_count",'r')
HTSeqFile27 = open("/Users/langyy/Desktop/Differential_expression/htseq_count/MV14_rep01_htseq_count",'r')
HTSeqFile28 = open("/Users/langyy/Desktop/Differential_expression/htseq_count/MV14_rep02_htseq_count",'r')
HTSeqFile29 = open("/Users/langyy/Desktop/Differential_expression/htseq_count/MV15_rep01_htseq_count",'r')
HTSeqFile30 = open("/Users/langyy/Desktop/Differential_expression/htseq_count/MV15_rep02_htseq_count",'r')
HTSeqFile31 = open("/Users/langyy/Desktop/Differential_expression/htseq_count/MV16_rep01_htseq_count",'r')
HTSeqFile32 = open("/Users/langyy/Desktop/Differential_expression/htseq_count/MV16_rep02_htseq_count",'r')
HTSeqFile33 = open("/Users/langyy/Desktop/Differential_expression/htseq_count/MV17_rep01_htseq_count",'r')
HTSeqFile34 = open("/Users/langyy/Desktop/Differential_expression/htseq_count/MV17_rep02_htseq_count",'r')
HTSeqFile35 = open("/Users/langyy/Desktop/Differential_expression/htseq_count/MV18_rep01_htseq_count",'r')
HTSeqFile36 = open("/Users/langyy/Desktop/Differential_expression/htseq_count/MV18_rep02_htseq_count",'r')
HTSeqFile37 = open("/Users/langyy/Desktop/Differential_expression/htseq_count/MV19_rep01_htseq_count",'r')
HTSeqFile38 = open("/Users/langyy/Desktop/Differential_expression/htseq_count/MV19_rep02_htseq_count",'r')
HTSeqFile39 = open("/Users/langyy/Desktop/Differential_expression/htseq_count/MV20_rep01_htseq_count",'r')
HTSeqFile40 = open("/Users/langyy/Desktop/Differential_expression/htseq_count/MV20_rep02_htseq_count",'r')

outputfile = open("/Users/Neurospora/Desktop/Joseph_Things/MergedHTSeq.txt", 'w')
outputfile2 = open("/Users/Neurospora/Desktop/Joseph_Things/GeneNames.txt", 'w')

MIGHTYDICTIONARY = {}

MightyList = [HTSeqFile1, HTSeqFile2, HTSeqFile3, HTSeqFile4, HTSeqFile5, HTSeqFile6, HTSeqFile7, HTSeqFile8, HTSeqFile9, HTSeqFile10,
              HTSeqFile11, HTSeqFile12, HTSeqFile13, HTSeqFile14, HTSeqFile15, HTSeqFile16, HTSeqFile17, HTSeqFile18, HTSeqFile19, HTSeqFile20,
              HTSeqFile21, HTSeqFile22, HTSeqFile23, HTSeqFile24, HTSeqFile25, HTSeqFile26, HTSeqFile27, HTSeqFile28, HTSeqFile29, HTSeqFile30,
              HTSeqFile31, HTSeqFile32, HTSeqFile33, HTSeqFile34, HTSeqFile35, HTSeqFile36, HTSeqFile37, HTSeqFile38, HTSeqFile39, HTSeqFile40]

GeneNames = []
for line in HTSeqFile1:
    word = line.split("\t")
    MIGHTYDICTIONARY[word[0]] = (word[1]).strip()
    GeneNames.append(word[0])

previousnumbers = str
for i in range(1,40):
    for line in MightyList[i]:
        word = line.split("\t")
        previousnumbers = MIGHTYDICTIONARY[word[0]] + "\t"
        MIGHTYDICTIONARY[word[0]] = previousnumbers +  (word[1]).strip()
print("VVD25DD0R1", "VVD25DD0R2", "VVD25DD4R1", "VVD25DD4R2", "VVD25DD8R1", "VVD25DD8R2", "VVD25DD12R1", "VVD25DD12R2", "VVD25DD16R1", "VVD25DD16R2",
                                    "VVD28DD0R1", "VVD28DD0R2", "VVD28DD4R1", "VVD28DD4R2", "VVD28DD8R1", "VVD28DD8R2", "VVD28DD12R1", "VVD28DD12R2", "VVD28DD16R1", "VVD28DD16R2",
                                    "WT25DD0R1", "WT25DD0R2", "WT25DD4R1", "WT25DD4R2", "WT25DD8R1", "WT25DD8R2", "WT25DD12R1", "WT25DD12R2", "WT25DD16R1", "WT25DD16R2",
                                    "WT28DD0R1", "WT28DD0R2", "WT28DD4R1", "WT28DD4R2", "WT28DD8R1", "WT28DD8R2", "WT28DD12R1", "WT28DD12R2", "WT28DD16R1", "WT28DD16R2", sep = "\t", file = outputfile)
for i in range(0, len(GeneNames)):
    print(GeneNames[i], file= outputfile2)
    print(MIGHTYDICTIONARY[GeneNames[i]], "\n", sep = '', file= outputfile)

