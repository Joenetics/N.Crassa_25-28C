import sys

try:
    #DEfile = open("/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_vvd_wt_DD0_25.txt",'r')       #This one goes to classic DESeq files.
    DEfile = open("/Users/Neurospora/Desktop/Joseph_Things/31shared.txt", 'r')
    upstreamfile = open("/Users/Neurospora/Desktop/Joseph_Things/NC_Upstream.fasta",'r')    #File with all the upstream 1kb regions.
    outputfile = open("31shared.fasta", "w")     #open("WT_VVD_DD0_25C_Upsteam.fasta", "w")   #different files to print to with differnt conditions ('upreg in VVD' / 'only if increase by 2fold'/...)
    outputfile2 = open("WT_VVD_DD0_25C_UpsteamVVDupreg.fasta", "w")
    outputfile3 = open("WT_VVD_DD0_25C_UpsteamWTupreg.fasta", "w")
except:
    print("Files not found: Probably wrong directory.")
    sys.exit()

numberlist = []
numberlist1 = []
numberlist2 = []
numba = "test"
sequence = "also test"
sequencedict = {}
number = str

###this is for only custom lists of gene accession numbers, one per line and nothing else. (e.g, "1"\n"2"\n...)
#for line in DEfile:
#    numberlist.append(line.strip())


for line in DEfile: #This stuff works for files from DESeq2 package, from RStudio. Otherwise, extract gene names differently
    if line.startswith("\"NCU"):
        data = line.split(" ")
        if float(data[2]) <= (999.0): # pretty much always print. This will use an output file with ALL genes names
            number = (data[0])[1:-1]
            numberlist.append(number)
        if float(data[2]) <= (0.0): # if expression change is less than 0 => upreg. => add to 'upreg. in WT' file
            number = (data[0])[1:-1]
            numberlist1.append(number)
        if float(data[2]) >= (0.0): # opposite to previous.
            number = (data[0])[1:-1]
            numberlist2.append(number)


##make dictionary for NCU numnbers and upstream elements (1KB). this will then be used later with list from before.
for line in upstreamfile:
    if line.startswith(">"):
        sequencedict[numba] = ''.join(sequence)
        sequence = []
        data = line.split(" ")
        numba = (data[1])[1:-1]
    else:
        sequence.append((line).strip())

counterforfailmatch = 0
failedmatchnumbers = []
for i in range(0, len(numberlist)):
    if numberlist[i] in sequencedict:
        print(">", str(numberlist[i]), "\n", sequencedict[str(numberlist[i])], sep = '', file = outputfile)
    else:
        counterforfailmatch += 1
        failedmatchnumbers.append(numberlist[i])
print("Number of upstream elements found successfully: ", len(numberlist) - counterforfailmatch)
print("Number of failed matches: ", counterforfailmatch, "\nAccession numbers of failed matches: ", failedmatchnumbers)


#these will only print to file.
for i in range(0, len(numberlist1)):
    if numberlist[i] in sequencedict:
        print(">", str(numberlist1[i]), "\n", sequencedict[str(numberlist1[i])], sep = '', file = outputfile2)

for i in range(0, len(numberlist2)):
    if numberlist2[i] in sequencedict:
        print(">", str(numberlist2[i]), "\n", sequencedict[str(numberlist2[i])], sep = '', file = outputfile3)
