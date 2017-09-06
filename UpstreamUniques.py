#This is redundant with UpstreamFinder.py and other RStudio scripts!
print('This is redundant with UpstreamFinder.py and other RStudio scripts!')

#same as UpstreamFinder.py, but only counts gene for vvd if up/downreg isnt seen in wt.
#Thereofre, select for uniques.


import sys

try:
    DEfile = open("/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_vvd25_0_4.txt",'r') #vvd here.
    DEfile2 = open("/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_wt25_0_4.txt",'r') #wt here.
    upstreamfile = open("/Users/Neurospora/Desktop/Joseph_Things/NC_Upstream.fasta",'r')
    outputfile = open("UpstreamFound_vvd25_0_4_corrected.txt", "w")                                         #Added 'corrected' to name of file
    namesoutput = open("UpstreamFound_vvd25_0_4_names.txt", "w")

except:
    print("Files not found: Probably wrong directory.")
    sys.exit()

validletters = "UBDubd"
upordown2 = input("Up or Down reg or both? (U/D/B) : ")
if upordown2 in validletters:
    print("Assuming ", upordown2.upper(), sep='')
else:
    print("Please only use single-letter option U/B/D.")
    sys.exit()
upordown = upordown2.upper()
numberlist = []
numberlist1 = []
numberlist2 = []
WholeFile1 = {}
WholeFile2 = {}
numba = "test"
sequence = "also test"
sequencedict = {}
number = str


for line in DEfile:
    if line.startswith("\"NCU"):
        data = line.split(" ")
        number = (data[0])[1:-1]
        numberlist1.append(number)
        WholeFile1[(data[0])[1:-1]] = data[2]
for line in DEfile2:
    if line.startswith("\"NCU"):
        data = line.split(" ")
        number = (data[0])[1:-1]
        numberlist2.append(number)
        WholeFile2[(data[0])[1:-1]] = data[2]

duplicatecounter1 = 0
duplicatecounter2 = 0
for u in range(0, len(numberlist1)):
    if numberlist1[u] in numberlist2:
        duplicatecounter1 += 1

        if (float(WholeFile1[numberlist1[u]])- float(WholeFile2[numberlist1[u]])) >= 1 and (upordown == 'U' or upordown == 'B'):
            print(numberlist1[u], " has been upregulated.", sep = '')
        elif (float(WholeFile1[numberlist1[u]])- float(WholeFile2[numberlist1[u]])) <= -1 and (upordown == 'D' or upordown == "B"):
            print(numberlist1[u], " has been downregulated.", sep='')
        numberlist.append(numberlist1[u])
    elif float(WholeFile1[numberlist1[u]]) >= 1 and (upordown == 'U' or upordown == 'B'):
        numberlist.append(numberlist1[u])
        duplicatecounter2 += 1
    elif float(WholeFile1[numberlist1[u]]) <= -1 and (upordown == 'D' or upordown == 'B'):
        numberlist.append(numberlist1[u])
        duplicatecounter2 += 1


##make dictionary for NCU numnbers and upstream elements.
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
        print(numberlist1[i], "\n", sep = '', file = namesoutput)
    else:
        counterforfailmatch += 1
        failedmatchnumbers.append(numberlist[i])
print("Reduced from ", len(numberlist1), " to ", len(numberlist))
print("Number of upstream elements found successfully: ", len(numberlist) - counterforfailmatch)
print("Number of failed matches: ", counterforfailmatch, "\nAccession numbers of failed matches: ", failedmatchnumbers)
