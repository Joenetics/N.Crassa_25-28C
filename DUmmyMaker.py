#Joseph Shepherd 2017 RP2
#This program creates 20 random DNA sequences of 1kb each. Each sequence has, at a random location, a selected motif thrown in.
#This motif is therefore present atleast 20 times; once per sequence. Then, (if it is longer than 4nt), serially shortened versions are
#tested; their occurrance should rise as they become less unique. You may also test a motif's presence in the randomly generated sequences.
#I.e, see if a test motif appears at all in the sequences.

import random
from random import randint
import sys
import re
import matplotlib.pyplot as plt

ValidBases = "ATGCatgc"
ValidBases2 = "ATGCatcg[]"
motifinput = input("Please input DNA string of more than 4nt for motif insertion:")
validnumbers = "0123456789."

if len(motifinput.strip("\[\]")) <4: # if motif too long/ short, don't test it
    print("String is incorrect length. Dna of more than 4 nt only.")
    sys.exit()

for i in range(0, len(motifinput)): #If non-DNA char, dont create program. Waste of output file time.
    if motifinput[i] not in ValidBases:
        print("Non-DNA letter detected in motif. Aborting...")
        sys.exit()

motifchecker2 = input("What motif do you want to test the presence of? (can insert square brackets): ")
for i in range(0, len(motifchecker2)):                                            #Check if all DNA-bases. Soon, add [AG] .
    if motifchecker2[i] not in ValidBases2:
        print("Non-DNA letter detected in tested motif. Aborting...")
        sys.exit()
GCContent = input("What is your GC content, in a percentage?: ")

for i in range(0, len(GCContent)):
    if GCContent[i] not in validnumbers:
        print("Only numbers allowed; including float.")
        sys.exit()
if float(GCContent) >= 0 and float(GCContent) <= 100:
    pass
else:
    print('GC Content must be between 0 and 100%.')
    sys.exit()

TAContent = 100 - (float(GCContent))
GCFraction = float(GCContent)/100
TAFraction = float(TAContent)/100

motifchecker = motifchecker2.upper()
outputfile = open("dummydata1.txt", 'w')  #Wipe old, create new output file once checks are passed.

#random sequences with motif thrown in
sequencelist = []
motif = motifinput.upper() #input motif in capitals to match rest of sequence. Looks neater
TwoDArray = []

for i in range(0,20): #no 2D vector cos im lazy
    tempseq = []
    for u in range (0,(1000 - len(motif))): #leave space for adding motif (than sequence will fill to 1000 nt)
        randomnumber = random.random()
        if randomnumber <= (float(GCFraction)/2):
            tempseq.append("G")
        elif randomnumber <= (float(GCFraction)):
            tempseq.append("C")
        elif randomnumber <= (float(GCFraction) + (float(TAFraction)/2)):
            tempseq.append("T")
        elif randomnumber <= (1):
            tempseq.append("A")
    randomposition = randint(0,(1000 - len(motif)))
    outputseq = (''.join(tempseq))[:randomposition] + motif + (''.join(tempseq))[randomposition:] #Place motif within string version of random sequence at random location
    TwoDArray.append(''.join(outputseq))
    print(">Sequence", i, "\n", outputseq, sep = '', file= outputfile)

counter1 = 0
counter2 = 0
counter3 = 0
counter4 = 0
counter5 = 0
counter6 = 0 #this one is checking dodgy motif
for i in range(0, len(TwoDArray)):
    if motif in TwoDArray[i]: #check if motif in sequence. might pop up twice +. one is minimum.
        counter1 += len(re.findall(motif, TwoDArray[i]))
    if len(motif) > 4 and motif[:-1] in TwoDArray[i]: #if longer than 4, check if reduced length motif is also in sequence.
        counter2 += len(re.findall(motif[:-1], TwoDArray[i]))
    if len(motif) > 5 and motif[:-2] in TwoDArray[i]: #if longer than 5, check if reduced length motif is also in sequence.
        counter3 += len(re.findall(motif[:-2], TwoDArray[i]))
    if len(motif) > 6 and motif[:-3] in TwoDArray[i]: #if longer than 6, check if reduced length motif is also in sequence.
        counter4 += len(re.findall(motif[:-3], TwoDArray[i]))
    if len(motif) > 7 and motif[:-4] in TwoDArray[i]: #if longer than 7, check if reduced length motif is also in sequence.
        counter5 += len(re.findall(motif[:-4], TwoDArray[i]))
    counter6 += len(re.findall(motifchecker, TwoDArray[i]))

print("Number of matches of motif shortened motif to all sequences is: ",  counter1,
      "\nNumber of matches of 1-nt shortened motif to all sequences is: ", counter2,
      "\nNumber of matches of 2-nt shortened motif to all sequences is: ", counter3,
      "\nNumber of matches of 3-nt shortened motif to all sequences is: ", counter4,
      "\nNumber of matches of 4-nt shortened motif to all sequences is: ", counter5,
      "\nNumber of matches of test motif to all sequences is: ", counter6)

y = [counter1, counter2, counter3, counter4, counter5, counter6]
N = len(y)
x = range(N)
width = 1/1.5
plt.bar(x, y, width, color="red")
plt.xlabel("Motifs")
plt.ylabel("Number of matches")
plt.title("Testing specificity of motifs on randomly generated sequences")

plt.show()