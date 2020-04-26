import matplotlib.pyplot as plt
import numpy as np
import math

def entropy_list(file_name, window):
    span=window
    file_length=int(0)
    y=int(0)
    entropy_list=[]
    decimal_places=3

    file=open(file_name, "r")

    if file.mode == 'r':
        RNAstring=file.read()
        RNAstring=RNAstring.replace("\n","")
        RNAstring=RNAstring.replace("\r","")
        RNAstring=RNAstring.replace(" ","")
        file_length=len(RNAstring)



        for x in range(file_length-window+1):
            numA=float(0)
            numC=float(0)
            numT=float(0)
            numG=float(0)

            while y<span:
                if RNAstring[y]=='A':
                    numA+=1
                if RNAstring[y]=='C':
                    numC+=1
                if RNAstring[y]=='T':
                    numT+=1
                if RNAstring[y]=='G':
                    numG+=1
                y+=1
            pA=(numA/window)
            pC=(numC/window)
            pT=(numT/window)
            pG=(numG/window)

            if pA==0:
                A=0
            else:
                A=(pA*math.log(pA,2))

            if pC==0:
                C=0
            else:
                C=(pC*math.log(pC,2))

            if pT==0:
                T=0
            else:
                T=(pT*math.log(pT,2))

            if pG==0:
                G=0
            else:
                G=(pG*math.log(pG,2))
            entropy=-(A+C+T+G)
            entropy_list.insert(x,round(entropy, decimal_places))
            y=x+1
            span+=1
    return entropy_list

def generateX(ylist):
    Ylength = len(ylist)
    xlist = []

    for i in range(Ylength):
        xlist.insert(i, i)
    return xlist

def main():
    #1. this part will be changed to allow for dynamic input from the webpage GUI
    # input parameters: RNA sequence, window length, number of RNA sequences to compare
    #2. and for connection to chosen DB.
    # below sequences were downloaded from https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/

    filename1 = 'Covid19-RNAsequence.txt'
    #filename2 = 'Coronavirus229E-RNAsequence.txt'
    label1 = 'covid-19 w8'
    #label2 = 'covid-229E'
    label2 = 'covid-19 w16'

    entropyY1 = entropy_list(filename1, 8)
    entropyX1 = generateX(entropyY1)

    entropyX1_1000 = entropyX1[:100]
    entropyY1_1000 = entropyY1[:100]

    entropyY2 = entropy_list(filename1, 16)
    entropyX2 = generateX(entropyY2)

    entropyX2_1000 = entropyX2[:100]
    entropyY2_1000 = entropyY2[:100]

    fig, ax = plt.subplots()
    ax.plot(entropyX1_1000, entropyY1_1000, label = label1)
    ax.plot(entropyX2_1000, entropyY2_1000, label = label2)
    ax.set_xlabel('no of points')
    ax.set_ylabel('entropy points')
    ax.set_title('Virus RNA entropy plot')
    ax.legend()


    plt.show()

main()
