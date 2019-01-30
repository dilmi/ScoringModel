__author__ = 'dilmiperera'
import csv
IDfile ='/Users/dilmiperera/Desktop/Lowy-Dilmi-PC/Scoring/Non-Path-7HUFK/Mutation-20bpExtended.bed'
motifScores='/Users/dilmiperera/Desktop/Lowy-Dilmi-PC/Scoring/Non-Path-7HUFK/MotifScore.bed'
outputFile=open('/Users/dilmiperera/Desktop/Lowy-Dilmi-PC/Scoring/Non-Path-7HUFK/MotifScore-mapped-to-ID.bed', 'w')
found=False

scores={}

with open(motifScores, 'rU') as csv_file2:
    csv_reader2 = csv.reader(csv_file2, delimiter='\t')
    for row in csv_reader2:
        scores[str(row[3])]=str(row[4])


with open(IDfile, 'rU') as csv_file:
    csv_reader1 = csv.reader(csv_file, delimiter='\t')
    for row in csv_reader1:
        if str(row[3]) in scores:
            outputFile.write(str(row[0])+"\t"+str(row[1])+"\t"+str(row[2])+"\t"+str(row[3])+"\t"+scores[str(row[3])]+"\n")

        else:
            content=str(row[0])+"\t"+str(row[1])+"\t"+str(row[2])+"\t"+str(row[3])+"\t"+"0"+"\n"
            outputFile.write(content)



outputFile.close()
