__author__ = 'dilmiperera'
import csv
normalMotifsCSV='/Users/dilmiperera/Desktop/Lowy-Dilmi-PC/Scoring/Non-Path-7HUFK/Motif1/Possum-for-scoring.txt'
mutantMotifsCSV='/Users/dilmiperera/Desktop/Lowy-Dilmi-PC/Scoring/Non-Path-7HUFK/Motif2/Possum-for-scoring.txt'
outputFile=open('/Users/dilmiperera/Desktop/Lowy-Dilmi-PC/Scoring/Non-Path-7HUFK/MotifScore.bed', 'w')
normalMotifs={}
mutantMotifs={}
seqs=[]

#Format files
#
# with open(normalMotifsCSV, 'rU') as csv_file:
#     csv_reader1 = csv.reader(csv_file, delimiter='\t')
#     seq=""
#     motif=""
#     start=""
#     end=""
#     strand=""
#     score=""
#     for row in csv_reader1:
#         if str(row[0]).__contains__(">"):
#             normalMotifs.append([seq,motif,start,end,strand,score])
#             seq=str(row[0])
#             motif=""
#             start=""
#             end=""
#             strand=""
#             score=""
#         else :
#             motif=motif+","+str(row[0])
#             start=start+","+str(row[1])
#             end=end+","+str(row[2])
#             strand=strand+","+str(row[3])
#             score=score+","+str(row[5])
#
# with open(mutantMotifsCSV, 'rU') as csv_file:
#     csv_reader1 = csv.reader(csv_file, delimiter='\t')
#     seq=""
#     motif=""
#     start=""
#     end=""
#     strand=""
#     score=""
#     for row in csv_reader1:
#         if str(row[0]).__contains__(">"):
#             mutantMotifs.append([seq,motif,start,end,strand,score])
#             seq=str(row[0])
#             motif=""
#             start=""
#             end=""
#             strand=""
#             score=""
#         else :
#             motif=motif+","+str(row[0])
#             start=start+","+str(row[1])
#             end=end+","+str(row[2])
#             strand=strand+","+str(row[3])
#             score=score+","+str(row[5])


with open(normalMotifsCSV, 'rU') as csv_file:
    csv_reader1 = csv.reader(csv_file, delimiter='\t')
    seq=""
    motifs=[[]]
    for row in csv_reader1:
        if str(row[0]).__contains__(">"):
            if seq != "":
                normalMotifs[seq]=motifs
                seqs.append(seq)
            seq=str(row[0])

            motifs=[]
        else :
            motifs.append([str(row[0]),str(row[1]),str(row[2]),str(row[3]),float(row[5])])


with open(mutantMotifsCSV, 'rU') as csv_file:
    csv_reader1 = csv.reader(csv_file, delimiter='\t')
    seq=""
    motifs=[[]]
    for row in csv_reader1:
        if str(row[0]).__contains__(">"):
            if seq != "":
                mutantMotifs[seq]=motifs

            seq=str(row[0])
            motifs=[]
        else :
            motifs.append([str(row[0]),str(row[1]),str(row[2]),str(row[3]),float(row[5])])

for seq in seqs:
    normalMotifList=normalMotifs.get(seq)
    mutantMotifList=mutantMotifs.get(seq)
    max=0
    if len(normalMotifs.get(seq))==0 and len(mutantMotifs.get(seq))==0:
        continue

    for motifN in normalMotifList:
        found=False
        for motifM in mutantMotifList:
            if motifN[0]==motifM[0] and motifN[1]==motifM[1] and motifN[2]==motifM[2] and motifN[3]==motifM[3]:
                diff=float(motifN[4]) - float(motifN[4])
                score=float(motifN[4]) * abs(diff)
                if max<score:
                    max=score
                found=True
                break;

        if found==False:
            score=float(motifN[4])*float(motifN[4])
            if max<score:
                max=score


    for motifM in mutantMotifList:
        found=False
        for motifN in normalMotifList:
            if motifN[0]==motifM[0] and motifN[1]==motifM[1] and motifN[2]==motifM[2] and motifN[3]==motifM[3]:
                found=True
                break;

        if found==False:
            score=float(motifM[4])*float(motifM[4])
            if max<score:
                max=score

    outputFile.write(seq+"\t"+str(max)+"\n")
    print(seq+"\t"+str(max))

outputFile.close()
