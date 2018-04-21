__author__ = 'majid'
import sys

res = open(sys.argv[1])
gt = open(sys.argv[2])

pref1 = sys.argv[3]
pref2 = sys.argv[4]
seperator = sys.argv[5]

numcorrect = 0
numincorrect = 0
numnotcovered = 0

allgt = set()
covered = set()
for line in gt:
    args = line.strip().split()
    allgt.add((args[0], args[1]))

for line in res:
    if "@prefix" in line:
        continue
    args = line.strip().split(seperator)
    dom1 = args[0].replace(pref1,"").replace('> .','').replace('<','').replace('>','')
    dom2 = args[1].replace(pref1,"").replace('> .','').replace('<','').replace('>','')
    if dom1 != dom2:
        if (dom1,dom2) in allgt:
            numcorrect += 1
            covered.add((dom1,dom2))
        else:
            print((dom1,dom2))
            numincorrect += 1

numnotcovered = len(allgt-covered)
print(numcorrect)
recall = float(numcorrect)/float(len(allgt))
precision = float(numcorrect)/float(numcorrect+numincorrect)
fscore = 2.0*precision*recall/(precision+recall)
print("precision: " + str(precision))
print("recall: " + str(recall))
print("F-score: " + str(fscore))
