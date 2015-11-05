import sys

if len(sys.argv) < 3:
    print "\nUSAGE: python splitTable.py [filename] [number of columns]\n"
    print "For instance: \npython splitTable.py /scratch/mmasciov/analysisCode_scan/analysis/EventYields_mc_PHYS14_v5_dummy_3fb/yieldtable_full_HT1000to1500.log 6\n"
    print "where the first argument is the text file containing the yield table, and 6 is the max number of TR's I want to display on a single line\n"
    sys.exit(1)

filename=str(sys.argv[1])
maxline=int(sys.argv[2])

lines=open(filename).readlines()

for index, line in enumerate(lines):
    if line=="" or line=="\n" or ("end" in line):
        continue
    elif "hline" in line or ("begin" in line) or ("centering") in line or ("caption") in line:
        print line
    elif "multicolumn" in line:
        htl=index+1
        print line
    else:
        yields=line.split("&")
        for i in range( maxline+1 ):
            if i < maxline:
                print yields[i], " & ",
            else:
                print yields[i], "\\\\\n",

print "\\hline"
for index, line in enumerate(lines):
    if index < htl+1 or line=="" or line=="\n":
        continue
    elif "hline" in line or ("begin" in line) or ("end" in line) or ("centering") in line or ("caption") in line:
        print line
    else:
        yields=line.split("&")
        print yields[0], " & ",
        for i in xrange(maxline+1, len(yields)+len(yields)/2-maxline+1):
            if i < len(yields)-1:
                print yields[i], " & ",
            elif i == len(yields)-1: 
                print yields[i].replace("\\\\\n", ""),
            else:
                print " & ",
        print "\\\\\n",
