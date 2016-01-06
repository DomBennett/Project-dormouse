
import re

with open("final_trees/combined_2segs.tre", 'rb') as f:
    newick = f.readline()

new_names = []
with open("2_alignments/supermatrix.phy", 'rb') as f:
    n = int(f.readline().split()[0])
    lines = f.readlines()[:n]

for line in lines:
    line = line.split()
    if line:
        new_names.append(line[0])

# do the ones printed manually

for name in new_names:
    sp = name.split('__')[0]
    finds = re.findall(sp, newick)
    if len(finds) > 1:
        print name
    else:
        robj = re.search(sp, newick)
        start = robj.start()
        robj = re.search(":", newick[start:])
        end = robj.end() + start - 1
        newick = newick[:start] + name + newick[end:]

with open("final_trees/combined_2segs_corrected.tre", 'w') as f:
    f.write(newick)
