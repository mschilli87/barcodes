import time

start_time = time.time()

topDict = dict()
bulkDict = dict()
bcToCollapse = dict()

with open('topBarcodes.txt', 'r') as tB:
	for line in tB:
		topDict[line.strip()] = line.strip()

with open('restBarcodes.txt', 'r') as rB:
	for line in rB:
		bulkDict[line.strip()] = line.strip()

lex = ['A', 'C', 'G', 'T']

for s in topDict:
	if 'N' in s:
		continue
	for i in range(12):
		lex.remove(s[i])
		q = [a + b for a, b in zip([s[0:i], s[0:i], s[0:i]], lex)]
		q = [a + b for a, b in zip(q, [s[i+1:12], s[i+1:12], s[i+1:12]])]
		lex.append(s[i])
		if q[0] in bulkDict and q[0] not in bcToCollapse:
			bcToCollapse[q[0]] = s
		elif q[0] in bulkDict and q[0] in bcToCollapse:
			del bcToCollapse[q[0]]
			del bulkDict[q[0]]
		if q[1] in bulkDict and q[1] not in bcToCollapse:
			bcToCollapse[q[1]] = s
		elif q[1] in bulkDict and q[1] in bcToCollapse:
			del bcToCollapse[q[1]]
			del bulkDict[q[1]]
		if q[2] in bulkDict and q[2] not in bcToCollapse:
			bcToCollapse[q[2]] = s
		elif q[2] in bulkDict and q[2] in bcToCollapse:
			del bcToCollapse[q[2]]
			del bulkDict[q[2]]

with open ('collapseThem.txt', 'w') as fw:
	fw.writelines('{}\t{}\n'.format(k, v) for k, v in bcToCollapse.items())

print len(bcToCollapse), 'barcodes to collapse in', time.time() - start_time, 'seconds'







