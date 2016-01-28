import sys
import time
import itertools
from collections import defaultdict

filename = sys.argv[1]

d = defaultdict(list)
keepThem = dict()
oneDeletion = dict()
twoDeletion = dict()
threeDeletion = dict()

start_time = time.time()

with open('all.txt', 'r') as allBC:
	for line in allBC:
		key, val = line.split()
		if key not in d:
			d[key] = [val]
		elif key in d:
			d[key].append(val)

print len(d), 'unique barcodes found'
print 'dictionary read in', round(time.time() - start_time, 2), 'seconds'

start_time = time.time()

for barcode in d:
	counts = 0
	for umi in d[barcode]:
		if umi[-1] == 'T':
			counts += 1
	if len(d[barcode]) >= 40 and float(counts)/len(d[barcode]) > 0.8:
		counts = 0
		for umi in d[barcode]:
			if umi[-2] == 'T':
				counts += 1
		if float(counts)/len(d[barcode]) > 0.8:
			counts = 0
			for umi in d[barcode]:
				if umi[-3] == 'T':
					counts += 1
			if float(counts)/len(d[barcode]) > 0.8:
				threeDeletion[barcode] = barcode
				continue
			else:
				twoDeletion[barcode] = barcode
				continue
		else:
			oneDeletion[barcode] = barcode
			continue
	keepThem[barcode] = barcode

print '\n'
print len(oneDeletion), 'barcodes with one deletion'
print len(twoDeletion), 'barcodes with two deletions'
print len(threeDeletion), 'barcodes with three deletions'
print len(keepThem), 'barcodes kept intact'
print '...analysis took', round(time.time() - start_time, 2), 'seconds'

start_time = time.time()

lex = ['A', 'C', 'G', 'T']

collapseThem = defaultdict(list)

for barcode in threeDeletion:
	notFoundinKeepThem = True
	notFoundinOneDeletion = True
	notFoundinTwoDeletion = True
	for threeBasesExtra in itertools.product(lex, repeat=2):
		if barcode[0:9] + ''.join(threeBasesExtra) in keepThem:
			collapseThem[barcode] = [barcode[0:10] + ''.join(threeBasesExtra), 3]
			notFoundinKeepThem = False
			break
	if notFoundinKeepThem:
		for threeBasesExtra in itertools.product(lex, repeat=2):
			bcTwoDel = barcode[0:9] + ''.join(threeBasesExtra)
			if bcTwoDel in oneDeletion:
				collapseThem[barcode] = [bcTwoDel[0:10] + 'GC', 3]
				notFoundinOneDeletion = False
				break
	if notFoundinKeepThem and notFoundinOneDeletion:
		for threeBasesExtra in itertools.product(lex, repeat=2):
			bcOneDel = barcode[0:9] + ''.join(threeBasesExtra)
			if bcOneDel in oneDeletion:
				collapseThem[barcode] = [bcOneDel[0:11] + 'C', 3]
				notFoundinOneDeletion = False
				break
	if notFoundinKeepThem and notFoundinTwoDeletion and notFoundinOneDeletion:
		collapseThem[barcode] = [barcode[0:9] + 'AGC', 3]

for barcode in twoDeletion:
	notFoundinKeepThem = True
	notFoundinOneDeletion = True
	for twoBasesExtra in itertools.product(lex, repeat=2):
		if barcode[0:10] + ''.join(twoBasesExtra) in keepThem:
			collapseThem[barcode] = [barcode[0:10] + ''.join(twoBasesExtra), 2]
			notFoundinKeepThem = False
			break
	if notFoundinKeepThem:
		for twoBasesExtra in itertools.product(lex, repeat=2):
			bcOneDel = barcode[0:10] + ''.join(twoBasesExtra)
			if bcOneDel in oneDeletion:
				collapseThem[barcode] = [bcOneDel[0:11] + 'C', 2]
				notFoundinOneDeletion = False
				break
	if notFoundinOneDeletion and notFoundinKeepThem:
		collapseThem[barcode] = [barcode[0:10] + 'GC', 2]

for barcode in oneDeletion:
	notFoundinKeepThem = True
	for oneBaseExtra in lex:
		if barcode[0:11] + oneBaseExtra in keepThem:
			collapseThem[barcode] = [barcode[0:11] + oneBaseExtra, 1]
			notFoundinKeepThem = False
			break
	if notFoundinKeepThem:
		collapseThem[barcode] = [barcode[0:11] + 'C', 1]


print '\n'
print len(collapseThem), 'barcodes collapsed'
print '...analysis took', round(time.time() - start_time, 2), 'seconds'


# sys.exit()

# with open('keepThem.txt', 'w') as f:
# 	for barcode in keepThem:
# 		f.write(barcode + '\n')

# with open('oneDeletion.txt', 'w') as f:
# 	for barcode in oneDeletion:
# 		f.write(barcode + '\n')

# with open('twoDeletion.txt', 'w') as f:
# 	for barcode in twoDeletion:
# 		f.write(barcode + '\n')


start_time = time.time()

print '\n'
print 'changing fastq file... takes some time...'

fi = open(filename + '_1.fastq', 'r')
fo = open(filename + '_collapsed_1.fastq', 'w')

for line in fi:
	if line[0:12] in collapseThem:
		newBarcode = collapseThem[line[0:12]][0]
		dels = collapseThem[line[0:12]][1]
		fo.write(newBarcode + line[(12-dels):(20-dels)] + '\n')
		continue
	fo.write(line)

print '...changing fastq file took', round(time.time() - start_time, 2), 'seconds'




