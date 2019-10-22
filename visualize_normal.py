import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import gridspec
from operator import add

roi_path = 'egfr.roi.txt'
path_dir = 'out'
normal_dir = 'out_normal'
file_list = os.listdir(path_dir)
normal_list = os.listdir(normal_dir)

roifile = open(roi_path, "r")

lines = roifile.readlines()

for i, line in enumerate(lines) :
	first = line.find(':')
	second = line.find('-')
	contig = line[:first]
	start = int(line[first+1:second])
	stop = int(line[second+1:])
	normal_coverage = np.zeros(stop-start+1)
	coverage = [[] for i in range(stop-start+1)]

	# Refseq
	print('reading refseq data...')
	refseq = pd.read_csv("visibleData.bed",
		sep='\t',
		names=['bin', 'name', 'chrom', 'strand',
			'txStart', 'txEnd', 'cdsStart', 'cdsEnd',
			'exonCount', 'exonStarts', 'exonEnds', 'score',
			'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames']
		)

	strands = refseq.strand.tolist()
	tx_s = refseq.txStart.tolist()
	tx_e = refseq.txEnd.tolist()
	cds_s = refseq.cdsStart.tolist()
	cds_e = refseq.cdsEnd.tolist()
	exon_s = refseq.exonStarts.tolist()
	exon_e = refseq.exonEnds.tolist()
	names = refseq.name2.tolist()

	# Normal Bam
	print('calculating average of coverages of normal bam samples...')
	for normal_name in normal_list :
		outlist = list(map(int, open(normal_dir+"/"+normal_name, "r").readlines()[i].split(',')))
		normal_coverage = list(map(add, normal_coverage, outlist))
	
	normal_coverage = [x / len(normal_list) for x in normal_coverage]

	# Cancer Bam
	for file_name in file_list :
		outlist = list(map(int, open(path_dir+"/"+file_name, "r").readlines()[i].split(',')))
		for j, out in enumerate(outlist) :
			cov = 1
			if normal_coverage[j] != 0 :
				cov = out/normal_coverage[j]
			coverage[j].append(cov)
			print('\r', file_name+", roi"+str(i+1), end='')
	
	
	print('\n'+"roi"+str(i+1)+" saving...")

	xticks = np.arange(start, stop+1)

	# Lineplot
	fig2 = plt.figure()
	gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[4, 1])

	df2 = pd.DataFrame(coverage, index=xticks, columns=None)
	ax_main = plt.subplot(gs[0])
	ax_main.plot(df2, color='black', alpha=0.1)
	plt.xticks(np.arange(start, stop+1, step=(stop-start+1)/5-1))

	xx, locs = plt.xticks()
	ll = ['%d' % a for a in xx]
	plt.xticks(xx, ll)

	reddot = np.ones(stop-start+1)
	ax_main.plot(xticks, reddot, 'r--')

	ax_bottom = plt.subplot(gs[1], xticklabels=[], yticklabels=[])
	ax_bottom.plot(xticks, reddot, 'black')
	ax_bottom.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False)

	print('tx, cds start')
	for j, cs in enumerate(cds_s) :
			print(cs, tx_s[j])
			rect = patches.Rectangle((cs, 0.975),tx_s[j]-cs,0.05,edgecolor='none',facecolor='black')
			ax_bottom.add_patch(rect)

	print('tx, cds end')
	for j, ce in enumerate(cds_e) :
			print(ce, tx_e[j])
			rect = patches.Rectangle((tx_e[j], 0.975),ce-tx_e[j],0.05,edgecolor='none',facecolor='black')
			ax_bottom.add_patch(rect)

	print('draw directions...')
	for j, ts in enumerate(tx_s) :
		strand = strands[j]
		if (ts > stop or tx_e[j] < start) :
			continue
		interval = int((stop-start+1)/40)
		a_s = ts if ts > start else start
		a_e = tx_e[j] if tx_e[j] < stop else stop 
		for k in range(a_s, a_e, interval) :
			if strand == '+' :
				ax_bottom.arrow(k, 1, interval, 0, head_width=0.03, head_length=interval/2, overhang=1)
			else :	
				ax_bottom.arrow(k, 1, interval*(-1), 0, head_width=0.03, head_length=interval/2, overhang=1)

	print('exons')
	for j, e_s in enumerate(exon_s) :
		ess = list(map(int, e_s[:-1].split(',')))
		ees = list(map(int, exon_e[j][:-1].split(',')))
		for k, es in enumerate(ess) :
			if (es > stop or ees[k] < start) :
				continue
			print(es, ees[k])
			rect = patches.Rectangle((es, 0.96),ees[k]-es,0.08,edgecolor='none',facecolor='black')
			ax_bottom.add_patch(rect)
			leftt = es if es > start else start
			rightt = ees[k] if ees[k] < stop else stop
			ax_bottom.text((leftt+rightt)/2, 1, names[k], horizontalalignment='center', verticalalignment='center', color='white')

	plt.savefig(roi_path+"_"+"roi"+str(i+1)+"_normal_lineplot"+'.png')
	plt.close(fig2)
	print(roi_path+"_"+"roi"+str(i+1)+"_normal_lineplot"+'.png saved!')


'''
	# Boxplot
	fig = plt.figure()
	xticks = np.arange(start, stop+1)
	df = pd.DataFrame(list(map(list, zip(*coverage))))
	boxplot = df.boxplot()
	plt.savefig(roi_path+"_"+"roi"+str(i+1)+"_normal_boxplot"+'.png')
	plt.close(fig)
	print(roi_path+"_"+"roi"+str(i+1)+"_normal_boxplot"+'.png saved!')
'''

roifile.close()


