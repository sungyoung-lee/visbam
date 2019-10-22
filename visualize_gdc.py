import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import gridspec
from operator import add
import pysam, argparse




'''
Argument Setting
'''

# 파일 이름과 span을 argument로 불러들인다.
parser = argparse.ArgumentParser()
parser.add_argument('bam_dir_path', help='bam파일을 읽어들일 디렉토리를 정합니다.')
parser.add_argument('sample_list_path', help='해당하는 sample 이름이 들어있는 경로를 지정합니다.')
parser.add_argument('normal_dir_path', help='normal sample이 들어있는 경로를 지정합니다.')
parser.add_argument('refseq_path', help='refseq 파일의 경로를 지정합니다.')
parser.add_argument('nmid_to_draw', help='사용할 NMID를 지정합니다.')
parser.add_argument('draw_span', type=int, help='사진을 몇 bp단위로 분할할 것인지 정합니다.')
parser.add_argument('output_prefix', help='output 파일명을 정합니다.')
parser.add_argument('-f','--flag', help='cds 주변만 그립니다.', action='store_true')

args = parser.parse_args()
bam_dir = args.bam_dir_path
sample_list_path = args.sample_list_path
normal_dir_path = args.normal_dir_path
refseq_path = args.refseq_path
nmid_to_draw = args.nmid_to_draw
draw_span = args.draw_span
output_prefix = args.output_prefix
flag = args.flag



'''
Reading Refseq Data
'''

# Refseq
print('reading refseq data...')
refseq = pd.read_csv(refseq_path,
	sep='\t',
	names=['bin', 'name', 'chrom', 'strand',
		'txStart', 'txEnd', 'cdsStart', 'cdsEnd',
		'exonCount', 'exonStarts', 'exonEnds', 'score',
		'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames']
	)

refseq_nm = refseq[refseq.name.str.contains(nmid_to_draw)]

chrom_nm = refseq_nm.chrom.tolist()
tx_s_nm = refseq_nm.txStart.tolist()
tx_e_nm = refseq_nm.txEnd.tolist()
exon_s_nm = refseq_nm.exonStarts.tolist()
exon_e_nm = refseq_nm.exonEnds.tolist()

contig = chrom_nm[0]
start = tx_s_nm[0]
stop = tx_e_nm[0]


refseq = refseq[refseq.name.str.contains("NM")]
refseq = refseq[refseq.txStart <= stop]
refseq = refseq[refseq.txEnd >= start]

chrom = refseq.chrom.tolist()
strands = refseq.strand.tolist()
tx_s = refseq.txStart.tolist()
tx_e = refseq.txEnd.tolist()
cds_s = refseq.cdsStart.tolist()
cds_e = refseq.cdsEnd.tolist()
exon_s = refseq.exonStarts.tolist()
exon_e = refseq.exonEnds.tolist()
nmids = refseq.name.tolist()
names = refseq.name2.tolist()


print('there are '+str(len(names))+' refseq datas')





'''
Bam Information Analysis
'''

normal_coverage = np.zeros(stop-start+1)
coverage = [[] for i in range(stop-start+1)]
samfile = None

# Normal Bam
print('analyzing normal bam information...')

bam_list = os.listdir(normal_dir_path)
bam_list = [file for file in bam_list if file.endswith(".bam")]

for bam in bam_list :

	print('\r', bam, end='')
	sys.stdout.write("\033[K")

	sam_path = normal_dir_path+'/'+bam
	
	if not os.path.isfile(sam_path) :
		continue

	samfile = pysam.AlignmentFile(sam_path)
	print('loaded')

	cv_original = np.array(samfile.count_coverage(contig, start=start, stop=stop+1))
	print('calculated')
	samfile.close()
	cv = cv_original.sum(axis=0) 

	normal_coverage = list(map(add, normal_coverage, cv))


normal_coverage = [x / len(bam_list) for x in normal_coverage]


# Cancer Bam
print('\nanalyzing cancer bam information...')

slfile = open(sample_list_path, 'r')
sl_ls = slfile.readlines()
bam_list = []
slfile.close()


sdfs = False

for sl_l in sl_ls :
	bam_list.append(sl_l[:-1]+'.bwamem.sorted.dedup.realn.recal.dedup.bam')



bamd_list = os.listdir(bam_dir)
print(bamd_list)
bamd_list = [file for file in bamd_list if os.path.isdir(bam_dir+'/'+file)]
print(bamd_list)
bam_list = []
for ddd in bamd_list :
	ffff = os.listdir(bam_dir+'/'+ddd)
	ffff = [ddd+'/'+file for file in ffff if file.endswith(".bam")]
	bam_list.extend(ffff)

print(bam_list)



for bamn, bam in enumerate(bam_list) :

#	if sdfs :
#		break

	sys.stderr.write('\r'+bam+':'+str(bamn)+'/'+str(len(bam_list))+': start...')
	sys.stderr.write("\033[K")
	sam_path = bam_dir+'/'+bam	

#	if not os.path.isfile(sam_path) :
#		continue

	sdfs = True
	sys.stderr.write('\r'+bam+':'+str(bamn)+'/'+str(len(bam_list))+': loading bam file...')
	samfile = pysam.AlignmentFile(sam_path, "rb")

	sys.stderr.write('\r'+bam+':'+str(bamn)+'/'+str(len(bam_list))+': coverage calculating...')
	assdfe = samfile.count_coverage(contig, start=start, stop=stop+1)

	sys.stderr.write('\r'+bam+':'+str(bamn)+'/'+str(len(bam_list))+': converting coverage to array...')
	cv_original = np.array(assdfe)
	samfile.close()
	sys.stderr.write('\r'+bam+':'+str(bamn)+'/'+str(len(bam_list))+': coverage adding...')
	cv = cv_original.sum(axis=0) 

	sys.stderr.write('\r'+bam+':'+str(bamn)+'/'+str(len(bam_list))+': coverage dividig by normal coverage...')
	for j, out in enumerate(cv) :
		cov = 1
		if normal_coverage[j] != 0 :
			cov = out/normal_coverage[j]
		coverage[j].append(cov)	



'''
Draw Lineplot
'''

draw_range = []
draw_range_e = []

if flag :	
	for j, e_s in enumerate(exon_s_nm) :
		ess = list(map(int, e_s[:-1].split(',')))
		ees = list(map(int, exon_e_nm[j][:-1].split(',')))
		for k, es in enumerate(ess) :
			if es-100 >= start :
				draw_range.append(es-100)
			else :
				draw_range.append(start)
			draw_range_e.append(ees[k]+100)
else :
	draw_range = range(start, stop+1, draw_span)



for n, st in enumerate(draw_range) :

	print('\n'+output_prefix+'_'+str(n)+' saving...')

	stop_n = 0

	if flag :
		stop_n = stop+1 if draw_range_e[n] >= stop+1 else draw_range_e[n]
	else :
		stop_n = stop+1 if st+draw_span >= stop+1 else st+draw_span

	xticks = np.arange(st, stop_n)

	refseq_r = refseq[refseq.txStart <= stop_n]
	refseq_r = refseq_r[refseq_r.txEnd >= st]
		
	chrom = refseq_r.chrom.tolist()
	strands = refseq_r.strand.tolist()
	tx_s = refseq_r.txStart.tolist()
	tx_e = refseq_r.txEnd.tolist()
	cds_s = refseq_r.cdsStart.tolist()
	cds_e = refseq_r.cdsEnd.tolist()
	exon_s = refseq_r.exonStarts.tolist()
	exon_e = refseq_r.exonEnds.tolist()
	nmids = refseq_r.name.tolist()
	names = refseq_r.name2.tolist()
	
	nl = len(nmids)
	# Lineplot
	fig2 = plt.figure(figsize=(30, 2+12+1*nl))
	gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[4+(nl-1)*0.2, nl])
	gs.update(wspace=0, hspace=0.05)

	print(start, st, stop_n)
	df2 = pd.DataFrame(coverage[st-start:stop_n-start], index=xticks, columns=None)
	ax_main = plt.subplot(gs[0])
	ax_main.plot(df2, color='black', alpha=0.1)
#	plt.ylim(0, 1.5)
#	ax_main.set_yscale('log')
	plt.xticks(np.arange(st, stop_n+1, step=(stop_n-st)/10))

	xx, locs = plt.xticks()
	ll = ['%d' % a for a in xx]
	plt.xticks(xx, ll)

	reddot = np.ones(stop_n-st)
	ax_main.plot(xticks, reddot, 'r--')

	byts = range(2-nl, 2)

	ytlbs = [aa+"\n"+names[aai] for aai, aa in enumerate(nmids)]

	ax_bottom = plt.subplot(gs[1], yticks=byts, xticklabels=[], yticklabels=list(reversed(ytlbs)))
	ax_bottom.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False)

	print('range of genes')
	for j, ts in enumerate(tx_s) :	
		a_s = ts if ts > st else st
		a_e = tx_e[j] if tx_e[j] < stop_n else stop_n
		if a_e-a_s < 0 :
			continue
		bxts = np.arange(a_s, a_e)
		blns = np.full(a_e-a_s,  1-j, dtype=int)
		ax_bottom.plot(bxts, blns, 'black')

	print('tx, cds start')
	for j, cs in enumerate(cds_s) :	
			if (tx_s[j] > stop_n or cs < st) :
				continue
			rect = patches.Rectangle((tx_s[j], 0.9-j),cs-tx_s[j],0.2,edgecolor='none',facecolor='black')
			ax_bottom.add_patch(rect)

	print('tx, cds end')
	for j, ce in enumerate(cds_e) :
			if (ce > stop_n or tx_e[j] < st) :
				continue
			rect = patches.Rectangle((ce, 0.9-j),tx_e[j]-ce,0.2,edgecolor='none',facecolor='black')
			ax_bottom.add_patch(rect)

	print('draw directions...')
	for j, ts in enumerate(tx_s) :
		strand = strands[j]
		if (ts > stop_n or tx_e[j] < st) :
			continue
		interval = int((stop_n-st)/60)
		a_s = ts if ts > st else st
		a_e = tx_e[j] if tx_e[j] < stop_n else stop_n 
		for k in range(a_s, a_e, interval) :
			if strand == '+' :
				ax_bottom.arrow(k, 1-j, interval, 0, head_width=0.2, head_length=interval/2, overhang=1)
			else :	
				ax_bottom.arrow(k, 1-j, interval*(-1), 0, head_width=0.2, head_length=interval/2, overhang=1)

	print('exons')
	for j, e_s in enumerate(exon_s) :
		ess = list(map(int, e_s[:-1].split(',')))
		ees = list(map(int, exon_e[j][:-1].split(',')))
		for k, es in enumerate(ess) :
			if (es > stop_n or ees[k] < st) :
				continue
			rect = patches.Rectangle((es, 0.8-j),ees[k]-es,0.4,edgecolor='none',facecolor='black')
			ax_bottom.add_patch(rect)
			leftt = es if es > st else st
			rightt = ees[k] if ees[k] < stop_n else stop_n
			ax_bottom.text((leftt+rightt)/2, 1-j, str(k+1), horizontalalignment='center', verticalalignment='center', color='white')

	plt.subplots_adjust(top = 0.95, bottom = 0.05, right = 0.9, left = 0.1, wspace=0, hspace=0)
	matplotlib.rcParams.update({'font.size': 22})

	plt.savefig(output_prefix+'_'+str(n)+'.pdf')
	plt.close(fig2)
	print(output_prefix+'_'+str(n)+'.pdf saved!')

'''

for i in range(55222713-start, 55223713-start, 100) :
	print(i+start, coverage[i])


for j, e_s in enumerate(exon_s) :
	ess = list(map(int, e_s[:-1].split(',')))
	ees = list(map(int, exon_e[j][:-1].split(',')))
	for k, es in enumerate(ess) :
		print(k, es, ees[k])
'''
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



