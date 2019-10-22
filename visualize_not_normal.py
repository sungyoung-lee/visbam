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
from sklearn.cluster import KMeans



'''
Argument Setting
'''

# 파일 이름과 span을 argument로 불러들인다.
parser = argparse.ArgumentParser()
parser.add_argument('bam_dir_path', help='bam파일을 읽어들일 디렉토리를 정합니다.')
parser.add_argument('sample_list_path', help='해당하는 sample 이름이 들어있는 경로를 지정합니다.')
parser.add_argument('refseq_path', help='refseq 파일의 경로를 지정합니다.')
parser.add_argument('nmid_to_draw', help='사용할 NMID를 지정합니다.')
parser.add_argument('draw_span', type=int, help='사진을 몇 bp단위로 분할할 것인지 정합니다.')
parser.add_argument('output_prefix', help='output 파일명을 정합니다.')
parser.add_argument('--color', help='선 색을 정합니다.(파일로 정의)')
parser.add_argument('--show_mean', help='주어진 sample들의 평균을 그립니다.', action='store_true')
parser.add_argument('--simplify', help='sample 전체를 그리지 않고 최댓값, 최솟값만 간단하게 표시합니다.', action='store_true')
parser.add_argument('-f','--flag', help='exon 주변의 100bp 부분만 그립니다.', action='store_true')
parser.add_argument('--flag3', help='curated genes', action='store_true')

args = parser.parse_args()
bam_dir = args.bam_dir_path
sample_list_path = args.sample_list_path
refseq_path = args.refseq_path
nmid_to_draw = args.nmid_to_draw
draw_span = args.draw_span
output_prefix = args.output_prefix
flag = args.flag
flag3 = args.flag3
color = args.color
show_mean = args.show_mean
simplify = args.simplify




# curated된 refseq만 분류하는 함수

def include(refseq) :
	# curated된 refseq 목록을 불러온다.
	curated_f = open('/data/allbam/asset/allCuratedGenes.txt', 'r')
	curated_r = curated_f.readlines()
	# 모든 값에 대해 false인 Series를 만든다.
	curated = refseq.contains("NNNNNNNNNNNNNNNNNNNNNN")
	# 한 줄씩 반복하면서 해당 유전자가 있는지 확인하고
	# OR 연산을 반복하여 return 한다.
	for cl in curated_r :
		crt = cl.split('\t')[1]
		crt = crt[:crt.find('.')]
		if crt == '' :
			continue
		curated = curated | refseq.contains(crt)
	return curated


'''
Reading Refseq Data
'''

# Refseq 불러오기
print('reading refseq data...')
refseq = pd.read_csv(refseq_path,
	sep='\t',
	names=['bin', 'name', 'chrom', 'strand',
		'txStart', 'txEnd', 'cdsStart', 'cdsEnd',
		'exonCount', 'exonStarts', 'exonEnds', 'score',
		'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames']
	)

# 사용자가 찾으려 하는 nmid의 refseq 정보만 불러온다.
refseq_nm = refseq[refseq.name.str.contains(nmid_to_draw)]

chrom_nm = refseq_nm.chrom.tolist()
tx_s_nm = refseq_nm.txStart.tolist()
tx_e_nm = refseq_nm.txEnd.tolist()
exon_s_nm = refseq_nm.exonStarts.tolist()
exon_e_nm = refseq_nm.exonEnds.tolist()

# 그 유전자를 토대로 그릴 범위를 정한다.
contig = chrom_nm[0]
start = tx_s_nm[0]
stop = tx_e_nm[0]


# 그릴 refseq 정보만 filtering 후 불러온다
refseq = refseq[refseq.name.str.contains("NM")]
refseq = refseq[refseq.txStart <= stop]
refseq = refseq[refseq.txEnd >= start]
if flag3 :
	refseq = refseq[include(refseq.name.str)]

# refseq 정보
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

# coverage를 저장할 데이터를 초기화한다.
coverage = [[] for i in range(stop-start+1)]
samfile = None


# Cancer Bam
print('\nanalyzing cancer bam information...')

# sample list가 들어있는 파일을 불러온다.
slfile = open(sample_list_path, 'r')
sl_ls = slfile.readlines()
bam_list = []
slfile.close()

real_bam_list = []

sdfs = False

# 목록을 불러온 뒤 유효한 파일 이름으로 바꾸어준다.
for sl_l in sl_ls :
	bam_list.append(sl_l[:-1]+'.bwamem.sorted.dedup.realn.recal.dedup.bam')




for bamn, bam in enumerate(bam_list) :

#	if sdfs :
#		break

	if bamn >= 20 :
		break

	sys.stderr.write('\r'+bam+':'+str(bamn+1)+'/'+str(len(bam_list))+': start...')
	sys.stderr.write("\033[K")
	# Cancer Bam 경로
	sam_path = bam_dir+'/'+bam	

	# Cancer Bam이 파일인지 확인
	if not os.path.isfile(sam_path) :
		continue
	
	real_bam_list.append(sl_ls[bamn][:-1])

	sdfs = True
	sys.stderr.write('\r'+bam+':'+str(bamn+1)+'/'+str(len(bam_list))+': loading bam file...')
	sys.stderr.write("\033[K")
	# Bam파일을 불러온다.
	samfile = pysam.AlignmentFile(sam_path, "rb")

	sys.stderr.write('\r'+bam+':'+str(bamn+1)+'/'+str(len(bam_list))+': coverage calculating...')
	sys.stderr.write("\033[K")
	# Coverage 계산
	assdfe = samfile.count_coverage(contig, start=start, stop=stop+1)

	sys.stderr.write('\r'+bam+':'+str(bamn+1)+'/'+str(len(bam_list))+': converting coverage to array...')
	sys.stderr.write("\033[K")
	cv_original = np.array(assdfe)
	samfile.close()
	sys.stderr.write('\r'+bam+':'+str(bamn+1)+'/'+str(len(bam_list))+': coverage adding...')
	sys.stderr.write("\033[K")
	# Coverage가  A, T, G, C 4성분으로 다 따로 출력되기 때문에,
	# 이를 합쳐주는 작업을 한다.
	cv = cv_original.sum(axis=0) 

	# coverage의 최댓값으로 나머지 값을 나누어주는 작업을 한다.
	# 그러면 전체 값의 범위가 0~1이 된다.
	# 전체 coverage의 최댓값
	maxcv = max(cv)	
	# 모든 위치에 대해 최댓값으로 나누어 준다.
	for j, out in enumerate(cv) :
		coverage[j].append(out/maxcv)




	

'''
Draw Lineplot
'''

# Lineplot을 그릴 위치를 정합니다.
# 각 그림의 시작
draw_range = []
# 각 그림의 끝
draw_range_e = []

# flag 속성이 켜져 있으면, 지정한 refseq의
# 각 exon 주위의 100bp 만큼의 범위를 그린다.
# 그렇지 않으면, refseq 전체 범위를 지정한 크기만큼 잘라 그린다.
if flag :	
	for j, e_s in enumerate(exon_s_nm) :
		# 각 exon의 start 부분
		ess = list(map(int, e_s[:-1].split(',')))
		# 각 exon의 stop 부분
		ees = list(map(int, exon_e_nm[j][:-1].split(',')))
		for k, es in enumerate(ess) :
			# refseq의 start와 비교해 그보다 작으면 start를 시작점으로 한다.
			if es-100 >= start :
				draw_range.append(es-100)
			else :
				draw_range.append(start)
			draw_range_e.append(ees[k]+100)
else :
	draw_range = range(start, stop+1, draw_span)

# 그래프 색 정보를 불러옵니다.
color_ls = []
if not color == None :
	color_f = open(color, 'r')
	color_ls = color_f.readlines()



# 위 배열에 해당하는 범위를 그래프로 그린다.
# for문을 시작점으로 돌린다.
for n, st in enumerate(draw_range) :

	print('\n'+output_prefix+'_'+str(n)+' saving...')

	# 그래프가 끝나는 지점을 정한다.
	# refseq가 끝나는 지점보다 크면 stop+1을 끝점으로 한다.
	stop_n = 0

	if flag :
		stop_n = stop+1 if draw_range_e[n] >= stop+1 else draw_range_e[n]
	else :
		stop_n = stop+1 if st+draw_span >= stop+1 else st+draw_span

	# x축의 값을 정해진 시작점과 끝점으로 한다.
	xticks = np.arange(st, stop_n)

	# refseq를 해당되는 부분에 포함되는 것만 선별한다.
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
	
	# refseq의 이름들
	nl = len(nmids)

	# Lineplot
	# figure 초기화(크기 조정, 간격 조정)
	fig2 = plt.figure(figsize=(30, 2+12+1*nl))
	gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[4+(nl-1)*0.2, nl])
	gs.update(wspace=0, hspace=0.05)


	# 색 설정
	print(real_bam_list)
	print('----------')
	colors = ['k']*len(real_bam_list)
	for color_l in color_ls :
		print(color_l)
		c_l = color_l.split('\t')
		if c_l[0] in real_bam_list :
			colors[real_bam_list.index(c_l[0])] = c_l[1][:-1]
	
	print(colors)

	# coverage의 dataframe을 만들고 plot의 윗 부분을 불러온다.
	df2 = pd.DataFrame(coverage[st-start:stop_n-start], index=xticks, columns=None)
	ax_main = plt.subplot(gs[0])
	

	'''
	Clustering Test
	'''
	
	df2_T = df2.T
	data_lines = df2_T.values
	kmeans = KMeans(n_clusters=2).fit(data_lines)

	df2_T['cluster_id'] = kmeans.labels_
	df2_c1 = df2_T[df2_T.cluster_id == 0].drop('cluster_id', 1).T
	df2_c2 = df2_T[df2_T.cluster_id == 1].drop('cluster_id', 1).T	

	colors = ['g' if lb == 1 else 'b' for lb in kmeans.labels_]

	# c1
	ax_main.fill_between(xticks, df2_c1.max(axis=1), df2_c1.min(axis=1), facecolor='b', alpha=0.5)
	ax_main.plot(df2_c1.max(axis=1), color='b', linewidth=3.0)
	ax_main.plot(df2_c1.min(axis=1), color='b', linewidth=3.0)
	ax_main.plot(df2_c1.mean(axis=1), color='red', linewidth=3.0)	

	# c2
	ax_main.fill_between(xticks, df2_c2.max(axis=1), df2_c2.min(axis=1), facecolor='g', alpha=0.5)
	ax_main.plot(df2_c2.max(axis=1), color='g', linewidth=3.0)
	ax_main.plot(df2_c2.min(axis=1), color='g', linewidth=3.0)
	ax_main.plot(df2_c2.mean(axis=1), color='red', linewidth=3.0)	
	
	# simplify 속성이 켜져 있으면
	# 전체 coverage 값에서 최댓값, 최솟값만 표시합니다.
	# 그렇지 않으면, coverage를 해당된 부분만 그린다.
#	if simplify :
#		ax_main.plot(df2.max(axis=1), color='black', linewidth=3.0)
#		ax_main.plot(df2.min(axis=1), color='black', linewidth=3.0)
#	else :
#		for rank, column in enumerate(df2.columns) :
#			ax_main.plot(df2[column], alpha=0.1, color=colors[rank])

	# show_mean 속성이 켜져 있을 경우 coverage의 평균을 그린다.
#	if show_mean :
#		ax_main.plot(df2.mean(axis=1), color='red', linewidth=3.0)	

	plt.ylim(0, 1) # coverage 표시 범위 설정
#	ax_main.set_yscale('log') # logscale 여부
	plt.xticks(np.arange(st, stop_n+1, step=(stop_n-st)/10)) # x축 값 표시.

	# bp 위치의 값이 너무 크기 때문에 e를 포함한 값으로 표시된다.
	# 따라서 아래와 같은 과정을 거쳐 자연수가 나오도록 한다.
	xx, locs = plt.xticks()
	ll = ['%d' % a for a in xx]
	plt.xticks(xx, ll)

	# 1 기준선을 그린다.
	reddot = np.ones(stop_n-st)
	ax_main.plot(xticks, reddot, 'r--')

	
	# refseq를 표시할 범위를 정한다.
	byts = range(2-nl, 2)

	# refseq가 표시될 이름을 정한다.
	ytlbs = [aa+"\n"+names[aai] for aai, aa in enumerate(nmids)]
	ax_bottom = plt.subplot(gs[1], yticks=byts, xticklabels=[], yticklabels=list(reversed(ytlbs)))
	ax_bottom.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False)

	# refseq의 범위를 검은 선으로 그린다.
	print('range of genes')
	for j, ts in enumerate(tx_s) :	
		a_s = ts if ts > st else st
		a_e = tx_e[j] if tx_e[j] < stop_n else stop_n
		if a_e-a_s < 0 :
			continue
		bxts = np.arange(a_s, a_e)
		blns = np.full(a_e-a_s,  1-j, dtype=int)
		ax_bottom.plot(bxts, blns, 'black')

	# tx, cds의 시작 부분을 얇은 네모로 그린다.
	print('tx, cds start')
	for j, cs in enumerate(cds_s) :	
			if (tx_s[j] > stop_n or cs < st) :
				continue
			rect = patches.Rectangle((tx_s[j], 0.9-j),cs-tx_s[j],0.2,edgecolor='none',facecolor='black')
			ax_bottom.add_patch(rect)

	# tx, cds의 끝 부분을 얇은 네모로 그린다.
	print('tx, cds end')
	for j, ce in enumerate(cds_e) :
			if (ce > stop_n or tx_e[j] < st) :
				continue
			rect = patches.Rectangle((ce, 0.9-j),tx_e[j]-ce,0.2,edgecolor='none',facecolor='black')
			ax_bottom.add_patch(rect)

	# 방향을 그린다.
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

	# 엑손을 그린다.
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

	# lineplot 간 여백 설정
	plt.subplots_adjust(top = 0.95, bottom = 0.05, right = 0.9, left = 0.1, wspace=0, hspace=0)
	# 폰트 사이즈 조정
	matplotlib.rcParams.update({'font.size': 22})

	# pdf로 저장한 뒤 닫는다.
	plt.savefig(output_prefix+'_'+str(n+1)+'.pdf')
	plt.close(fig2)
	print(output_prefix+'_'+str(n+1)+'.pdf saved!')




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



