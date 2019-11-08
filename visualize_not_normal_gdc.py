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
parser.add_argument('refseq_path', help='refseq 파일의 경로를 지정합니다.')
parser.add_argument('nmid_to_draw', help='사용할 NMID를 지정합니다.')
parser.add_argument('draw_span', type=int, help='사진을 몇 bp단위로 분할할 것인지 정합니다.')
parser.add_argument('output_prefix', help='output 파일명을 정합니다.')
parser.add_argument('-f','--flag', help='exon 주변의 100bp 부분만 그립니다.', action='store_true')
parser.add_argument('--flag3', help='curated genes', action='store_true')
parser.add_argument('--exclude', default='', help='주어진 exon을 제외하고 계산&출력합니다.(1,2,3,4,...)')
parser.add_argument('--combine', help='모든 exon을 붙여 출력합니다.', action='store_true')
parser.add_argument('--view_mode', type=int, default=0,
			help='0 : 아무 표시도 하지 않음\n1 : 중앙값 표시\n2 : select exon에서 떨어지는 sample 표시')
parser.add_argument('--select', default='', help='view_mode에서 select할 exon 정의')
parser.add_argument('--threshold', type=float, default=1.0, help='view_mode 2에서 일정 범위 이상의 sample 제외(0~1)')

args = parser.parse_args()
bam_dir = args.bam_dir_path
refseq_path = args.refseq_path
nmid_to_draw = args.nmid_to_draw
draw_span = args.draw_span
output_prefix = args.output_prefix
flag = args.flag
flag3 = args.flag3
exclude = args.exclude.strip()
combine = args.combine
view_mode = args.view_mode
threshold = args.threshold
select = args.select

exclude_list = [] if exclude == '' else reversed(sorted(list(map(int, exclude.split(',')))))
select_list = [1, 7] if select == '' else list(map(int, select.split(',')))

# curated된 refseq만 분류하는 함수

def include(refseq) :
	# curated된 refseq 목록을 불러온다.
	curated_f = open('allCuratedGenes.txt', 'r')
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

name_nm = refseq_nm.name.tolist()
chrom_nm = refseq_nm.chrom.tolist()
tx_s_nm = refseq_nm.txStart.tolist()
tx_e_nm = refseq_nm.txEnd.tolist()
exon_s_nm = refseq_nm.exonStarts.tolist()
exon_e_nm = refseq_nm.exonEnds.tolist()

# 그 유전자를 토대로 그릴 범위를 정한다.
name_n = name_nm[0]
contig = chrom_nm[0]
start = tx_s_nm[0]
stop = tx_e_nm[0]
ess_nm = list(map(int, exon_s_nm[0].split(',')[:-1]))
ees_nm = list(map(int, exon_e_nm[0].split(',')[:-1]))

#include_list = [i for i in range(len(ess_nm)) if not i in exclude_list]

if not exclude == '' :
	stop = 0
	for i in exclude_list :
		if i <= len(ess_nm) :
			del ess_nm[i-1]
			del ees_nm[i-1]
	stop = ees_nm[-1]+100	


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


bamd_list = os.listdir(bam_dir)
bamd_list = [file for file in bamd_list if os.path.isdir(bam_dir+'/'+file)]
bam_list = []
for ddd in bamd_list :
	ffff = os.listdir(bam_dir+'/'+ddd)
	ffff = [ddd+'/'+file for file in ffff if file.endswith(".bam")]
	bam_list.extend(ffff)


for bamn, bam in enumerate(bam_list) :

#	if sdfs :
#		break


	sys.stderr.write('\r'+bam+':'+str(bamn+1)+'/'+str(len(bam_list))+': start...')
	sys.stderr.write("\033[K")
	# Cancer Bam 경로
	sam_path = bam_dir+'/'+bam	

	# Cancer Bam이 파일인지 확인
	if not os.path.isfile(sam_path) :
		continue

	# cache path
	bam_fn = os.path.basename(bam)
	cache_path = 'cache/'+bam_fn+'_'+name_n+'_'+contig+'_'+str(start)+'_'+str(stop)
	
	cv = []

	if not os.path.isfile(cache_path+'.npy') :

		# Bam파일을 불러온다.
		sys.stderr.write('\r'+bam+':'+str(bamn+1)+'/'+str(len(bam_list))+': loading bam file...')
		sys.stderr.write("\033[K")
		samfile = pysam.AlignmentFile(sam_path, "rb")

		# Coverage 계산
		sys.stderr.write('\r'+bam+':'+str(bamn+1)+'/'+str(len(bam_list))+': coverage calculating...')
		sys.stderr.write("\033[K")
		assdfe = samfile.count_coverage(contig, start=start, stop=stop+1)

		sys.stderr.write('\r'+bam+':'+str(bamn+1)+'/'+str(len(bam_list))+': converting coverage to array...')
		sys.stderr.write("\033[K")
		cv_original = np.array(assdfe)
		samfile.close()

		# Coverage가  A, T, G, C 4성분으로 다 따로 출력되기 때문에,
		# 이를 합쳐주는 작업을 한다.
		sys.stderr.write('\r'+bam+':'+str(bamn+1)+'/'+str(len(bam_list))+': coverage adding...')
		sys.stderr.write("\033[K")
		cv = cv_original.sum(axis=0) 

		# caching
		sys.stderr.write('\r'+bam+':'+str(bamn+1)+'/'+str(len(bam_list))+': caching...')
		sys.stderr.write("\033[K")
		cv_np = np.array(cv)
		np.save(cache_path, cv_np)
	else :
		sys.stderr.write('\r'+bam+':'+str(bamn+1)+'/'+str(len(bam_list))+': load cache...')
		sys.stderr.write("\033[K")
		cv = list(np.load(cache_path+'.npy'))
	

	# coverage의 최댓값으로 나머지 값을 나누어주는 작업을 한다.
	# 그러면 전체 값의 범위가 0~1이 된다.
	# 전체 coverage의 최댓값
	maxcv = max(cv)	
	# 모든 위치에 대해 최댓값으로 나누어 준다.
	for j, out in enumerate(cv) :
		coverage[j].append(out/maxcv)
#		coverage[j].append(out)



print()	



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
if combine :
	for k, es in enumerate(ess_nm) :
		# refseq의 start와 비교해 그보다 작으면 start를 시작점으로 한다.
		draw_range.append(es)
		draw_range_e.append(ees_nm[k])
	draw_range[0] = start if draw_range[0]-100 < start else draw_range[0]-100
	draw_range_e[-1] = stop+1 if draw_range_e[-1]+100 > stop+1 else draw_range_e[-1]+100
elif flag :	
	for k, es in enumerate(ess_nm) :
		# refseq의 start와 비교해 그보다 작으면 start를 시작점으로 한다.
		if es-100 >= start :
			draw_range.append(es-100)
		else :
			draw_range.append(start)
		draw_range_e.append(ees_nm[k]+100)
else :
	draw_range = range(start, stop+1, draw_span)


# 위 배열에 해당하는 범위를 그래프로 그린다.
# for문을 시작점으로 돌린다.
if combine :

	cv_whole = []
	
	for n, st in enumerate(draw_range) :
		stop_n = draw_range_e[n]
		cv_whole += coverage[st-start:stop_n-start]

	max_whole = max(list(map(max, cv_whole)))

	'''
	Clustering Test
	'''
	
	df = pd.DataFrame(cv_whole)

	df_T = df.T
	data_lines = df_T.values
	kmeans = KMeans(n_clusters=2).fit(data_lines)

	df_T['cluster_id'] = kmeans.labels_
	df_c1 = df_T[df_T.cluster_id == 0].drop('cluster_id', 1).T
	df_c2 = df_T[df_T.cluster_id == 1].drop('cluster_id', 1).T	

	df_c1_mean = np.mean(df_c1.values)
	df_c2_mean = np.mean(df_c2.values)
	
	c1c = ''
	c2c = ''


	if df_c1_mean < df_c2_mean :
		c1c = 'g'
		c2c = 'b'
	else :
		c1c = 'b'
		c2c = 'g'




	drops = []
	rises = []
	boths = []
	drop_means = []
	rise_means = []

	if view_mode == 2 :

#		df_drop = pd.DataFrame(cv_whole)
		xticks = np.arange(start, stop+1)
		df_drop = pd.DataFrame(coverage, index=xticks)

		ci = 1.96
		
		drop_n = select_list[0]-1
		rise_n = select_list[1]-1
		dintv = 30



		for cn, c in enumerate(df_drop.columns) :

#			e1_e = draw_range[0] - start
#			for rn, range_e in enumerate(draw_range_e) :
#				if rn > drop_n:
#					break
#				range_s = draw_range[rn]
#				e1_e += range_e - range_s + 1

			e1_e = draw_range_e[drop_n] - start 
			e2_s = draw_range[drop_n+1] - start 

			df_cl = df_drop.iloc[e1_e-dintv:e1_e+1][c]

			cl_means = df_cl.mean()
			cl_std = df_cl.std()

			cl_high = cl_means + ci*cl_std
			cl_low = cl_means - ci*cl_std

			df_cr = df_drop.iloc[e2_s:e2_s+dintv+1][c]
			cr_means = df_cr.mean()


			if cr_means < threshold :
				if cl_high < cr_means or cl_low > cr_means :
					drops.append(cn)
					drop_means.append((cn, df_cl, df_cr))		


		for cn, c in enumerate(df_drop.columns) :

#			e1_e = draw_range[0] - start
#			for rn, range_e in enumerate(draw_range_e) :
#				if rn > rise_n:
#					break
#				range_s = draw_range[rn]
#				e1_e += range_e - range_s + 1

			e1_e = draw_range_e[rise_n] - start
			e2_s = draw_range[rise_n+1] - start


			df_cl = df_drop.iloc[e1_e-dintv:e1_e+1][c]
			cl_means = df_cl.mean()

			df_cr = df_drop.iloc[e2_s:e2_s+dintv+1][c]
			cr_means = df_cr.mean()
			cr_std = df_cr.std()

			cr_high = cr_means + ci*cr_std
			cr_low = cr_means - ci*cr_std

			if cr_high < cl_means or cr_low > cl_means :
				rises.append(cn)
				rise_means.append((cn, df_cl, df_cr))		

		print('drops : '+str(len(drops)))
		print('rises : '+str(len(rises)))
		
#		drop_means = reversed(sorted(drop_means, key=lambda means: means[1]))
#		print(drop_means)
#		for dn, dm in enumerate(drop_means) :
#			if dn > 3 :
#				break
#			if dm[0] in drops :
#				drops.remove(dm[0])
		for dn in drops :
			if dn in rises :
				print(bam_list[dn])
				boths.append(dn)

		print('boths : '+str(len(boths)))

		

#		for dn, dm in enumerate(drop_means) :
#			if dm[0] in boths :
#				print(dm[1], dm[2])
#				print('----------------------------------------------------------------')




	refseq_r = refseq[refseq.txStart <= draw_range_e[-1]]
	refseq_r = refseq_r[refseq_r.txEnd >= draw_range[0]]
	
	nmids_whole = refseq_r.name.tolist()

	nl_whole = len(nmids_whole)

	# figure 설정
	# 여러 칸으로 나눈다.
	fig = plt.figure(figsize=(120, 2+12+1*nl_whole))
	gs = gridspec.GridSpec(nrows=2, ncols=len(ess_nm), height_ratios=[4+(nl_whole-1)*0.2, nl_whole])
	gs.update(wspace=0, hspace=0.05)

	# 첫번째 칸부터 차례로 채워나간다.
	for n, st in enumerate(draw_range) :
		sys.stderr.write('\r'+'exon '+str(n+1)+' drawing...')
		sys.stderr.write("\033[K")

		stop_n = draw_range_e[n]
	
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



		# coverage의 dataframe을 만들고 plot의 윗 부분을 불러온다.
		df2 = pd.DataFrame(coverage[st-start:stop_n-start], index=xticks)
		ax_main = plt.subplot(gs[n])
		if not n == 0 :
			ax_main.get_yaxis().set_ticks([])
			ax_main.spines['left'].set_visible(False)
		if not n == len(draw_range)-1 :
			ax_main.spines['right'].set_visible(False)

		'''
		Clustering Test
		'''
		
		df2_T = df2.T
		data_lines = df2_T.values

		df2_T['cluster_id'] = kmeans.labels_
		df2_c1 = df2_T[df2_T.cluster_id == 0].drop('cluster_id', 1).T
		df2_c2 = df2_T[df2_T.cluster_id == 1].drop('cluster_id', 1).T	

		df2_c1_mean = np.mean(df2_c1.values)
		df2_c2_mean = np.mean(df2_c2.values)

		labels = []


	#	colors = ['g' if lb == 1 else 'b' for lb in labels]


		# c1
		ax_main.plot(df2_c1, color=c1c, alpha=0.1)
	#	ax_main.fill_between(xticks, df2_c1.max(axis=1), df2_c1.min(axis=1), facecolor=c1c, alpha=0.5)
	#	ax_main.plot(df2_c1.max(axis=1), color=c1c, linewidth=3.0)
	#	ax_main.plot(df2_c1.min(axis=1), color=c1c, linewidth=3.0)
		if view_mode == 1 :
			ax_main.plot(df2_c1.mean(axis=1), color='red', linewidth=3.0)	

		# c2
		ax_main.plot(df2_c2, color=c2c, alpha=0.1)
	#	ax_main.fill_between(xticks, df2_c2.max(axis=1), df2_c2.min(axis=1), facecolor=c2c, alpha=0.5)
	#	ax_main.plot(df2_c2.max(axis=1), color=c2c, linewidth=3.0)
	#	ax_main.plot(df2_c2.min(axis=1), color=c2c, linewidth=3.0)
		if view_mode == 1 :
			ax_main.plot(df2_c2.mean(axis=1), color='red', linewidth=3.0)	
		ax_main.set_ylim([0, max_whole])
		ax_main.set_xlim([st, stop_n])

		if view_mode == 2 :
#			ax_main.plot(df2.iloc[:, drops], color='red', linewidth=0.5)	
#			ax_main.plot(df2.iloc[:, rises], color='blue', linewidth=0.5)	
			if len(boths) > 0:
				ax_main.plot(df2.iloc[:, boths], color='red', linewidth=1.0)	

	#	plt.ylim(0, 1) # coverage 표시 범위 설정
	#	ax_main.set_yscale('log') # logscale 여부
		plt.xticks(np.arange(st, stop_n+1, step=(stop_n-st)/2)) # x축 값 표시.

		# bp 위치의 값이 너무 크기 때문에 e를 포함한 값으로 표시된다.
		# 따라서 아래와 같은 과정을 거쳐 자연수가 나오도록 한다.
		xx, locs = plt.xticks()
		ll = ['%d' % a for a in xx]
		plt.xticks(xx, ll)

		# 1 기준선을 그린다.
		# reddot = np.ones(stop_n-st)
		# ax_main.plot(xticks, reddot, 'r--')

		# refseq를 표시할 범위를 정한다.
		byts = range(2-nl, 2)

		# refseq가 표시될 이름을 정한다.
		ytlbs = [aa+"\n"+names[aai] for aai, aa in enumerate(nmids)]
		ax_bottom = plt.subplot(gs[n+len(ess_nm)], yticks=byts, xticklabels=[], yticklabels=list(reversed(ytlbs)))
		ax_bottom.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False)
		
		if not n == 0 :
			ax_bottom.get_yaxis().set_ticks([])
			ax_bottom.spines['left'].set_visible(False)
		if not n == len(draw_range)-1 :
			ax_bottom.spines['right'].set_visible(False)
		ax_bottom.set_xlim([st, stop_n])
		

		# refseq의 범위를 검은 선으로 그린다.
		for j, ts in enumerate(tx_s) :	
			a_s = ts if ts > st else st
			a_e = tx_e[j] if tx_e[j] < stop_n else stop_n
			if a_e-a_s < 0 :
				continue
			bxts = np.arange(a_s, a_e)
			blns = np.full(a_e-a_s,  1-j, dtype=int)
			ax_bottom.plot(bxts, blns, 'black')

		# tx, cds의 시작 부분을 얇은 네모로 그린다.
		for j, cs in enumerate(cds_s) :	
				if (tx_s[j] > stop_n or cs < st) :
					continue
				rect = patches.Rectangle((tx_s[j], 0.9-j),cs-tx_s[j],0.2,edgecolor='none',facecolor='black')
				ax_bottom.add_patch(rect)

		# tx, cds의 끝 부분을 얇은 네모로 그린다.
		for j, ce in enumerate(cds_e) :
				if (ce > stop_n or tx_e[j] < st) :
					continue
				rect = patches.Rectangle((ce, 0.9-j),tx_e[j]-ce,0.2,edgecolor='none',facecolor='black')
				ax_bottom.add_patch(rect)

		# 방향을 그린다.
#		print('draw directions...')
#		for j, ts in enumerate(tx_s) :
#			strand = strands[j]
#			if (ts > stop_n or tx_e[j] < st) :
#				continue
#			interval = int((stop_n-st)/3)
#			a_s = ts if ts > st else st
#			a_e = tx_e[j] if tx_e[j] < stop_n else stop_n 
#			for k in range(a_s, a_e, interval) :
#				if strand == '+' :
#					ax_bottom.arrow(k, 1-j, interval, 0, head_width=0.2, head_length=interval/2, overhang=1)
#				else :	
#					ax_bottom.arrow(k, 1-j, interval*(-1), 0, head_width=0.2, head_length=interval/2, overhang=1)

		# 엑손을 그린다.
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

	print()
	# lineplot 간 여백 설정
	plt.subplots_adjust(top = 0.95, bottom = 0.05, right = 0.97, left = 0.03, wspace=0, hspace=0)
	# 폰트 사이즈 조정
	matplotlib.rcParams.update({'font.size': 22})

	# pdf로 저장한 뒤 닫는다.
	plt.savefig(output_prefix+'.pdf')
	plt.close(fig)
	print(output_prefix+'.pdf saved!')





else : 
	for n, st in enumerate(draw_range) :

		print('\n'+output_prefix+'_'+str(n+1)+' saving...')

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

		df2_c1_mean = np.mean(df2_c1.values)
		df2_c2_mean = np.mean(df2_c2.values)

		labels = []

		c1c = ''
		c2c = ''

		if df2_c1_mean < df2_c2_mean :
			c1c = 'g'
			c2c = 'b'
		else :
			c1c = 'b'
			c2c = 'g'

	#	colors = ['g' if lb == 1 else 'b' for lb in labels]

		# c1
		ax_main.plot(df2_c1, color=c1c, alpha=0.1)
	#	ax_main.fill_between(xticks, df2_c1.max(axis=1), df2_c1.min(axis=1), facecolor=c1c, alpha=0.5)
	#	ax_main.plot(df2_c1.max(axis=1), color=c1c, linewidth=3.0)
	#	ax_main.plot(df2_c1.min(axis=1), color=c1c, linewidth=3.0)
		ax_main.plot(df2_c1.mean(axis=1), color='red', linewidth=3.0)	

		# c2
		ax_main.plot(df2_c2, color=c2c, alpha=0.1)
	#	ax_main.fill_between(xticks, df2_c2.max(axis=1), df2_c2.min(axis=1), facecolor=c2c, alpha=0.5)
	#	ax_main.plot(df2_c2.max(axis=1), color=c2c, linewidth=3.0)
	#	ax_main.plot(df2_c2.min(axis=1), color=c2c, linewidth=3.0)
		ax_main.plot(df2_c2.mean(axis=1), color='red', linewidth=3.0)	
		

	#	plt.ylim(0, 1) # coverage 표시 범위 설정
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




