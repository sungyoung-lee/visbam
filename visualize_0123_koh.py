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
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import f1_score, balanced_accuracy_score
from sklearn.metrics import silhouette_score, silhouette_samples
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from skmisc import loess
from sklearn.decomposition import NMF





## Argument Setting ##


# 파일 이름과 span을 argument로 불러들인다.
parser = argparse.ArgumentParser()
# required argument
parser.add_argument('--bam_path', help='bam파일을 읽어들일 디렉토리를 정합니다.')
parser.add_argument('--sample_path', help='해당하는 sample 이름이 들어있는 경로를 지정합니다.')
parser.add_argument('--normal_path', help='normal sample이 들어있는 경로를 지정합니다.')
parser.add_argument('--refseq_path', help='refseq 파일의 경로를 지정합니다.')
parser.add_argument('--variants_path', help='Generic Variants 파일 경로를 지정합니다.')
parser.add_argument('--refseq', help='사용할 NMID를 지정합니다.')
parser.add_argument('--prefix', help='output 파일명을 정합니다.')
# optional argument
# draw plot
parser.add_argument('--exon_sliced', help='exon 주변의 100bp 부분만 그립니다.', action='store_true')
parser.add_argument('--curated_genes', help='curated genes')
parser.add_argument('--exclude_exon', default='', help='주어진 exon을 제외하고 계산&출력합니다.(1,2,3,4,...)')
parser.add_argument('--combine_slices', help='모든 분리된 그림을 붙여 출력합니다.', action='store_true')
parser.add_argument('--draw_average_line', action='store_true', help='그래프에 평균값을 빨간 선으로 표시합니다.')
parser.add_argument('--draw_span', type=int, default=10000, help='사진을 몇 bp단위로 분할할 것인지 정합니다.')
parser.add_argument('--smoothing', default='average', help='그래프를 smoothing할 모드를 지정합니다. (average, loess)')
parser.add_argument('--average', type=int, default=0, help='주어진 정수값 span만큼 Moving Average를 적용합니다.')
parser.add_argument('--fill', help='moving average과정에서 바깥의 값을 가져올지 경계의 값으로 채울지 결정합니다.', action='store_true')
parser.add_argument('--font_size', type=int, default=7, help='본 그래프의 fontsize를 조정합니다. (단위:pt)')
parser.add_argument('--marker_size', type=int, default=9, help='variant marker의 size를 조정합니다. (단위:pt)')
parser.add_argument('--ylim', type=int, default=0, help='표시할 y축의 최댓값을 정합니다.')
parser.add_argument('--exon_space', type=int, default=0, help='exon_sliced일 때 표시할 exon 주위 간격을 설정합니다.')
parser.add_argument('--min_max', action='store_true', help='그래프의 최댓값과 최솟값만 표시합니다.')

parser.add_argument('--exclude_sample', default='', help='출력에서 제외할 sample을 입력합니다.')

# clustering
parser.add_argument('--clustering', action='store_true', help='주어진 그래프를 두 그룹으로 clustering합니다.')
parser.add_argument('--clustering_mode', default='silhouette', help='view_mode 2에서 filtering할 method 설정(silhouette, nmf, splice_site)')
parser.add_argument('--select_exon', default='', help='clustering에서 select할 exon 정의')

# silhouette clustering 
parser.add_argument('--score_plot_width', type=int, default=12, help='clustering에서score scatter plot의 width를 조정합니다. (단위:inch)')
parser.add_argument('--score_plot_height', type=int, default=12, help='clustering에서score scatter plot의 height를 조정합니다. (단위:inch)')
parser.add_argument('--limit_tau', type=float, default=None, help='clustering에서 tau limit 값 조정')
parser.add_argument('--limit_tau_low', type=float, default=None, help='clustering에서 low tau limit 값 조정')
parser.add_argument('--silhouette_dintv', type=int, default=30, help='clustering silhouette 모드에서 exon 계산 범위 (단위:bp)')

parser.add_argument('--select_CI', default='', help='clustering silhouette 모드에서 신뢰구간(CI)을 선택합니다.')
parser.add_argument('--select_tau', default='', help='silouette clustering일 때 tau값을 고정합니다.')
parser.add_argument('--variant_pos', type=int, default=None, help='clustering mode가 silhouette일때, clustering ratio를 확인할 variant의위치를 정합니다.')

# scatter clustering
parser.add_argument('--input_space',type=int, default=None, help='scatter plot을 그릴 때 계산하는 간격')
parser.add_argument('--x_division', type=float, default=None, help='scatter로 clustering을 할 때 x 기준점을 정합니다.')
parser.add_argument('--y_division', type=float, default=None, help='scatter로 clustering을 할 때 y 기준점을 정합니다.')


args = parser.parse_args()
bam_dir = args.bam_path
sample_list_path = args.sample_path
normal_dir_path = args.normal_path
refseq_path = args.refseq_path
nmid_to_draw = args.refseq
variants_dir = args.variants_path
draw_span = args.draw_span
output_prefix = args.prefix
flag = args.exon_sliced
flag3 = args.curated_genes
exclude = args.exclude_exon.strip()
combine = args.combine_slices
max_whole = args.ylim
min_max = args.min_max
exon_space = args.exon_space
view_mode = args.clustering
filt_mode = args.clustering_mode
draw_average_line = args.draw_average_line
select = args.select_exon
average = args.average
smooth = args.smoothing
fill = args.fill
score_plot_width = args.score_plot_width
score_plot_height = args.score_plot_height
limit_tau = args.limit_tau
limit_tau_low = args.limit_tau_low
font_size = args.font_size
marker_size = args.marker_size
silhouette_dintv = args.silhouette_dintv

input_space = args.input_space
select_tau = args.select_tau
select_CI = args.select_CI
variant_pos = args.variant_pos
x_div = args.x_division
y_div = args.y_division
exclude_sample=args.exclude_sample

# set title

title = 'graph'



# string argument를 배열로 변환
exc_samples = [] if exclude_sample == '' else list(exclude_sample.split(','))
exclude_list = [] if exclude == '' else reversed(sorted(list(map(int, exclude.split(',')))))
select_list = [1, 7] if select == '' else list(map(int, select.split(',')))
CI_list = [] if select_CI == '' else list(map(float,select_CI.split(',')))
Tau_list = [] if select_tau == '' else list(map(float,select_tau.split(',')))



# curated된 refseq만 분류하는 함수

def include(refseq) :
	# curated된 refseq 목록을 불러온다.
	curated_f = open(flag3, 'r')
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





## Reading Refseq Data ##


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
strand_nm = refseq_nm.strand.tolist()
tx_s_nm = refseq_nm.txStart.tolist()
tx_e_nm = refseq_nm.txEnd.tolist()
cds_s_nm = refseq_nm.cdsStart.tolist()
cds_e_nm = refseq_nm.cdsEnd.tolist()
exon_s_nm = refseq_nm.exonStarts.tolist()
exon_e_nm = refseq_nm.exonEnds.tolist()
name2_nm = refseq_nm.name2.tolist()

# 그 유전자를 토대로 그릴 범위를 정한다.
name_n = name_nm[0]
contig = chrom_nm[0]
start = tx_s_nm[0]
stop = tx_e_nm[0]
ess_nm = list(map(int, exon_s_nm[0].split(',')[:-1]))
ees_nm = list(map(int, exon_e_nm[0].split(',')[:-1]))

Exon_s_list = ess_nm[:]
Exon_e_list = ees_nm[:]


# 그릴 refseq 정보만 filtering 후 불러온다
refseq = refseq[refseq.name.str.contains("NM")]
refseq = refseq[refseq.txStart <= stop]
refseq = refseq[refseq.txEnd >= start]
if not flag3 == None :
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





## Bam Information Analysis ##

# cache folder 있는지 확인
if not os.path.isdir('cache'):
	os.mkdir('cache')


# coverage를 저장할 데이터를 초기화한다.
coverage = [[] for i in range(stop-start+1)]
normal_coverage = np.zeros(stop-start+1)
samfile = None


## Normal Bam Coverage ##
print('Analyzing normal bam information...')

# Normal Bam의 파일 리스트
Nbam_list = os.listdir(normal_dir_path)
Nbam_list = [file for file in Nbam_list if file.endswith(".bam")]


for Nbamn, Nbam in enumerate(Nbam_list) :

	print('\r', Nbam, end='')
	sys.stdout.write("\033[K")

	# Normal Bam 경로
	sam_path = normal_dir_path+'/'+Nbam
	
	# Normal Bam이 파일인지 확인
	if not os.path.isfile(sam_path) :
		continue

	nm_fn = normal_dir_path[normal_dir_path.rfind('/')+1:]

	# cache path
	cache_path = 'cache/'+nm_fn+'_'+Nbam+'_'+name_n+'_'+contig+'_'+str(start)+'_'+str(stop)
	print(cache_path)	

	cv = []

	# 해당 캐시가 없을 경우 coverage 계산
	if not os.path.isfile(cache_path+'.npy') :
		# Bam파일을 불러온다.
		samfile = pysam.AlignmentFile(sam_path)
		print('loaded')

		# Coverage 계산
		cv_original = np.array(samfile.count_coverage(contig, start=start, stop=stop+1))
		print('calculated')
		samfile.close()
		# Coverage가  A, T, G, C 으로 다 따로 출력되기 때문에,
		# 이를 합쳐주는 작업을 한다.
		cv = cv_original.sum(axis=0) 

		# caching
		cv_np = np.array(cv)
		np.save(cache_path, cv_np)
	# 캐시가 있으면 캐시를 불러온다.
	else :
		cv = list(np.load(cache_path+'.npy'))

	normal_coverage = list(map(add, normal_coverage, cv))

normal_coverage = [x / len(Nbam_list) for x in normal_coverage]



## Cancer Bam Coverage ##

print('\nAnalyzing cancer bam information...')

# sample list가 들어있는 파일을 불러온다.
sample_file = open(sample_list_path, 'r')
sample_list = sample_file.read().splitlines()

if not exc_samples == '':
	for n, s in enumerate(sample_list):
		if s in exc_samples:
			sample_list.remove(s)


bam_list = []
sample_file.close()
real_bam_list = []

# 목록을 불러온 뒤 유효한 파일 이름으로 바꾸어준다.
for sl in sample_list:
	bam_list.append(sl+'.bwamem.sorted.dedup.realn.recal.dedup.bam')



for bamn, bam in enumerate(bam_list) :

	sys.stderr.write('\r'+bam+':'+str(bamn+1)+'/'+str(len(bam_list))+': start...')
	sys.stderr.write("\033[K")
	# Cancer Bam 경로
	sam_path = bam_dir+'/'+bam	

	# Cancer Bam이 파일인지(존재하는지)  확인
	if not os.path.isfile(sam_path) :
		continue

	# cache path
	real_bam_list.append(bam)
	bam_fn = os.path.basename(bam)
	cache_path = 'cache/'+bam_fn+'_'+name_n+'_'+contig+'_'+str(start)+'_'+str(stop)
	
	cv = []

	# 해당 캐시가 없을 경우 coverage 계산
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

		# Coverage가  A, T, G, C 으로 다 따로 출력되기 때문에,
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
		# 캐시가 있으면 캐시를 불러온다.
		sys.stderr.write('\r'+bam+':'+str(bamn+1)+'/'+str(len(bam_list))+': load cache...')
		sys.stderr.write("\033[K")
		cv = list(np.load(cache_path+'.npy'))
	
	
	# Normal Coverage 평균으로 Normalize
	for j, out in enumerate(cv) :
		# normal_coverage[j]가 0이면 1로 처리한다.
		cov = 1
		if normal_coverage[j] != 0 :
			cov = out/normal_coverage[j]
		coverage[j].append(cov)	

print()



## Special Sample Coverage ##

print('Loading special sample...')

coverage_special = [[] for i in range(stop-start+1)]

# special sample의 파일 이름
bam = 'MS190000674_S13.bwamem.sorted.dedup.realn.recal.dedup.bam'

# special sample의 경로
sam_path = '200116_work/'+bam

# cache path
bam_fn = os.path.basename(bam)
cache_path = 'cache/'+bam_fn+'_'+name_n+'_'+contig+'_'+str(start)+'_'+str(stop)

cv = []

# Cache가 있는지 없으면 coverage 계산
if not os.path.isfile(cache_path+'.npy') :
	# Bam파일을 불러온다.
	samfile = pysam.AlignmentFile(sam_path, "rb")

	# Coverage 계산
	assdfe = samfile.count_coverage(contig, start=start, stop=stop+1)

	cv_original = np.array(assdfe)
	samfile.close()

	# Coverage가  A, T, G, C 4성분으로 다 따로 출력되기 때문에,
	# 이를 합쳐주는 작업을 한다.
	cv = cv_original.sum(axis=0) 

	# caching
	cv_np = np.array(cv)
	np.save(cache_path, cv_np)
# Cache 불러오기
else :
	cv = list(np.load(cache_path+'.npy'))

# special coverage를 normal sample로 normalize
for j, out in enumerate(cv) :
	# normal_coverage[j]가 0이면 1로 처리한다.
	cov = 1
	if normal_coverage[j] != 0 :
		cov = out/normal_coverage[j]
	coverage_special[j].append(cov)	





## Load Genetic Variants Data ##

bam_all_list = os.listdir(bam_dir)
bam_all_list = [file for file in bam_list if file.endswith(".bam")]


# variants를 모두 불러온다.
var_list = os.listdir(variants_dir)
var_files = [file for file in var_list if file.endswith(".txt")]

# variants 데이터를 저장할 dataframe 생성
df_var = pd.DataFrame(columns=['index','MSID','pos','effect'])

# 실제 불러온 bam list의 sample 이름을 추출
bam_names = [b[:b.find('.')] for b in real_bam_list]

# 마커 모양 지정
var_markers = ['D', '*', 'v', 's', 'P', 'h', 'x']

# 마커 색 지정
var_colors = ['blue', 'yellow', 'magenta', 'black']

# Variants 데이터 수집
for varf in var_files :
	sys.stderr.write('\r'+varf+' variants loading.')
	sys.stderr.write("\033[K")


	# 파일 csv로 불러오기
	df_var_list = pd.DataFrame.from_csv(variants_dir+'/'+varf, sep='\t')
	varf_name = varf[:varf.find('.')] # var file명에서 sample명 추출 

	# 그러나 해당 varinats의 sample이 앞서 불러온 sample에 없으면 넘김
	if varf_name in bam_names :
		bam_index = bam_names.index(varf_name)
	else :
		continue
	
	# 데이터 저장
	for i, var in df_var_list.iterrows() :
		# Refseq 정보 읽어들이기
		var_refseq = var['Refseq']
		
		if var_refseq == var_refseq :
			var_refseq = var_refseq[:var_refseq.find('.')]

			# 선택한 NMID의 variants만 표시
			if var_refseq.strip() == nmid_to_draw.strip() :
				# Effect 앞부분만 표시
				var_effect = var['Effect']
				if not var_effect.find('+') == -1 :
					var_effect = var_effect[:var_effect.find('+')]
				# sample 번호, position, effect 정보를 한 row에 저장
				df_var = df_var.append({'index':bam_index,'MSID':varf_name,'pos' : var['Pos'], 'effect' : var_effect}, ignore_index=True)


#pos별 variant를 가진sample 개수 저장하는 dataframe 저장
pos_var = pd.DataFrame(columns=['pos','cnt'])

for pos in range(start,stop+1):
	cnt = 0
	for var_pos in df_var['pos']:
		if pos == var_pos:
			cnt=cnt+1		
	if cnt == 0:
		continue
	
	pos_var = pos_var.append({'pos':pos, 'cnt':cnt},ignore_index=True)

 
# Generic Variants의 Effect 종류를 모두 모아놓은 배열
effect_list = sorted(list(set(np.array(df_var['effect'].tolist()).squeeze())))

# Generic Variants가 있는 bam을 모두 모아놓은 배열
bam_num_list = list(set(np.array(df_var['index'].tolist()).squeeze()))
bam_var_list = set(df_var['MSID'].tolist())

print('\n----------------------------------------------')
for b in bam_num_list :
	print(real_bam_list[int(b)])
print('----------------------------------------------')





## Draw Range ##

# Line plot을 그릴 위치를 정합니다.
# 각 그림의 시작
# 별다른 옵션이 없으면 txStart부터 그림을 그릴거야
draw_start = start
draw_range = []

# 각 그림의 끝
# 별다른 옵션이 없으면 txEnd까지 그림을 그릴거야
draw_stop = stop
draw_range_e = []

# exclude_exon 옵션이 켜져있을 때 draw할 exon을 조정해주자
exons = list(np.arange(1,len(ess_nm)+1))

if not exclude == '':
	draw_stop = 0
	for i in exclude_list:
		if i <= len(ess_nm):
			del ess_nm[i-1]
			del ees_nm[i-1]
			exons.remove(i)

	draw_stop = ees_nm[-1]


# flag 속성이 켜져 있으면, 지정한 refseq의 각 exon 주위의 100bp 만큼의 범위를 그린다.
# 그렇지 않으면, refseq 전체 범위를 지정한 크기만큼 잘라 그린다.
e_bp = exon_space

if flag : 
	#그림그리는 가장 시작과 끝을 정해줍니다.
	for i, ess in enumerate(ess_nm) :
		# draw_range 정하기
		# refseq의 start와 비교해 그보다 작으면 start를 시작점으로 한다.
		if i>0 and ess-e_bp < ees_nm[i-1] :
			draw_range.append(ees_nm[i-1]+1)
		else :
			draw_range.append(ess-e_bp)

		# draw_range_e 정하기
		if i<len(ees_nm)-1 and ees_nm[i] + e_bp > ess_nm[i+1] :
			draw_range_e.append(ess_nm[i+1]-1)
		else :
			draw_range_e.append(ees_nm[i]+e_bp)

	draw_range[0] = draw_start if draw_range[0] <= draw_start else ess_nm[0]
	draw_range_e[-1] = draw_stop if draw_range_e[-1] < draw_stop else draw_range_e[-1]

else :
	draw_range = list(range(draw_start, draw_stop+1, draw_span))
	draw_range_e = list(range(draw_start+draw_span, draw_stop+1, draw_span))

	if draw_stop+1 in draw_range :
		draw_range.remove(draw_stop+1)

	if not draw_stop+1 in draw_range_e :
		draw_range_e.append(draw_stop+1)


# 위 배열에 해당하는 범위를 그래프로 그린다.
# for문을 시작점으로 돌린다.

cv_whole = []

for n, st in enumerate(draw_range) :
	stop_n = draw_range_e[n]
	cv_whole += coverage[(st-start):(stop_n-start)]





## Smoothing(Simplification) ##

# df에 해당하는 xticks 설정
xticks_whole = np.arange(start, start+len(coverage))

# smoothing에 사용할 dataframe 생성
df2_whole = pd.DataFrame(coverage, index=xticks_whole)
df2_mean = None


if smooth == 'average' and average > 1 :
	# average 에서 exon 앞 뒤 간격을 주고 smoothing 할지 말지 결정
	if fill :
		coverage_df2 = df2_whole.values.tolist()
		# exon별로 따로 처리한다
		for n, st in enumerate(draw_range) :
			sys.stderr.write('\r'+'Exon '+str(n+1)+' smoothing...')
			sys.stderr.write("\033[K")
			stop_n = draw_range_e[n]-1
			exon_list = [coverage_df2[st-start]]*average
			exon_list += coverage_df2[(st-start):(stop_n-start)]
			exon_list += [coverage_df2[stop_n-start]]*average
			df2_exon = pd.DataFrame(exon_list)
			df2_exon = df2_exon.rolling(window=average*2, min_periods=1, center=True).mean()
			exon_list = df2_exon.values.tolist()
			# smoothing이 완료된 exon을 그 부분만 잘라 배열에 저장
			coverage_df2[(st-start):(stop_n-start)] = exon_list[average:(average+(stop_n-st))]

		df2_mean = pd.DataFrame(coverage_df2, index=xticks_whole)
	else :
		# 그래프 전체에 대해 rolling 실시
		df2_mean = df2_whole.rolling(window=average*2, min_periods=1, center=True).mean()

elif smooth == 'loess' :
	# 처리 시간상 loess를 전체 그래프에 대해 실행 할 수 없어
	# sample별로, exon별로 따로 처리한다.
	loess_list = []
	for di in range(len(df2_whole.columns)) :
		# 한 샘플의 dataframe
		loess_di = df2_whole.iloc[:, di].tolist()
		for n, st in enumerate(draw_range) :
			sys.stderr.write('\r'+'DNA '+str(di+1)+' / '+str(len(df2_whole.columns))+', Exon '+str(n+1)+' smoothing...')
			sys.stderr.write("\033[K")
			stop_n = draw_range_e[n]
			loess_xticks = xticks_whole[(st-start):(stop_n-start)]
			# 한 샘플, 한 exon에 대한 배열로 loess 실행
			loess_value = loess.loess(loess_xticks, df2_whole.iloc[(st-start):(stop_n-start), di].tolist())
			loess_value.fit()
			loess_predict =	loess_value.predict(loess_xticks, stderror=True).values
			# 환성된 결과를 원래 배열에 저장
			loess_di[(st-start):(stop_n-start)] = loess_predict
		# 별도에 배열에 한 sample에 대한 결과 저장
		loess_list.append(loess_di)
	df2_mean = pd.DataFrame(loess_list)
	df2_mean = df2_mean.T
	df2_mean = df2_mean.set_index(xticks_whole)
else : 
	# 처리 과정이 없거나 이름이 잘못된 경우 smoothing 처리 하지 않음
	df2_mean = df2_whole


xticks = np.arange(start, stop+1)

cv_mean = pd.DataFrame(df2_mean.mean(axis=1),index=xticks)
df = pd.DataFrame(cv_whole)
df_special = pd.DataFrame(coverage_special, index=xticks) # Special Plot의 dataframe

# 표시된 부분을 Dataframe으로 만든다.
df_T = df.T
data_lines = df_T.values

print(df2_mean.shape)





## Clustering ##
# clustering 최종 결과 저장
drops = []
rises = []
boths = []
drop_means = []
rise_means = []
boths_01 = np.zeros(len(coverage[0]))
CI = []

if view_mode :

	# 최종 optimal 결과 저장
	highest_score = -1	
	highest_ci = 0
	highest_tau = 0
	highest_tnum = 0
	highest_ratio = 0

	# 최종optimal 결과 저장
	ci_list = []
	tau_list = []
	score_list = []
	ratio_list = []
	ratio_list_t = []
	green_list = []
	red_list = []
	variant_list=[]
	var_index = []


	# 선택된 exon number 불러오기
	# 배열에 들어갈 값이므로 1을 빼줍니다. (0부터 시작)
	drop_n = select_list[0]-1
	rise_n = select_list[1]-1
	# exon 계산 범위 앞뒤 간격
	dintv = silhouette_dintv

	# 계산을 진행할 dataframe
	xticks = np.arange(start, stop+1)
	df_drop = pd.DataFrame(coverage, index=xticks)

	e1_l = Exon_e_list[drop_n]-start
	e1_r = Exon_s_list[drop_n+1]-start
	
	e2_l = Exon_e_list[rise_n]-start
	e2_r = Exon_s_list[rise_n+1]-start

	

	# 1번째 방법 silhouette
	if filt_mode == 'silhouette' : 
		
		if not select_CI == '' :
			CI = CI_list			
		else:
			CI=[0.5,1,1.5]+list(np.arange(2,26,0.5))	
		print("Confidence interval is",CI)


		for ci in CI:
			
			# 임시로 저장할 배열
			drops_t = []
			rises_t = []
			boths_t = []
			drop_means_t = []
			rise_means_t = []	
		
			# Drops 계산
			for cn, c in enumerate(df_drop.columns) :

				# 1번째 exon 끝에서 dintv 만큼 간격의 값 계산(left) 
				e1_cl = df_drop.iloc[(e1_l-dintv):(e1_l+1)][c]

				# 평균, 표준편차 계산 후 오차범위 계산
				cl_means = e1_cl.mean()
				cl_std = e1_cl.std()
				
				cl_high = cl_means+ci*cl_std
				cl_low = cl_means-ci*cl_std
				
				# 1번째 exon+1 시작점에서 dintv 만큼 간격의 값 계산(right) 
				e1_cr = df_drop.iloc[e1_r:(e1_r+dintv+1)][c]

				# 그 간격의 평균
				cr_means = e1_cr.mean()

				# 오차범위에서 벗어나면 drops에 추가
				if cl_high < cr_means or cl_low > cr_means :
					drops_t.append(cn)
					# tau 값을 계산할 정보를 저장
					drop_means_t.append((cn, cl_means, cr_means))		

			# Rises 계산
			for cn, c in enumerate(df_drop.columns) :

				# 2번째 exon 끝에서 dintv 만큼 간격의 값 계산(left) 
				e2_cl = df_drop.iloc[(e2_l-dintv):(e2_l+1)][c]

				# 그 간격의 평균
				cl_means = e2_cl.mean()

				# w번째 exon+1 시작점에서 dintv 만큼 간격의 값 계산(right) 
				e2_cr = df_drop.iloc[e2_r:(e2_r+dintv+1)][c]

				# 평균, 표준편차 계산 후 오차범위 계산
				cr_means = e2_cr.mean()
				cr_std = e2_cr.std()

				cr_high = cr_means + ci*cr_std
				cr_low = cr_means - ci*cr_std

				# 오차범위에서 벗어나면 rises에 추가
				if cr_high < cl_means or cr_low > cl_means :
					rises_t.append(cn)
					# tau 값을 계산할 정보 저장
					rise_means_t.append((cn, cl_means, cr_means))		
			
			# 둘 중 하나라도 없으면 break
			if len(drops_t) < 1 or len(rises_t) < 1 :
				break

			# (1번째 엑손)+1 엑손의 오른쪽 값 기준으로 정렬시킨다.			
			drop_means_t = sorted(drop_means_t, key=lambda means: means[2])
			drop_means_t = list(reversed(drop_means_t))

			print(len(drops_t), len(rises_t))
	

			#tau list가 정해져 있을 때
			if not select_tau == '':
				if not(limit_tau == None or limit_tau_low):
					print("::Error:: --limit_tau, --limit_tau_low and select_tau can't be used together.")
					exit()
			
				#tau list의 크기만큼 반복
				for tn, tau in enumerate(Tau_list):
					drops_cp = drops_t.copy()
					drop_means_cp = drop_means_t.copy()	
					
					# drop과 tau 비교# drop_means_t = [[index,cl_mean,cr_mean],...]
					for dn, drop_set in enumerate(drop_means_cp):
						drop_index = drop_set[0] 
						drop_mean = drop_set[2]

						#drop의 cr_mean과 tau를 비교해서 cr_mean이 크면 drop 아님
						if drop_mean < tau:
							drops_cp.remove(drop_index)
							del drop_means_cp[dn]

					boths_t=[]
					boths_01 = np.zeros(len(coverage[0]))
					# drop과 rise에 공통적으로 있는 것만 both에 추가
					for dn in drops_cp:
						if dn in rises_t:
							boths_t.append(dn)
							boths_01[dn]=1
																
					if len(boths_t) < 1:
						continue	
										
					# silhouette score를 계산할 부분 추출
					df_silhouette = df_drop.iloc[(e1_l-dintv):(e1_l+1)] 
					df_silhouette = df_silhouette.append(df_drop.iloc[e1_r:(e1_r+dintv+1)]) 
					df_silhouette = df_silhouette.append(df_drop.iloc[(e2_l-dintv):(e2_l+1)])
					df_silhouette = df_silhouette.append(df_drop.iloc[e2_r:(e2_r+dintv+1)])
					df_silhouette = df_silhouette.T.values
					
					#각 구간의 silhoutte score 계산/ silhouette_score(sample array,label)
					score = silhouette_score(df_silhouette, boths_01)

					# 그래프를 그리기 위해 해당 결과 저장
					ci_list.append(ci)
					tau_list.append(tau)
					score_list.append(score)

					
					# 입력위치에서 cluster 그룹에 속하는 variant의 비율을 알고 싶을 때
					if not variant_pos == None:
						#variant 비율 계산을 위한 작업
						cluster_var = 0 #cluster 그룹에 속하는 특정위치의 variant 개수
						cnt_var = int(pos_var[pos_var['pos']==variant_pos]['cnt']) # 입력위치에 variant를 가진 sample의 수
							
						# 입력위치에  variant가 있는 sample을 하나씩 볼거에요
						for i, var in df_var[df_var['pos']==variant_pos].iterrows() :
							# cluster 된 sample일 때
							if boths_01[int(var['index'])] == 1 :						
								cluster_var = cluster_var+1
								var_index.append(int(var['index']))
								variant_list.append(var['MSID'])

						var_ratio = cluster_var/cnt_var
					
						ratio_list.append(var_ratio)
						ratio_list_t.append("{0:.2f}".format(var_ratio))
						print(ci,tau,var_ratio,score)
					
					else:
						ratio_list.append(len(boths_t)/len(coverage[0]))
						ratio_list_t.append("{0:.2f}".format(len(boths_t)/len(coverage[0])))
						print(ci,tau,len(boths_t)/len(coverage[0]),score)
					
					# silhouette score가 highest score보다 높다면 결과 갱신
					if score > highest_score :
						print("Update highest score")
						drops = drops_t[:]
						rises = rises_t[:]
						boths = boths_t[:]

						highest_ci = ci
						highest_tau = tau
						highest_tnum = tn
						highest_ratio = len(boths)/len(coverage[0]) if variant_pos == None else var_ratio
			
		
			# 입력된 tau 값이 없을 때	
			else:
				# Tau 계산
				for tnum in range(0, 41, 1) :
					tau = 0
				
					# drop sample 중 평균이 가장 높은 값부터 Tau를 설정.
					if tnum > 0 and len(drop_means_t) > 1 :
						# 오른쪽 값이 제일 큰 sample 선택 후 삭제
						hdmt = drop_means_t[0]
						del drop_means_t[0]
						
						# 해당 오른쪽 값을 tau로 지정후 drops에서 삭제
						if hdmt[0] in drops_t :
							drops_t.remove(hdmt[0])
							tau = hdmt[2]

					# tau의 값이 제한에 걸리지 않으면 pass, 제한이 걸려 있으면 해당 부분 생략
					if not limit_tau == None and (tnum == 0 or tau > limit_tau) :
						continue
					if not limit_tau_low == None and (tnum == 0 or tau < limit_tau_low) :
						continue

	
					# 임시 sample로 boths 계산
					boths_t = []
					boths_01 = np.zeros(len(coverage[0]))	
					# drops_t, rises_t를 비교하여 일치하면 boths_t로 분류
					for dn in drops_t :
						if dn in rises_t :
							boths_t.append(dn)
							boths_01[dn] = 1

					if len(boths_t) < 1 :
						continue
	

					# silhouette score를 계산할 부분 추출
					df_silhouette = df_drop.iloc[(e1_l-dintv):(e1_l+1)] 
					df_silhouette = df_silhouette.append(df_drop.iloc[e1_r:(e1_r+dintv+1)]) 
					df_silhouette = df_silhouette.append(df_drop.iloc[(e2_l-dintv):(e2_l+1)])
					df_silhouette = df_silhouette.append(df_drop.iloc[e2_r:(e2_r+dintv+1)])
					df_silhouette = df_silhouette.T.values
					
					#각 구간의 silhoutte score 계산/ silhouette_score(sample array,label)
					score = silhouette_score(df_silhouette, boths_01)

					# 그래프를 그리기 위해 해당 결과 저장
					ci_list.append(ci)
					tau_list.append(tau)
					score_list.append(score)

					#cluster 그룹에 속하는 variant의 비율을 알고 싶을 때
					if not variant_pos == None:
						#variant 비율 계산을 위한 작업
						cluster_var = 0 #cluster 그룹에 속하는 특정위치의 variant 개수
						cnt_var = int(pos_var[pos_var['pos']==variant_pos]['cnt']) # 입력위치에 variant를 가진 sample의 수
				       
						 # 입력위치에  variant가 있는 sample을 하나씩 볼거에요
						for i, var in df_var[df_var['pos']==variant_pos].iterrows() :
							# cluster 된 sample일 때
							if boths_01[int(var['index'])] == 1 :
								cluster_var = cluster_var+1
								var_index.append(int(var['index']))
								variant_list.append(var['MSID'])

						var_ratio = cluster_var/cnt_var
	
						ratio_list.append(var_ratio)
						ratio_list_t.append("{0:.2f}".format(var_ratio))
	
						print(ci,tau,var_ratio,score)

					else:
						ratio_list.append(len(boths_t)/len(coverage[0]))
						ratio_list_t.append("{0:.2f}".format(len(boths_t)/len(coverage[0])))

						print(ci,tau,len(boths_t)/len(coverage[0]), score)


					# silhouette score가 highest score보다 높다면 결과 갱신
					if score > highest_score :
						print("Update higest_score")
						drops = drops_t[:]
						rises = rises_t[:]
						boths = boths_t[:]

						highest_score = score
						highest_ci = ci
						highest_tau = tau
						highest_tnum = tnum
						highest_ratio = len(boths)/len(coverage[0]) if variant_pos == None else var_ratio


		# for 문 끝
		print('-------------------------------------------------')

		# 최종 결과 출력
		print(highest_ci, highest_tnum, highest_score, highest_tau, highest_ratio)

		print('drops : '+str(len(drops)))
		print('rises : '+str(len(rises)))
		print('boths : '+str(len(boths)))


		# boths 01 계산
		boths_01 = np.zeros(len(coverage[0]))
		for dn in boths :
			boths_01[dn] = 1
	

		# Make ci/tau/score scatter plot

		print('silhouette plot drawing...')

		# plot size & font
		plt.figure(figsize=(score_plot_width, score_plot_height))
		score_fontsize = np.sqrt(score_plot_width*score_plot_height)
		score_fontsize_big = np.sqrt(score_plot_width*score_plot_height)*2

		# scatter size
		size_cmap = plt.cm.viridis
		size_list = [pow((s+1), 4)*100 for s in score_list]


		# make scatter plot
		plt.scatter(ci_list, tau_list, s=size_list, c=score_list, cmap=size_cmap, alpha=0.5)
		# 우측에 color bar 생성 후 폰트 크기 조정
		cbar = plt.colorbar()
		cbar.ax.tick_params(labelsize=score_fontsize)

		# 선택된scatter plot의 dot에 강조표시함
		plt.scatter([highest_ci], [highest_tau], s=[pow((highest_score+1), 4)*100], c=['None'], edgecolor='red', linewidth='1.5')


		# make texts
		for tn, text in enumerate(ratio_list):
			plt.text(ci_list[tn], tau_list[tn], ratio_list_t[tn], color='r',
				fontsize=score_fontsize/3, ha='center', va='center')


		# scatter x & y label
		plt.xlabel('CI', fontsize=score_fontsize_big)
		plt.ylabel('Tau', fontsize=score_fontsize_big)
		
		plt.xticks(range(0,int(max(ci_list))+1), fontsize=score_fontsize)
		plt.yticks(fontsize=score_fontsize)
		# initialize legends
		legends = []
		scores = []


		# make legends
		for lnum in range(5) :
			# size와 score list 불러오기
			size_sort = sorted(size_list[:])
			score_sort = [round(pow(s/100, 1/4)-1, 2) for s in size_sort]
	
			# legends를 5분위로 나누어 저장
			lnum_size = size_sort[int((len(size_sort)-1)/4*lnum)]
			lnum_score = score_sort[int((len(score_sort)-1)/4*lnum)]
			scores.append(str(lnum_score))
			legends.append(plt.scatter([], [], s=lnum_size, c=size_cmap((1/4)*lnum), alpha=0.5))
		
		# legends 출력
		plt.legend(tuple(legends), tuple(scores), loc='upper right', title='Scores', fontsize=score_fontsize)
		
		# adjust space
		plt.subplots_adjust(wspace=0, hspace=0)


		# save ci/tau/ratio scatter plot
		plt.savefig(output_prefix+'_silhouette_plot.pdf')
		plt.close()
		print(output_prefix+'_silhouette_plot.pdf saved!')

		print()

		
		#Make ratio(proportion)/score scatter plot
	
		print('ratio/score scatter plot drawing...')

		plt.scatter(ratio_list, score_list, s=9)
		plt.xlabel('ratio')
		plt.ylabel('score')
		plt.xticks(np.linspace(0.0, 1.0, num=10))
		plt.savefig(output_prefix+'_ratio_score_plot.pdf')
		plt.close()
		print(output_prefix+'_ratio_score_plot.pdf saved!')
		
		print()


	# nmf 계산
	elif filt_mode == 'nmf':

		# nmf는 exon 사이 간격으로 최적화 시킨다.
		highest_dintv = 0
		highest_score = -1

		for dintv in range(5, 31, 1) :
			boths_t = []
			boths_01_t = []

			# silhouette score 그릴 부분 추출
			df_silhouette = df_drop.iloc[(e1_l-dintv):(e1_l+1)]
			df_silhouette = df_silhouette.append(df_drop.iloc[e1_r:(e1_r+dintv+1)])
			df_silhouette = df_silhouette.append(df_drop.iloc[(e2_l-dintv):(e2_l+1)])
			df_silhouette = df_silhouette.append(df_drop.iloc[e2_r:(e2_r+dintv+1)])
			df_silhouette = df_silhouette.T.values

			# NMF 실행
			model = NMF(n_components=2, init='random', random_state=0)
			W = model.fit_transform(df_silhouette.T)
			H = model.components_

			# W(n*2) 행렬에서 1 column이 2column보다 작으면 clustering.
			for wi, wl in enumerate(W):
				if wl[0] < wl [1] :
					boths_t.append(wi)

			# boths 01 생성
			boths_01_t = np.zeros(len(coverage[0]))
			for dn in boths_t :
				boths_01_t[dn] = 1
				
			# silhoutte score
			score = silhouette_score(df_silhouette, boths_01_t)

			print(dintv, " : ", score)

			# score가 highest score보다 크면 갱신
			if score > highest_score :
				highest_score = score
				highest_dintv = dintv
				boths = boths_t[:]
				boths_01 = boths_01_t[:]

		if not variant_pos == None:
			# 입력위치에  variant가 있는 sample을 하나씩 볼거에요
			for i, var in df_var[df_var['pos'] == variant_pos].iterrows():
				# cluster 된 sample일 때
				if boths_01[int(var['index'])] == 1:
					var_index.append(int(var['index']))
					variant_list.append(var['MSID'])

		print()
		print(highest_dintv, highest_score)
		print(len(boths)/len(boths_01))


	# splice site로 실행
	elif filt_mode == 'splice_site':
		df_var_p = df_var[:]

		m_pos_list = []
		m_mar_list = []
		
		# splice site중 effect에 Splice Site가 포함되어 있으면
		# 해당하는 bam을 추가
		for vp_i, vp in df_var_p.iterrows() :
			if 'SPLICE_SITE' in str(vp['effect']) :
				boths.append(int(vp['index']))
		
		# boths 생성
		boths = sorted(list(set(boths)))
		boths_01 = np.zeros(len(coverage[0]))
		for dn in boths :
			boths_01[dn] = 1

		if not variant_pos == None:
			# 입력위치에  variant가 있는 sample을 하나씩 볼거에요
			for i, var in df_var[df_var['pos'] == variant_pos].iterrows():
				# cluster 된 sample일 때
				if boths_01[int(var['index'])] == 1:
					var_index.append(int(var['index']))
					variant_list.append(var['MSID'])

		print(len(boths),len(boths_01))


	# scatter plot을 그려서 clustering
	elif filt_mode == 'scatter':
		scatter_tau =[]
		scatter_ci = []
		v_tau = []
		v_ci = []

		boths_01 = np.zeros(len(coverage[0]))
		space = input_space
		
		# 첫번째 선택 exon에 대해 ci(gap), tau(cr_mean) 계산
		for cn, c in enumerate(df_drop.columns):
			cluster1 = df_drop.iloc[(e1_l-space):(e1_l+1)][c]
			cluster2 = df_drop.iloc[e1_r:(e1_r+space+1)][c]
			
			cl_mean = cluster1.mean()
			cr_mean = cluster2.mean()
			
			gap = cl_mean-cr_mean

			scatter_tau.append(cr_mean)
			scatter_ci.append(gap)
		
			if (gap >= x_div) & (cr_mean <= y_div):
				boths.append(cn)
				boths_01[cn] = 1

		if not variant_pos == None:
			# 입력위치에  variant가 있는 sample을 하나씩 볼거에요
			for i, var in df_var[df_var['pos'] == variant_pos].iterrows():
				v_index = int(var['index'])
				# cluster 된 sample일 때
				if boths_01[v_index] == 1:
					var_index.append(v_index)
					variant_list.append(var['MSID'])
					v_tau.append(scatter_tau[v_index])
					v_ci.append(scatter_ci[v_index])					

		
		df_silhouette = df_drop.iloc[(e1_l-dintv):(e1_l+1)]
		df_silhouette = df_silhouette.append(df_drop.iloc[e1_r:(e1_r+dintv+1)])
		df_silhouette = df_silhouette.append(df_drop.iloc[(e2_l-dintv):(e2_l+1)])
		df_silhouette = df_silhouette.append(df_drop.iloc[e2_r:(e2_r+dintv+1)])
		df_silhouette = df_silhouette.T.values

		score = silhouette_score(df_silhouette, boths_01)
	
		print("silhouette score :", score)
		print("number of boths :", len(boths), "ratio of boths :", len(boths)/len(coverage[0]))


		# ci/tau scatter plot을 그려보자
		print('scatter plot drawing...')

		# plot size
		plt.figure(figsize=(score_plot_width, score_plot_height))

		# scatter size
		size_cmap = plt.cm.viridis

		# make scatter plot
		plt.scatter(scatter_ci, scatter_tau, s=marker_size*10 , c='#7ac5cd', edgecolor='green', linewidth=1.5, cmap=size_cmap, alpha=0.5)

		# variant sample의 plot을 강조표시함
		plt.scatter(v_ci, v_tau, s=marker_size*10, c='#e52a59')

		# scatter x & y label
		plt.xlabel('CI', fontsize=font_size*2.5)
		plt.ylabel('Tau', fontsize=font_size*2.5)

		plt.xticks(np.arange(int(min(scatter_ci)) - 1, int(max(scatter_ci)) + 1, 1), fontsize=font_size*1.5)
		plt.yticks(np.arange(int(min(scatter_tau)) - 1, int(max(scatter_tau)) + 1, 1), fontsize=font_size*1.5)

		# plot 저장
		plt.savefig(output_prefix + '_scatter_plot.pdf', bbox_inches='tight', pad_inches=3)
		plt.close()
		print(output_prefix + '_sctter_plot.pdf saved!')
		print()


	# filter mode를 잘못 입력했을 때
	else :
		print('Wrong Filtering Mode : '+filt_mode)
		exit()


	# variant_pos 입력이 있을 때 해당 variant에 대한 정보 출력
	if not variant_pos == None:
		v = df_var[df_var['pos']==variant_pos]['MSID'].values.tolist()

		print("\n*List of samples with variant on "+str(variant_pos)+' :', v)
		print("*Variant list in red cluster : ", variant_list,'\n')


	# 각 sample의 info(MSID, cluster group,exon2 mutation,exon2 splice)를 저장할 list 생성
	sample_info = [[] for i in range(len(bam_names))]

	for i, c in enumerate(boths_01):
		if c == 1:
			sample_info[i].extend([str(bam_names[i]), 'Red'])
		else:
			sample_info[i].extend([str(bam_names[i]), 'Green'])

	# 각 sample이 exon2에서 mutation, splice_site을 갖는지 조사
	for i, sample in enumerate(sample_info):
		mutant = 0
		splice = 0

		for n, var in df_var[df_var['MSID']==bam_names[i]].iterrows():
			# exon 2 앞뒤 100bp 범위로 mutation 조사 
			if (int(Exon_s_list[1])-100 <= int(var['pos']) and int(var['pos']) <= int(Exon_e_list[1])+100):
				mutant = mutant+1
				# splice site variant 조사
				if 'splice' in str(var['effect']).lower():
					splice = splice+1

		sample.extend([mutant,splice])
	
	df_sample_info= pd.DataFrame(sample_info, columns=['MSID','Cluster group','mutation on exon2','splice on exon2'])
	df_sample_info.to_excel(output_prefix+'_sample_info.xlsx', index=False)
	print(output_prefix+"_sample_info.xlsx saved!")
	print()
	

	
	
## Draw Plots of coverage and refseq ##

# 최종 출력 refseq 갯수
nl = len(name_nm)

# 전체적인 font size 설정
matplotlib.rcParams.update({'font.size': font_size})

# max_whole(coverage plot의 y축)
c_max=[]

if max_whole == 0:
	for i, s in enumerate(draw_range):
		e = draw_range_e[i]
		df2_mean_t = df2_mean[s-draw_start:e-draw_start].values.tolist()
		c_max.append(max(map(max,df2_mean_t)))
	
	max_whole = max(c_max)



# combine_slices 옵션에서 전체 plot 그리기
if combine :
	# figure 설정
	# 여러 칸으로 나눈다.
	fig = plt.figure(figsize=(50, 20))
	gs = gridspec.GridSpec(nrows=2, ncols=len(draw_range), height_ratios=[4+(nl-1)*0.2, nl])
	gs.update(wspace=0, hspace=0.1)
	fig.suptitle(title)

	for n, st in enumerate(draw_range) :
		sys.stderr.write('\r'+'exon '+str(n+1)+' drawing...')
		sys.stderr.write("\033[K")

		stop_n = draw_range_e[n]
	
		# x축의 값을 정해진 시작점과 끝점으로 한다.
		xticks = np.arange(st, stop_n)

		ax_main = plt.subplot(gs[n])
		ax_main.yaxis.set_tick_params(labelsize=font_size*2)

		if not n == 0 :
			ax_main.get_yaxis().set_ticks([])
			ax_main.spines['left'].set_visible(False)
		if not n == len(draw_range)-1 :
			ax_main.spines['right'].set_color('#CDCDCD')

		
		# set average window
		df2_mean_p = df2_mean[st-draw_start:stop_n-draw_start]
		df2_T = df2_mean_p.T
		data_lines = df2_T.values





		## Line plot ##
		# clustering 안했을 시 전체 그래프를 그린다.
		if not view_mode :	
			if min_max :
				ax_main.fill_between(xticks, df2_mean_p.max(axis=1), df2_mean_p.min(axis=1), facecolor='g', alpha=0.5)
				ax_main.plot(df2_mean_p.max(axis=1), color='g', linewidth=3.0)
				ax_main.plot(df2_mean_p.min(axis=1), color='g', linewidth=3.0)
			else :
				ax_main.plot(df2_mean_p, color='g', alpha=0.5)
			if draw_average_line :
				ax_main.plot(df2_mean_p.mean(axis=1), color='yellow',linewidth=3.0)
		
			# x,y 축 설정
			ax_main.set_ylim([0,max_whole])
			ax_main.set_xlim([st,stop_n])

		else :
			# boths가 아닌 부분 draw
			df_nb = df2_mean_p.iloc[:, [b for b in range(len(coverage[0])) if not b in boths]]

			if min_max :
				ax_main.fill_between(xticks, df_nb.max(axis=1), df_nb.min(axis=1), facecolor='g', alpha=0.5)
				ax_main.plot(df_nb.max(axis=1), color='g', linewidth=3.0)
				ax_main.plot(df_nb.min(axis=1), color='g', linewidth=3.0)
			else :
				ax_main.plot(df_nb, color='g', alpha=0.5)	
			if draw_average_line:
				ax_main.plot(df_nb.mean(axis=1),color='#006600',linewidth=5.0,linestyle='--')

			# y축과 x축 표시할 부분 지정
			ax_main.set_ylim([0, max_whole])
			ax_main.set_xlim([st, stop_n])

			real_bam_list_s = [r[:r.find('.')] for r in real_bam_list]
						
			if len(boths) > 0:
			# boths인 부분 draw
				df_b = df2_mean_p.iloc[:, boths]
			
				if min_max :
					ax_main.fill_between(xticks, df_b.max(axis=1), df_b.min(axis=1), facecolor='red', alpha=0.5)
					ax_main.plot(df_b.max(axis=1), color='red', linewidth=3.0)
					ax_main.plot(df_b.min(axis=1), color='red', linewidth=3.0)
				else :
					ax_main.plot(df_b, color='r', alpha=0.5)
				if draw_average_line:
					ax_main.plot(df_b.mean(axis=1), color='#cc0033', linewidth=5.0,linestyle='--')
		




		## special Plot print ##
		
#		ax_main.plot(df_special.iloc[st-draw_start:stop_n-draw_start], color='orange', linewidth = 3.0)
		




		## Draw Genetic Variants Data ##

		df_var_p = df_var[df_var['pos'] >= st]
		df_var_p = df_var_p[df_var_p['pos'] < stop_n]

		c_m = '#A0ffA0'
		marker_colors = ['red', 'g']
		m_pos_list = []
		m_col_list = []
		
		# genetic variants 그리고 각 데이터 위치와 색상 저장
		for vp_i, vp in df_var_p.iterrows() :

			# 마커 모양 지정 
			var_mark = var_markers[effect_list.index(vp['effect'])%len(var_markers)]
			# 마커 색상 지정
			var_col = ''
			if boths_01[int(vp['index'])] == 1 :
				var_col = 'red'
			else :
				var_col = 'g'
			# 마커 위치와 색상 저장
			m_pos_list.append(vp['pos'])
			m_col_list.append(var_col)
			# 마커 그리기
			ax_main.plot(vp['pos'], coverage[int(vp['pos']-start)][int(vp['index'])], marker=var_mark,
						markeredgecolor='black', markeredgewidth=2, color=var_col, markersize=marker_size, alpha=0.5)

		# 마커 위치 중복 제거 후 색상 비율을 결정하는 배열 생성
		marker_positions = list(set(m_pos_list))
		marker_ratio = [np.zeros(len(marker_colors)) for i in marker_positions]

		# 비율 표시 위치 지정
		marker_y = ax_main.get_ylim()[1]/100	

		# 각 마커의 정보를 marker_ratio배열에 저장
		for m_i, m_pos in enumerate(m_pos_list) :		
			# 해당 variant의 색상 불러오기
			m_col = m_col_list[m_i]

			# 위치와 색상 지정
			if m_pos in marker_positions :
				index_p = marker_positions.index(m_pos)

			if m_col in marker_colors :
				index_c = marker_colors.index(m_col)
			
			marker_ratio[index_p][index_c] += 1

		# 배열 정보에 따라 비율 표시
		for r_i, ratio_list in enumerate(marker_ratio) :
			r1 = ratio_list[0] / sum(ratio_list)

			# non-boths 마커 모양 설정
			x = [0] + np.cos(np.linspace(0, 2 * np.pi * r1, 10)).tolist()
			y = [0] + np.sin(np.linspace(0, 2 * np.pi * r1, 10)).tolist()
			xy1 = np.column_stack([x, y])

			# boths 마커 모양 설정
			x = [0] + np.cos(np.linspace(2 * np.pi * r1, 2 * np.pi, 10)).tolist()
			y = [0] + np.sin(np.linspace(2 * np.pi * r1, 2 * np.pi, 10)).tolist()
			xy2 = np.column_stack([x, y])

			# 마커 표시
			if not ratio_list[0] == 0 : 
				ax_main.plot(marker_positions[r_i], marker_y, marker=xy1, markersize=marker_size, color='#ffA0A0')
			if not ratio_list[1] == 0 : 
				ax_main.plot(marker_positions[r_i], marker_y, marker=xy2, markersize=marker_size, color=c_m)

			# 마커 테두리 표시
			ax_main.plot(marker_positions[r_i], marker_y, color='None', marker='o', markersize=marker_size, markeredgecolor='#888888', markeredgewidth=2)


		xt_step = (stop_n-st)/4
		plt.xticks(np.arange(st+xt_step, stop_n, step=xt_step)) # x축 값 표시.

		# bp 위치의 값이 너무 크기 때문에 e를 포함한 값으로 표시된다.
		# 따라서 아래와 같은 과정을 거쳐 자연수가 나오도록 한다.
		xx, locs = plt.xticks()
		ll = ['%d' % a for a in xx]
		plt.xticks(xx, ll)




		
		## 각 구간의 Clustering evaluation plot ##
		if view_mode :
			clustered_list = boths_01

			# silhouette를 각 sample에 대해 계산
			silhouette_avg = silhouette_score(data_lines, clustered_list)
			sample_silhouette_values = silhouette_samples(data_lines, clustered_list)

			# inset 그래프 설정
			ax_inset = inset_axes(ax_main, width="30%", height="50%", loc=2)
			ax_inset.set_ylim([0, len(data_lines) + (2 + 1) * 10])
			
			# 각 분류끼리의 간격 설정
			y_lower = 10

			# clustering evaluation drawing
			for i in range(2):
				# i번째 분류에 대한 silhouette 계산 및 정렬
				ith_cluster_silhouette_values = sample_silhouette_values[clustered_list == i]
				ith_cluster_silhouette_values.sort()

				# sample 갯수 계산 뒤 그릴 범위지정
				size_cluster_i = ith_cluster_silhouette_values.shape[0]
				y_upper = y_lower + size_cluster_i

				# i가 1이면 boths이다.
				color = 'red' if i == 1 else 'g'
				# 해당 분류의 그래프를 채운다.
				ax_inset.fill_betweenx(np.arange(y_lower, y_upper),
							0, ith_cluster_silhouette_values,
							facecolor=color, edgecolor=color, alpha=0.5)

				# 해당 분류에 텍스트를 단다.
				ax_inset.text(-0.05, y_lower + 0.5 * size_cluster_i, 'clustered' if i == 1 else 'non-clustered')

				y_lower = y_upper + 10

			# silhouette 값의 평균을 표시하고 0 값을 표시한다.
			ax_inset.axvline(x=silhouette_avg, color="red", linestyle="--")
			ax_inset.axvline(x=0, color="black", linestyle="-")

			ax_inset.set_yticks([])  # Clear the yaxis labels / ticks





		## Draw Refseq plot ##

		# refseq를 표시할 범위를 정한다.
		byts = range(2-nl, 2)

		# refseq가 표시될 이름을 정한다.
		ytlbs = [name_n+"\n"+name2_nm[0]]
		ax_bottom = plt.subplot(gs[n+len(draw_range)], yticks=byts, xticklabels=[], yticklabels=list(reversed(ytlbs)))
		ax_bottom.yaxis.set_tick_params(labelsize=font_size*1.5)
		ax_bottom.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False)
		ax_bottom.set_ylim(2-nl-0.4, 1+0.4)

		if not n == 0 :
			ax_bottom.get_yaxis().set_ticks([])
			ax_bottom.spines['left'].set_visible(False)
		if not n == len(draw_range)-1 :
			ax_bottom.spines['right'].set_color('#CDCDCD')
		ax_bottom.set_xlim([st, stop_n])
		
		# refseq의 범위를 검은 선으로 그린다.
		bxts = np.arange(draw_start, draw_stop)
		blns = np.full(draw_stop-draw_start,  1, dtype=int)
		ax_bottom.plot(bxts, blns, 'black')

		# tx, cds의 시작 부분을 얇은 네모로 그린다.
		cds_s = cds_s_nm[0]
		if (st <= start and cds_s<=stop_n):
			rect = patches.Rectangle((start, 0.9),cds_s-start,0.2,edgecolor='none', facecolor='black')		
			ax_bottom.add_patch(rect)
		
		# tx, cds의 끝 부분을 얇은 네모로 그린다.
		cds_e = cds_e_nm[0]
		if (st < cds_e and stop <= stop_n):
			rect = patches.Rectangle((cds_e,0.9), stop-cds_e, 0.2, edgecolor = 'none', facecolor='black')
			ax_bottom.add_patch(rect)
	

		# 엑손을 그린다.	
		for k, ess in enumerate(ess_nm) :
			if (ess > stop_n or ees_nm[k] < st) :
				continue
			rect = patches.Rectangle((ess, 0.8),ees_nm[k]-ess,0.4,edgecolor='none',facecolor='black')
			ax_bottom.add_patch(rect)
			leftt = ess if ess > st else st
			rightt = ees_nm[k] if ees_nm[k] < stop_n else stop_n
			ax_bottom.text((leftt+rightt)/2, 1, str(exons[k]), horizontalalignment='center', verticalalignment='center', color='white', fontsize=font_size*3)





	## Genetic variants legend 표시##
	var_legends = []
	var_names = []

	# 각 variant의 모양과 이름 수집 후 legend에 저장
	for ef_i, ef in enumerate(effect_list):
		var_mark = var_markers[ef_i%len(var_markers)]
		var_col = var_colors[ef_i%len(var_colors)]
		var_names.append(ef)
		var_legends.append(ax_main.scatter([], [], marker=var_mark, color=var_col, edgecolor='black', s=marker_size**2, linewidth='2', alpha=0.5))

	# 저장 된 것을 바탕으로 legend 표시
	ax_main.legend(tuple(var_legends), tuple(var_names), loc='upper right', title='Genetic Variants')

	print()

	# lineplot 간 여백 설정
	plt.subplots_adjust(top = 0.95, bottom = 0.05, right = 0.97, left = 0.03, wspace=0, hspace=0)





	## pdf로 저장##
	plt.savefig(output_prefix+'.pdf',bbox_inches='tight',pad_inches=3)
	plt.close(fig)

	print(output_prefix+'.pdf saved!')
	




## combine_slices 옵션없이 plot 그리기 ##
else :  

	## Draw range 범위만큼 slice 출력##
	for n, st in enumerate(draw_range) :

		print('\n'+output_prefix+'_'+str(n+1)+' saving...')

		# 그래프가 끝나는 지점을 정한다.
		# refseq가 끝나는 지점보다 크면 stop+1을 끝점으로 한다.
		stop_n = 0

		if flag :
			stop_n = draw_stop+1 if draw_range_e[n] >= draw_stop+1 else draw_range_e[n]
		else :
			stop_n = draw_stop+1 if st+draw_span >= draw_stop+1 else st+draw_span

		# x축의 값을 정해진 시작점과 끝점으로 한다.
		xticks = np.arange(st, stop_n)

		# figure 초기화(크기 조정, 간격 조정)
		fig2 = plt.figure(figsize=(30, 2+12+1*nl))
		gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[4+(nl-1)*0.2, nl])
		gs.update(wspace=0, hspace=2)

		ax_main = plt.subplot(gs[0])	

		# set average window

		df2_mean_p = df2_mean[st-draw_start:stop_n-draw_start]
		df2_T = df2_mean_p.T
		data_lines = df2_T.values





		## Line plot ##
		# clustering 안했을 시 전체 그래프를 그린다.
		if not view_mode :
			if min_max :
				ax_main.fill_between(xticks, df2_mean_p.max(axis=1), df2_mean_p.min(axis=1), facecolor='g', alpha=0.5)
				ax_main.plot(df2_mean_p.max(axis=1), color='g', linewidth=3.0)
				ax_main.plot(df2_mean_p.min(axis=1), color='g', linewidth=3.0)
			else :
				ax_main.plot(df2_mean_p, color='g', alpha=0.5)
			if draw_average_line :
				ax_main.plot(cv_mean[st-draw_start:stop_n-draw_start],color='yellow',linewidth=3.0)
						
			ax_main.set_xlim([st,stop_n])			

		else :
			# both가 아닌부분 draw
			df_nb = df2_mean_p.iloc[:, [b for b in range(len(coverage[0])) if not b in boths]]

			if min_max : 
				ax_main.fill_between(xticks, df_nb.max(axis=1), df_nb.min(axis=1), facecolor='g', alpha=0.5)
				ax_main.plot(df_nb.max(axis=1), color='g', linewidth=3.0)
				ax_main.plot(df_nb.min(axis=1), color='g', linewidth=3.0)
			else :
				ax_main.plot(df_nb, color='g', alpha=0.5)	
			if draw_average_line :
				ax_main.plot(df_nb.green(axis=1),color='#006600',linewidth=3.0)

			# y축과 x축 표시할 부분 지정
			ax_main.set_xlim([st, stop_n])

			if len(boths) > 0:
			# boths인 부분 draw
				df_b = df2_mean_p.iloc[:, boths]

				if min_max :
					ax_main.fill_between(xticks, df_b.max(axis=1), df_b.min(axis=1), facecolor='red', alpha=0.5)
					ax_main.plot(df_b.max(axis=1), color='red', linewidth=3.0)
					ax_main.plot(df_b.min(axis=1), color='red', linewidth=3.0)
				else :
					ax_main.plot(df_b, color='r', alpha=0.5)
				if draw_average_line:
					ax_main.plot(df_b.mean(axis=1), color='#cc0033', linewidth=3.0)			





		## Special Plot print ##

#		ax_main.plot(df_special.iloc[st-draw_start:stop_n-draw_start], color='cyan', linewidth = 3.0)


		


		## Draw Genetic Variants Data ##
		df_var_p = df_var[df_var['pos'] >= st]
		df_var_p = df_var_p[df_var_p['pos'] < stop_n]

		c_m = '#A0ffA0'
		marker_colors = ['red', 'g']
		m_pos_list = []
		m_col_list = []
		
		# genetic variants 그리고 각 데이터 위치와 색상 저장
		for vp_i, vp in df_var_p.iterrows() :
			# 마커 모양 지정 
			var_mark = var_markers[effect_list.index(vp['effect'])%len(var_markers)]

			# 마커 색상 지정
			var_col = ''
			if boths_01[int(vp['index'])] == 1 :
				var_col = 'red'
			else :
				var_col = 'g'
		
			# 마커 위치와 색상 저장
			m_pos_list.append(vp['pos'])
			m_col_list.append(var_col)
			# 마커 그리기
			ax_main.plot(vp['pos'], coverage[int(vp['pos']-draw_start)][int(vp['index'])], marker=var_mark,
						markeredgecolor='black', markeredgewidth=2, color=var_col, markersize=marker_size, alpha=0.5)

		# 마커 위치 중복 제거 후 색상 비율을 결정하는 배열 생성
		marker_positions = list(set(m_pos_list))
		marker_ratio = [np.zeros(len(marker_colors)) for i in marker_positions]

		# 비율 표시 위치 지정
		marker_y = ax_main.get_ylim()[1]/100	

		# 각 마커의 정보를 marker_ratio배열에 저장
		for m_i, m_pos in enumerate(m_pos_list) :
			# 해당 variant의 색상 불러오기
			m_col = m_col_list[m_i]

			# 위치와 색상 지정
			if m_pos in marker_positions :
				index_p = marker_positions.index(m_pos)

			if m_col in marker_colors :
				index_c = marker_colors.index(m_col)
			
			marker_ratio[index_p][index_c] += 1

		# 배열 정보에 따라 비율 표시
		for r_i, ratio_list in enumerate(marker_ratio) :
			r1 = ratio_list[0] / sum(ratio_list)

			# non-boths 마커 모양 설정
			x = [0] + np.cos(np.linspace(0, 2 * np.pi * r1, 10)).tolist()
			y = [0] + np.sin(np.linspace(0, 2 * np.pi * r1, 10)).tolist()
			xy1 = np.column_stack([x, y])

			# boths 마커 모양 설정
			x = [0] + np.cos(np.linspace(2 * np.pi * r1, 2 * np.pi, 10)).tolist()
			y = [0] + np.sin(np.linspace(2 * np.pi * r1, 2 * np.pi, 10)).tolist()
			xy2 = np.column_stack([x, y])

			# 마커 표시
			if not ratio_list[0] == 0 : 
				ax_main.plot(marker_positions[r_i], marker_y, marker=xy1, markersize=marker_size, color='#ffA0A0')
			if not ratio_list[1] == 0 : 
				ax_main.plot(marker_positions[r_i], marker_y, marker=xy2, markersize=marker_size, color=c_m)

			# 마커 테두리 표시
			ax_main.plot(marker_positions[r_i], marker_y, color='None', marker='o', markersize=marker_size, markeredgecolor='#888888', markeredgewidth=2)
	

		plt.xticks(np.arange(st, stop_n+1, step=(stop_n-st)/10)) # x축 값 표시.

		# bp 위치의 값이 너무 크기 때문에 e를 포함한 값으로 표시된다.
		# 따라서 아래와 같은 과정을 거쳐 자연수가 나오도록 한다.
		xx, locs = plt.xticks()
		ll = ['%d' % a for a in xx]
		plt.xticks(xx, ll)





		## 각 구간의 Clustering evaluation plot ##
		if view_mode :
			clustered_list = boths_01

			# silhouette를 각 sample에 대해 계산
			silhouette_avg = silhouette_score(data_lines, clustered_list)
			sample_silhouette_values = silhouette_samples(data_lines, clustered_list)

			# inset 그래프 설정
			ax_inset = inset_axes(ax_main, width="30%", height="50%", loc=2)
			ax_inset.set_ylim([0, len(data_lines) + (2 + 1) * 10])

			# 각 분류끼리의 간격 설정
			y_lower = 10

			# clustering evaluation drawing
			for i in range(2):
				# i번째 분류에 대한 silhouette 계산 및 정렬
				ith_cluster_silhouette_values = sample_silhouette_values[clustered_list == i]
				ith_cluster_silhouette_values.sort()

				# sample 갯수 계산 뒤 그릴 범위지정
				size_cluster_i = ith_cluster_silhouette_values.shape[0]
				y_upper = y_lower + size_cluster_i

				# i가 1이면 boths이다.
				color = 'red' if i == 1 else 'g'
				# 해당 분류의 그래프를 채운다.
				ax_inset.fill_betweenx(np.arange(y_lower, y_upper),
					0, ith_cluster_silhouette_values,
					facecolor=color, edgecolor=color, alpha=0.5)

				# 해당 분류에 텍스트를 단다.
				ax_inset.text(-0.05, y_lower + 0.5 * size_cluster_i, 'clustered' if i == 1 else 'non-clustered')

				y_lower = y_upper + 10

			# silhouette 값의 평균을 표시하고 0 값을 표시한다.
			ax_inset.axvline(x=silhouette_avg, color="red", linestyle="--")
			ax_inset.axvline(x=0, color="black", linestyle="-")

			ax_inset.set_yticks([])  # Clear the yaxis labels / ticks



	

		## Draw Refseqs ##
		# refseq를 표시할 범위를 정한다.
		byts = range(2-nl, 2)

		# refseq가 표시될 이름을 정한다.
		ytlbs = [name_n+"\n"+name2_nm[0]]
		ax_bottom = plt.subplot(gs[1], yticks=byts, xticklabels=[], yticklabels=list(reversed(ytlbs)))
		ax_bottom.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False)

		# refseq의 범위를 검은 선으로 그린다.
		a_s = draw_start if draw_start > st else st
		a_e = draw_stop if draw_stop < stop_n else stop_n

		if a_e-a_s < 0:
			continue

		bxts = np.arange(a_s,a_e)
		blns = np.full(a_e-a_s,1,dtype=int)
		ax_bottom.plot(bxts,blns,'black')

		# tx, cds의 시작 부분을 얇은 네모로 그린다.
		cds_s = cds_s_nm[0]
		if (st<=start and cds_s <= stop_n):
			rect = patches.Rectangle((start,0.9),cds_s-start,0.2,edgecolor='none',facecolor='black')
			ax_bottom.add_patch(rect)

		# tx, cds의 끝 부분을 얇은 네모로 그린다.
		cds_e = cds_e_nm[0]
		if (st<=cds_e and stop <= stop_n) :
			rect = patches.Rectangle((cds_e, 0.9), stop-cds_e, 0.2, edgecolor='none', facecolor='black')
			ax_bottom.add_patch(rect)	

		# 방향을 그린다
		if(start > stop_n or stop < st):
			continue
		interval = int((stop_n-st)/60)
		a_s = start if start > st else st
		a_e = stop if stop < stop_n else stop_n
		if (a_e > a_s) and (stop_n > st) and interval > 0:
			for k in range(a_s, a_e, interval) : 
				if strand_nm == '+' :
					ax_bottom.arrow(k, 1, interval, 0, head_width=0.02, head_length=interval/2, overhang=1)
				else:
					ax_bottom.arrow(k, 1, interval*(-1),0,head_width=0.02, head_length=interval/2, overhang=1)


		# 엑손을 그린다.
		print('exons')

		for j, e_s in enumerate(exon_s_nm) :
			for k, es in enumerate(ess_nm) :
				if (es > stop_n or ees_nm[k] < st) :
					continue
				rect = patches.Rectangle((es, 0.8-j),ees_nm[k]-es,0.4,edgecolor='none',facecolor='black')
				ax_bottom.add_patch(rect)
				leftt = es if es > st else st
				rightt = ees_nm[k] if ees_nm[k] < stop_n else stop_n
				ax_bottom.text((leftt+rightt)/2, 1-j, str(k+1), horizontalalignment='center', verticalalignment='center', color='white')





		## Genetic variants legend 표시##
		var_legends = []
		var_names = []

		# 각 variant의 모양과 이름 수집 후 legend에 저장
		for ef_i, ef in enumerate(effect_list):
			var_mark = var_markers[ef_i%len(var_markers)]
			var_col = var_colors[ef_i%len(var_colors)]
			var_names.append(ef)
			var_legends.append(ax_main.scatter([], [], marker=var_mark, color=var_col, edgecolor='black', s=marker_size**2, linewidth='2', alpha=0.5))
		
		# 저장 된 것을 바탕으로 legend 표시
		ax_main.legend(tuple(var_legends), tuple(var_names), loc='upper right', title='Generic Variants')

		# lineplot 간 여백 설정
		plt.subplots_adjust(top = 0.95, bottom = 0.05, right = 0.9, left = 0.1, wspace=0, hspace=0)





		## pdf로 저장 ##
		plt.savefig(output_prefix+'_'+str(n+1)+'.pdf')
		plt.close(fig2)

		print(output_prefix+'_'+str(n+1)+'.pdf saved!')








