visbam
======
여러 BAM파일들을 시각화하는 python 코드입니다.

# 사용법

## visualize_btp5.py

bam 파일을 visualize 하는 기본 코드입니다.

### usage

visualize_btp5.py [-h] [--color COLOR] [--show_mean] [--draw_mode DRAW_MODE] [--simplify] [-f] [--flag2] [--flag3] bam_dir_path sample_list_path normal_dir_path refseq_path nmid_to_draw draw_span output_prefix

### positional arguments
#### bam_dir_path          
bam파일을 읽어들일 디렉토리를 정합니다.
#### sample_list_path      
해당하는 sample 이름이 들어있는 경로를 지정합니다.
#### normal_dir_path       
normal sample이 들어있는 경로를 지정합니다.
#### refseq_path           
refseq 파일의 경로를 지정합니다.
#### nmid_to_draw          
사용할 NMID를 지정합니다.
#### draw_span             
사진을 몇 bp단위로 분할할 것인지 정합니다.
#### output_prefix         
output 파일명을 정합니다.

### optional arguments
#### -h, --help            
show this help message and exit
#### --color COLOR         
선 색을 정합니다.(파일로 정의)
#### --show_mean           
주어진 sample들의 평균을 그립니다.
#### --draw_mode DRAW_MODE
normal로 standardization시 어떤 값을 사용할지 정합니다. 0=평균, 1=최솟값, 2=최댓값(기본값 0)
#### --simplify            
sample 전체를 그리지 않고 최댓값, 최솟값만 간단하게 표시합니다.
#### -f, --flag            
exon 주변의 100bp 부분만 그립니다.
#### --flag2               
최댓값으로 나누어 그립니다.
#### --flag3               
curated genes만 그립니다.

## visualize_clustering.py

bam 파일을 K-Means clustering을 통해 두 그룹으로 나눈 후 그립니다. 사용법은 visualize_btp5.py 와 같습니다.

## visualize_not_normal.py

bam 파일을 normal bam 파일 없이 그립니다. 그리고 bam 파일을 K-Means clustering을 통해 두 그룹으로 나눈 후 그립니다.

### usage

visualize_gdc.py [-h] [--color COLOR] [--show_mean] [--simplify] [-f] [--flag3] bam_dir_path sample_list_path refseq_path nmid_to_draw draw_span output_prefix

### positional arguments
#### bam_dir_path          
bam파일을 읽어들일 디렉토리를 정합니다.
#### sample_list_path      
해당하는 sample 이름이 들어있는 경로를 지정합니다.
#### refseq_path           
refseq 파일의 경로를 지정합니다.
#### nmid_to_draw          
사용할 NMID를 지정합니다.
#### draw_span             
사진을 몇 bp단위로 분할할 것인지 정합니다.
#### output_prefix         
output 파일명을 정합니다.

### optional arguments
#### -h, --help            
show this help message and exit
#### --color COLOR         
선 색을 정합니다.(파일로 정의)
#### --show_mean           
주어진 sample들의 평균을 그립니다.
#### --simplify            
sample 전체를 그리지 않고 최댓값, 최솟값만 간단하게 표시합니다.
#### -f, --flag            
exon 주변의 100bp 부분만 그립니다.
#### --flag3               
curated genes만 그립니다.
