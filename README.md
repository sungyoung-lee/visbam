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







## visualize_not_normal_gdc.py

gdc samples에서의 clustering과 exon 사이의 값을 확인하기 위한 프로그램입니다. bam 파일을 normal bam 파일 없이 그립니다. 그리고 bam 파일을 K-Means clustering을 통해 두 그룹으로 나눈 후 그립니다. combine 모드를 켜면 모든 exon부위들을 한번에 모아 한 파일로 출력할 수 있습니다. 또한 view_mode에서 2번 선택 시 선택된 exon 사이의 값이 많이 차이나는 sample을 강조할 수 있는 기능을 추가했습니다. 그리고 그러한 sample과 그렇지 않은 sample을 각각 positive(1), negative(0)로 정의하여 logistic regression을 하여 그 결과를 (output_prefix)_summary.txt 로 저장하는 기능을 추가하였습니다.

### usage

visualize_not_normal.py [-h] [-f] [--flag3] [--exclude EXCLUDE] [--combine] [--view_mode VIEW_MODE] [--select SELECT] [--train TRAIN] [--test TEST] [--threshold THRESHOLD] bam_dir_path refseq_path nmid_to_draw draw_span output_prefix

### positional arguments
#### bam_dir_path          
bam파일을 읽어들일 디렉토리를 정합니다.
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
#### -f, --flag            
exon 주변의 100bp 부분만 그립니다.
#### --exclude EXCLUDE              
curated genes만 그립니다.
#### --combine            
모든 exon 부위들을 붙여서 그립니다.
#### --view_mode VIEW_MODE               
##### view_mode 0
아무 표시도 하지 않고 그립니다.
##### view_mode 1
각 region의 평균을 계산하여 한 붉은 선으로 표시합니다.
##### view_mode 2
이 모드는 현재 combine 속성이 켜 있을 때만 적용됩니다. 두 엑손을 선택하여 두 엑손과 두 엑손의 각각 다음 엑손 사이의 30bp 평균이 모두 95% 오차 범위를 넘어가는 sample들만 붉은 선으로 표시합니다. 두 엑손의 선택은 select 속성에서 합니다.
#### --select SELECT num1,num2
combine 속성이 켜져 있고 view_mode가 2일때만 적용합니다. 두 엑손의 번호 num1과 num2를 선택합니다. (1부터 시작, 자연수)
#### --train TRAIN num1              
view_mode 2에서는 표시된 sample과 그렇지 않은 sample을 구분하는 logistic regression을 만듭니다. 여기서 train할 sample의 범위를 지정하기 위하여 특정 exon 번호를 지정해 줍니다. 단일 exon만 가능합니다.
#### --test TEST num2              
logistic regression을 test할 exon 번호를 지정해 줍니다. 단일 exon만 가능합니다.
#### --threshold THRESHOLD               
combine 속성이 켜져 있고 view_mode가 2일때만 적용합니다. view_mode 2 에서 표시된 sample중 특정 값 이상을 넘어가는 sample을 제외합니다.
