Positional Arguments
=====================

Visbam을 실행할 때 필수로 입력해야 하는 필드들입니다.
이 필드들을 입력하지 않으면 실행이 되지 않습니다.
그리고 반드시 아래와 같은 순서로 입력하셔야 합니다.

python3 visualize.py ``bam_dir_path`` ``sample_list_path`` ``normal_dir_path``
``refseq_path`` ``draw_span`` ``output_prefix``


bam_dir_path
------------

Bam파일들이 있는 경로를 지정하여 줍니다.
디렉토리 내의 폴더 안에 있는 Bam파일들은 검색하지 않습니다.
이 디렉토리에 있는 Bam파일들중
``sample_list_path`` 파일의 목록에 있는 Bam파일들만 Coverage가 계산되며
그래프로 결과가 출력됩니다. 


sample_list_path 
----------------

``bam_dir_path`` 중 Coverage를 계산할 Bam파일 목록을
저장한 파일의 경로를 지정합니다.

.. warning::
    파일명은 맨 앞의 '.'이 위치한 바로 앞부분까지만 적어야 합니다.
   
    예시 :

    .. code::
        
       MS190000066_S8.bwamem.sorted.dedup.realn.recal.dedup.bam
       
       -> MS190000066_S8


normal_dir_path 
---------------

그래프를 그릴 Bam 파일을 Normalize해줄 Normal Bam 파일들의 경로입니다.
디렉토리 내의 폴더 안에 있는 Bam 파일들은 검색하지 않습니다.
이 Bam 파일들은 ``sample_list_path`` Bam 파일의 Coverage를 Normalize 할 때 사용됩니다.
따라서 이 디렉토리에 있는 Bam파일들은 그래프로 출력되지 않습니다.
Normalize에 관한 내용은 Reading_Files_ 문서를 참조하십시오.

.. _Reading_Files: https://visbam.readthedocs.io/en/latest/process/read_files.html#normal-bam


refseq_path
-----------

Refseq데이터를 불러옵니다.
Refseq데이터는 Tsv(Tab-Separated Values)파일 형식이어야 합니다.
NR을 제외한 NM만 불러옵니다. 


nmid_to_draw
------------

그래프를 그릴 Refseq의 NMID를 적어야 합니다.
이 NMID의 시작과 끝이 Coverage를 추출할 범위입니다.
NR을 제외한 NM만 가능합니다.

.. warning::
    뒤에 버전(. 뒷부분)은 떼고 적습니다.
   
    예시 :

    .. code::
      
       NM_001005484.1  
       
       -> NM_001005484


draw_span
----------

.. note::
    곧 Optional로 수정될 Argument입니다.

옵션 중 exon_sliced_ 가 활성화 되지 않았을 때,
그래프를 나눌 bp의 크기를 정합니다.

.. _exon_sliced: optional.html#exon-sliced


output_prefix
-------------

output을 출력할 경로와 파일명을 지정합니다.


.. warning::
    코드가 존재하는 경로를 기준으로
    확장자를 제외한
    파일 경로와 파일명까지만 적습니다.

    예시 :

    .. code::
      
       코드 경로가 user/unknown/code/visualize.py 이고,

       저장할 경로가 user/unknown/pdf/coverage.pdf 이면

       ../pdf/coverage 를 적으면 됩니다.
