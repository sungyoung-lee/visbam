Reading Files
==============

Visbam 프로그램에서 제일 처음에 진행되는 작업입니다.
본격적인 Visualize System을 시작하기 전, 파일을 불러와 초기화 시켜줍니다.
불러오는 파일은 Bam 파일, Generic Variants 데이터, Refseq 파일이 있습니다.


BAM Files
---------

먼저 pysam 라이브러리를 이용하여 BAM 파일들을 불러와줍니다.
Normal Bam 파일부터 불러온 뒤, Visualize 할 Bam들을 불러와줍니다.

Normal Bam
~~~~~~~~~~

normal_dir_path(링크) 내에 있는 BAM 파일들을 모두 불러와줍니다.
Bam파일들을 하나씩 불러 온 뒤 해당 파일과 구간의 Cache 파일이
있는지 없는지 확인합니다.
Cache 파일이 존재하면 Cache 파일을 불러옵니다.
없으면 각 파일의 Coverage를 pysam라이브러리를 이용하여 구해 줍니다.
그리고 Coverage를 Caching하여 저장합니다.

Cache File은 해당 코드가 있는 폴더에 ``cache/``
디렉토리를 새로 만든 뒤 그 안에 파일들이 저장됩니다. 

구해진 Normal Bam file들의 평균을 최종 결과값으로 저장합니다.


Cancer Bam
~~~~~~~~~~

Coverage 계산까지 Normal Bam과 동일하게 진행됩니다.
그리고 구해진 Cancer Bam에 position별로 각각
해당 positon의 Normal Bam의 평균으로 나누어 줍니다.
이 결과를 최종적으로 전부 저장하여 Visualize합니다.


Generic Variants
----------------

``variants_dir_path`` 폴더에 있는 모든 ``.txt`` 파일들을 불러옵니다.
이때 txt는 Bam File의 파일명 중 온점 앞까지의 값과 같아야 합니다.

예 :

.. code::
    MS180000237_S10.bwamem.sorted.dedup.realn.recal.dedup.bam
 
    -> MS180000237_S10.txt

그 뒤 ``nmid_to_draw`` 의 nmid와 일치하는 값이며,
Visualize 하는 Sample에 해당되는 값만 불러옵니다.


Refseq
------

``refseq_path`` 에 있는 Refseq 파일을 불러옵니다.
이 후 ``nmid_to_draw`` 와 일치하는 Refseq 중 첫번째 값을 불러옵니다.
이 값이 그래프의 처음과 끝, Exon을 결정하는 Refseq이 됩니다.

그리고 refseq 중 ``nmid_to_draw`` 의 범위에 들어오는 모든 Refseq을 불러옵니다.

``--curated_genes`` 파일이 지정되어 있으면 ``--curated_genes`` 파일에 있는 Refseq만 걸러냅니다.
