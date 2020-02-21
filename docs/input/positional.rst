Required arguments
==================

Visbam을 실행할 때 필수로 입력해야 하는 필드들입니다.
이 필드들을 입력하지 않으면 실행이 되지 않습니다.


bam_path
--------

.. code::

   --bam_path <bam_path>

Bam_ 파일들이 있는 경로를 지정하여 줍니다.
디렉토리 내의 폴더 안에 있는 Bam_ 파일들은 검색하지 않습니다.
이 디렉토리에 있는 Bam_ 파일들중
sample_path_ 파일의 목록에 있는 Bam_ 파일들만 coverage가 계산되며
그래프로 결과가 출력됩니다. 

.. _Bam : https://en.wikipedia.org/wiki/Binary_Alignment_Map

sample_path 
-----------

.. code::

   --sample_path <sample_path>

bam_path_ 중 coverage를 계산할 Bam_ 파일 목록을
저장한 파일의 경로를 지정합니다.
sample_path 파일에는 각 Bam_ 파일명의
'.'앞까지의 이름이 엔터로 구분되어 있어야 합니다.

.. warning::
    파일명은 맨 앞의 '.'이 위치한 바로 앞부분까지만 적어야 합니다.
   
    예시 :

    .. code::
        
       MS190000066_S8.bwamem.sorted.dedup.realn.recal.dedup.bam
       
       -> MS190000066_S8

.. _Bam : https://en.wikipedia.org/wiki/Binary_Alignment_Map

normal_path 
-----------

.. code::

   --normal_path <normal_path>

그래프를 그릴 Bam_ 파일을 normalize해줄 Normal Bam_ 파일들의 경로입니다.
디렉토리 내의 폴더 안에 있는 Bam_ 파일들은 검색하지 않습니다.
이 Bam_ 파일들은 sample_path_ Bam_ 파일의 coverage를 normalize 할 때 사용됩니다.
따라서 이 디렉토리에 있는 Bam_ 파일들은 그래프로 출력되지 않습니다.
normalize에 관한 내용은 Reading_Files_ 문서를 참조하십시오.

.. _Reading_Files: https://visbam.readthedocs.io/en/latest/process/read_files.html#normal-bam

.. _Bam : https://en.wikipedia.org/wiki/Binary_Alignment_Map

refseq_path
-----------

.. code::

   --refseq_path <refseq_path>

RefSeq_ 데이터를 불러옵니다.
RefSeq_ 데이터는 `TSV(Tab-Separated Values)`_ 파일 형식이어야 합니다.
아래와 같은 column이 포함되어 있어야 합니다.

.. code::

   'bin', 'name', 'chrom', 'strand',
   'txStart', 'txEnd', 'cdsStart', 'cdsEnd',
   'exonCount', 'exonStarts', 'exonEnds', 'score',
   'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames'

.. _RefSeq : https://en.wikipedia.org/wiki/RefSeq
.. _`TSV(Tab-Separated Values)` : https://en.wikipedia.org/wiki/Tab-separated_values

또, NR을 제외한 NM만 불러옵니다. 



variant_path
------------

.. code::

   --variant_path <variant_path>

`Genetic variants`_ 의 데이터가 있는 경로를 설정해줍니다.
각 파일의 파일명의 시작은 Bam_ 파일의 파일명('.'앞까지)이어야 합니다.
`Genetic variants`_ 는 `TSV(Tab-Seperated Values)`_ 형식의 txt 파일로
열에 ``RefSeq``, ``Pos``, ``Effect`` 가 포함되어 있어야 합니다.

.. _`Genetic variants` : https://en.wikipedia.org/wiki/Genetic_variant
.. _Bam : https://en.wikipedia.org/wiki/Binary_Alignment_Map
.. _`TSV(Tab-Separated Values)` : https://en.wikipedia.org/wiki/Tab-separated_values


refseq
------

.. code::

   --refseq <refseq>

그래프를 그릴 RefSeq_ 의 NMID를 적어야 합니다.
이 NMID의 시작과 끝이 coverage를 추출할 범위입니다.
NR을 제외한 NM만 가능합니다.

.. warning::
    뒤에 버전(. 뒷부분)은 떼고 적습니다.
   
    예시 :

    .. code::
      
       NM_001005484.1  
       
       -> NM_001005484


.. _RefSeq : https://en.wikipedia.org/wiki/RefSeq


prefix
------

.. code::

   --prefix <prefix>

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
