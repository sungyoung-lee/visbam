Reading Files
==============

Visbam 프로그램에서 제일 처음에 진행되는 작업입니다.
본격적인 visualize system을 시작하기 전, 파일을 불러와 초기화 시켜줍니다.
불러오는 파일은 Bam_ 파일, `Genetic Variants`_ 데이터, RefSeq_ 파일이 있습니다.

.. _Bam : https://en.wikipedia.org/wiki/Binary_Alignment_Map
.. _RefSeq : https://en.wikipedia.org/wiki/RefSeq
.. _`Genetic Variants` : https://en.wikipedia.org/wiki/Genetic_variant

BAM Files
---------

먼저 pysam 라이브러리를 이용하여 Bam_ 파일들을 불러와줍니다.
Normal Bam_ 파일부터 불러온 뒤, visualize 할 Bam_ 들을 불러와줍니다.

.. _Bam : https://en.wikipedia.org/wiki/Binary_Alignment_Map

Normal Bam
~~~~~~~~~~

normal_dir_path_ 내에 있는 Bam_ 파일들을 모두 불러와줍니다.
Bam_ 파일들을 하나씩 불러 온 뒤 해당 파일과 구간의 cache 파일이
있는지 없는지 확인합니다.
cache 파일이 존재하면 cache 파일을 불러옵니다.
없으면 각 파일의 coverage를 Pysam_ 라이브러리를 이용하여 구해 줍니다.
그리고 coverage를 caching하여 저장합니다.

Cache file은 해당 코드가 있는 폴더에 ``cache/``
디렉토리를 새로 만든 뒤 그 안에 파일들이 저장됩니다. 

구해진 Normal bam file들의 평균을 최종 결과값으로 저장합니다.


.. _normal_dir_path: https://visbam.readthedocs.io/en/latest/input/positional.html#normal-dir-path

.. _Bam : https://en.wikipedia.org/wiki/Binary_Alignment_Map
.. _Pysam : https://pysam.readthedocs.io/en/latest/index.html

Cancer Bam
~~~~~~~~~~

Coverage 계산까지 Normal Bam_ 과 동일하게 진행됩니다.
그리고 구해진 Cancer Bam_ 에 position별로 각각
해당 positon의 Normal Bam_ 의 평균으로 나누어 줍니다.
이 결과를 최종적으로 전부 저장하여 visualize합니다.

.. _Bam : https://en.wikipedia.org/wiki/Binary_Alignment_Map

Generic Variants
----------------

variants_dir_path_ 폴더에 있는 모든 ``.txt`` 파일들을 불러옵니다.
이때 txt 파일명은 Bam_ file의 파일명 중 온점 앞까지의 값과 같아야 합니다.

예 :

.. code::

    MS180000237_S10.bwamem.sorted.dedup.realn.recal.dedup.bam
 
    -> MS180000237_S10.txt

그 뒤 nmid_to_draw_ 의 nmid와 일치하는 값이며,
visualize 하는 sample에 해당되는 값만 불러옵니다.

.. _variants_dir_path : https://visbam.readthedocs.io/en/latest/input/positional.html#variants-dir-path

.. _nmid_to_draw : https://visbam.readthedocs.io/en/latest/input/positional.html#nmid-to-draw

.. _Bam : https://en.wikipedia.org/wiki/Binary_Alignment_Map

RefSeq
------

refseq_path_ 에 있는 RefSeq_ 파일을 불러옵니다.
이 후 nmid_to_draw_ 와 일치하는 RefSeq_ 중 첫번째 값을 불러옵니다.
이 값이 그래프의 처음과 끝, exon을 결정하는 RefSeq_ 이 됩니다.

그리고 RefSeq_ 중 nmid_to_draw_ 의 범위에 들어오는 모든 RefSeq_ 을 불러옵니다.

`\\-\\-curated_genes`_ 파일이 지정되어 있으면 `\\-\\-curated_genes`_ 파일에 있는 RefSeq_ 만 걸러냅니다.


.. _nmid_to_draw : https://visbam.readthedocs.io/en/latest/input/positional.html#nmid-to-draw

.. _refseq_path : https://visbam.readthedocs.io/en/latest/input/positional.html#refseq-path

.. _\\-\\-curated_genes : https://visbam.readthedocs.io/en/latest/input/optional.html#curated-genes
.. _RefSeq : https://en.wikipedia.org/wiki/RefSeq
