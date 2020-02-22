Visbam
======

.. meta::
   :description lang=ko: Visualize coverages of multiple bam files

Visbam은 여러 Bam_ file들의 coverage를 계산하여 visualize하여 주는 프로그램입니다.
주요 기능은 다음과 같습니다.

Coverage Visualize
    여러 Bam_ 파일들의 coverage를 계산 후 하나의 line plot으로 표현합니다.
    계산 후 결과를 caching하여 한번 계산한 후에는
    더욱 빠르게 visualize 할 수 있습니다.

RefSeq Visualize
    Bam_ 파일의 coverage가 선택한 dna의 어느 위치에 있는지 알기 쉽도록
    coverage 계산 결과 하단에 RefSeq_ data를 토대로
    exon정보 등을 표시합니다.

Sample Clustering
    Bam_ 파일들을 coverage의 계산 결과에 따라 두 그룹으로 나누어 줍니다.
    각 그룹은 붉은색과 초록색으로 표시됩니다.
    Coverage 하는 알고리즘을 선택하고 설정을 바꾸어,
    최적의 결과를 표시할 수 있도록 하였습니다.
    현재 3가지 알고리즘을 지원하고 있습니다.

.. _Bam : https://en.wikipedia.org/wiki/Binary_Alignment_Map
.. _RefSeq : https://en.wikipedia.org/wiki/RefSeq

Prerequisites
------------

Visbam을 실행하기 전 설치되어야 할 파이썬 패키지들이 있습니다.
아래 패키지들을 설치한 후 Visbam을 실행할 수 있습니다.
Python 3 버전만 지원하며, 3.5.2 이상 버전을 권장합니다.

* NumPy_ : 여러 수학적 계산을 위해 필요합니다.
* pandas_ : Coverage 데이터 분석을 위해 필요합니다.
* Matplotlib_ : 그래프를 그리기 위해 사용됩니다.
* pysam_ : Bam 파일을 읽어들여 coverage를 구하기 위해 사용합니다.
* scikit-learn_ : NMF clutering과 clustering 과정에서 Silhouette score를 계산하기 위해 사용됩니다. 
* scikit-misc_ : Loess smoothing 을 위해 필요합니다.


.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Prerequisites

.. _NumPy : https://numpy.org/
.. _pandas : https://pandas.pydata.org/
.. _Matplotlib : https://matplotlib.org/
.. _pysam : https://pysam.readthedocs.io/en/latest/index.html
.. _scikit-learn : https://scikit-learn.org/stable/
.. _scikit-misc : https://has2k1.github.io/scikit-misc/installation.html





Inputs
------

Visbam을 실행시키려면 여러 input이 필요합니다.
필수로 넣어야 하는 `Required arguments`_ 와 선택 옵션인 `Optional arguments`_ 이 있습니다.
`Optional arguments`_ 중에서도 특정 option을 입력할 시
필수로 입력해야 하는 필드가 있으니 주의하시기 바랍니다.
전체 명령어의 목록은 다음과 같습니다.
더 자세한 설명은 아래 문서를 참고하십시오.

* :doc:`Required arguments <input/positional>`
* :doc:`Optional arguments <input/optional>`



.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Inputs

   input/positional
   input/optional

.. _`Required arguments` : input/positional
.. _`Optional arguments` : input/optional


Process
--------

Visbam의 실행 프로세스와 알고리즘을 적어 놓은 문서입니다.
전체적인 과정을 정리하면 아래와 같습니다.

.. figure:: img/visbam_whole_process.png
    :align: center
    :figwidth: 100%
    :target: img/visbam_whole_process.png

    전체적인 Visbam Process 요약도

각 단계별로 개별적으로 문서를 정리하였습니다.

* :doc:`Reading files <process/read_files>`
* :doc:`Select range to display <process/select_range>`
* :doc:`Smoothing the graph <process/smoothing>`
* :doc:`Clustering samples <process/clustering>`
* :doc:`Drawing the final graph <process/draw_graph>`

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Process

   process/read_files
   process/select_range
   process/smoothing
   process/clustering
   process/draw_graph


Outputs
--------

Visbam의 출력 결과를 정리한 문서입니다.
clustering 과정에서 제공하는 중간 결과 그래프와
최종적으로 그려지는 coverage의 line plot이 있습니다.
단계별로 나올 수 있는 output들을 따로 정리하였습니다.

* :doc:`Result of clustering <output/clustering>`
* :doc:`Final graph <output/graph>`

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Outputs

   output/clustering
   output/graph


References
----------

전체 references를 정리해 놓은 페이지입니다.

* :doc:`References <references>`

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: References

   references
