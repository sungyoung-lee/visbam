Optional Arguments
==================

Optional Arguments를 정리해 놓은 문서입니다.
Visbam에 필요한 유용한 기능들이 많이 있습니다.
대부분 필수적으로 적어야 하지 않지만, 특정 Option에 따라
설정해줘야 하는 Input들도 있으니 참조하시기 바랍니다.
각 Option이 담당하는 기능별로 모았습니다.


Drawing Graph
-------------

그래프를 그리는 것과 관련된 Option들입니다.


combine_slices
~~~~~~~~~~~~~~

엑손별로, 혹은 각 bp별로 나눠진 graph들을 합쳐주는 옵션입니다.
합쳤을 때와 그렇지 않을 떄의 차이는 (링크)에서 확인할 수 있습니다.


curated_genes
~~~~~~~~~~~~~~

.. note::

    이 옵션을 이용하려면 소스 코드와 같은 폴더 내에
    tsv형식으로 되어 있는 ``allCuratedGenes.txt`` 파일이 있어야 합니다.
    ``allCuratedGenes.txt`` 에 대해서는 (링크)를 참조하시기 바랍니다.

Refseq 데이터 중 allCuratedGenes.txt에 포함되어 있는 Refseq 데이터만 표시합니다.


draw_average_line
~~~~~~~~~~~~~~~~~

전체 샘플의 bp별 평균을 Line Plot으로 표시합니다.
(그림)

exclude_exon
~~~~~~~~~~~~~


.. note::

    이 옵션을 이용하려면 ``exon_sliced`` 옵션이 활성화되어야 합니다.

일부 엑손을 제외하고 표시합니다.
엑손을 여러개를 선택하려면 쉼표로 구분하여 표시합니다.
예 : --exclude_exon 1,2,3


exon_sliced
~~~~~~~~~~~

exon_space
~~~~~~~~~~~

font_size
~~~~~~~~~~

marker_size
~~~~~~~~~~~

ylim
~~~~






Smoothing
---------

그래프를 Smoothing하는 것과 관련된 Option들입니다.


average
~~~~~~~~

fill
~~~~~

smoothing
~~~~~~~~~~





Clustering
---------

Sample들을 Clustering하는 것과 관련된 Option들입니다.


clustering
~~~~~~~~~~

clustering_mode 
~~~~~~~~~~~~~~~

select_exon
~~~~~~~~~~~

threshold
~~~~~~~~~~

score_plot_width
~~~~~~~~~~~~~~~~

score_plot_height
~~~~~~~~~~~~~~~~~

limit_tau
~~~~~~~~~~

limit_tau_low
~~~~~~~~~~~~~

silhouette_dintv
~~~~~~~~~~~~~~~~


