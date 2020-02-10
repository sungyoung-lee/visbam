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

예 :

.. code::

    --exclude_exon 1,2,3


exon_sliced
~~~~~~~~~~~

이 옵션을 켤 시, exon별로 그래프를 그릴 구간을 나누게 됩니다.
그렇지 않으면, draw_span_ 에 따라 나누게 됩니다. 
각 그래프의 차이는 (링크)에서 확인할 수 있습니다.

.. _draw_span: positional.html#draw_span


exon_space
~~~~~~~~~~~

.. note::

    이 옵션을 이용하려면 ``exon_sliced`` 옵션이 활성화되어야 합니다.

``exon_sliced`` 옵션에서 exon 앞뒤의 간격을 bp단위로 설정하여 줍니다.


font_size
~~~~~~~~~~

폰트 크기를 설정합니다. 단위는 pt입니다.


marker_size
~~~~~~~~~~~

Generic Variants를 Visualize할 때 Marker의 크기를 조정합니다.
단위는 pt입니다.


ylim
~~~~

그래프를 표시할 Coverage의 최댓값을 설정합니다.
이 옵션이 없으면 모든 Sample의 Coverage 중
제일 높을 값으로 설정됩니다.




Smoothing
---------

그래프를 Smoothing하는 것과 관련된 Option들입니다.


smoothing
~~~~~~~~~~

어떤 속성으로 Smoothing을 할 지 설정합니다.
설정할 수 있는 속성은 다음과 같습니다.


* ``average``

* ``loess``

Smoothing 속성에 대한 자세한 정보는 (Process 링크)를 참조하십시오.


average
~~~~~~~~

.. note::

    이 옵션을 이용하려면 ``smoothing`` 옵션이 ``average`` 이어야 합니다.

Smoothing이 average일 때, average를 적용할 앞 뒤 bp간격을 설정합니다.
average가 1이면, 앞과 뒤에 각각 1bp가 적용되어 계산됩니다.

fill
~~~~~

.. note::

    이 옵션을 이용하려면 ``smoothing`` 옵션이 ``average`` 이어야 합니다.

Smoothing이 average일 때, 앞 뒤로 ``average`` 만큼 간격을 더 주어
그 간격에서 Moving Average를 계산합니다.






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


