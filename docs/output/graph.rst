Final graph
==========

최종 결과물 그래프에 대해서 설명을 하는 페이지입니다.

.. figure:: ../img/combined_graph.png
    :align: center
    :figwidth: 100%
    :target: ../img/combined_graph.png


`\\-\\-combine_slices`_ 옵션이 켜져 있을 시 하나의 PDF로 출력됩니다.


.. figure:: ../img/sliced_graph.png
    :align: center
    :figwidth: 100%
    :target: ../img/sliced_graph.png

그렇지 않으면 각 slice별로 출력됩니다.
각 슬라이스의 파일 명 뒤에는 slice 된 번호
(혹은 exon 번호)가 붙습니다.

그래프는 크게 coverage 부분, RefSeq 부분으로 나누어집니다.
coverage 부분에서는 sample별로 계산한 coverage들을 표시합니다.
RefSeq부분에서는 불러온 RefSeq을 전부 표시합니다.

.. _\\-\\-combine_slices : https://visbam.readthedocs.io/en/latest/input/optional.html#combine-slices

Coverage
--------

.. figure:: ../img/coverage.png
    :align: center
    :figwidth: 100%
    :target: ../img/coverage.png

최종적으로 계산한 coverage를 표시하는 부분입니다.
여기에 `Genetic variants`_ 도 함께 표시됩니다.

:guilabel:`Coverage` 는 line plot으로 그려집니다.
``Clustering`` 을 하지 않을 시 모든 그래프가 반투명한 초록색으로 그려집니다.
``Clustering`` 을 할 시에는 그래프가 빨간색, 초록색 두 그룹으로 나뉘어져 그려집니다.

``Clustering`` 을 할 시 각 slice별로 좌측 상단에 :guilabel:`Clustering Evaluation` 그래프를 표시합니다.
:guilabel:`Clustering Evaluation` 은 silhouette score 기반으로 그려지며 값은 각각 -1 이상입니다.
또 두 그룹으로 나누어서 해당하는 group의 색으로 표시하여 줍니다.

Coverage line plot 위에 :guilabel:`Genetic variants` 를 표시합니다.
각 :guilabel:`Genetic variants` 의 ``Effect`` 별로 모양이 달라집니다.
모양이 어떤 ``Effect`` 를 의미하는 지는 우측 상단에 :guilabel:`Legend` 를 통해서 확인 할 수 있습니다.
또 색은 :guilabel:`Genetic variants` 의 Bam이 어떤 group에 ``Clustering`` 되었는지에 따라 달라집니다.

그리고 그래프 하단에 :guilabel:`Genetic variants` 가 존재하는 position에 한해
position별로 :guilabel:`Pie Graph` 가 표시됩니다.
어떤 group의 Bam의 :guilabel:`Genetic variants` 가 그 position에 있는지 한 눈에 확인할 수 있습니다.

.. _`Genetic variants` : https://en.wikipedia.org/wiki/Genetic_variant


RefSeq
------

.. figure:: ../img/refseq.png
    :align: center
    :figwidth: 100%
    :target: ../img/refseq.png

선택된 RefSeq들을 coverage 부분의 position에 따라 표시해 줍니다.
tx, cds와 direction, exon을 그려줍니다.

tx, cds
~~~~~~~

tx, cds 부분은 RefSeq의 시작과 끝을 표시합니다.

direction
~~~~~~~~~

Direction은 해당 RefSeq의 전사 방향을 표시하여 줍니다.
left면 왼쪽, right면 오른쪽으로 화살표 모양으로 표시됩니다.

exon
~~~~

Exon은 넓은 직사각형 모양으로 표시됩니다.
Exon의 중앙에 각 exon의 번호를 표시합니다.
