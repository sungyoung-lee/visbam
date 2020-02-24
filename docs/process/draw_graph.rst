Drawing the final graph
=======================

최종 Final graph를 출력하는 과정입니다.
이전 과정에서 계산했던 smoothing 된(또는 그렇지 않은) coverage,
`Genetic variants`_ , RefSeq_ 값을 이용하여 그래프에 표시합니다.

그래프는 `\\-\\-combine_slices`_  값의 활성화 여부에 따라 한 파일로 출력될 지
여러 파일로 출력될 지 결정됩니다.

Output 그림에 관한 자세한 설명은 Outputs_ 문서를 참조하십시오.


.. figure:: ../img/draw_graph.png
    :align: center
    :figwidth: 100%

.. _Outputs: https://visbam.readthedocs.io/en/latest/output/graph.html

.. _`\\-\\-combine_slices` : https://visbam.readthedocs.io/en/latest/input/optional.html#combine-slices
.. _`Genetic variants` : https://en.wikipedia.org/wiki/Genetic_variant
.. _RefSeq : https://en.wikipedia.org/wiki/RefSeq


Coverages
---------

먼저 sample들을 Line plot으로 표시하여 줍니다.
Clustering이 진행되지 않은 경우 모든 sample이 초록색으로 표시됩니다.
Clustering이 진행 된 경우, clustering 된 sample은 빨간색으로 표시됩니다.

그리고 해당 clustering이 얼마나 잘 되었는지 확인할 수 있는
Clustering Evaluation 그래프를 각 exon별로 그려줍니다.



Genetic variants
----------------

`\\-\\-variant_path`_ 로 불러왔던 `Genetic variants`_ 를
일치하는 sample과 position에 표시하여 줍니다.
그리고 effect 별로 다른 모양을 표시하여 줍니다.


또 그래프 하단에 variant가 존재하는 position에 한해
Clustering 된 두 그룹의 비율을 pie marker로 표시하여 줍니다.


.. _`\\-\\-variant_path` : https://visbam.readthedocs.io/en/latest/input/positional.html#variant-path 
.. _`Genetic variants` : https://en.wikipedia.org/wiki/Genetic_variant

RefSeq
-------

그래프 하단에 RefSeq_ 를 표시하여줍니다.
Reading files과정에서 불러왔던 RefSeq_ 들을 모두 아래쪽 하단에 표시합니다.

.. _RefSeq : https://en.wikipedia.org/wiki/RefSeq
