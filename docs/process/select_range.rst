Selecting Range to Display
==========================

그래프가 계산되고 보여질 범위를 설정합니다.
이 과정에서는 `--exon_sliced`_ 옵션에 따라 크게 결정됩니다.

.. _--exon_sliced : https://visbam.readthedocs.io/en/latest/input/optional.html#exon-sliced

Exon_Sliced
-----------

Exon으로 나눌 시, Exon 앞뒤 `--exon_space`_ 값을 더한 범위가 선택됩니다.
첫 Exon과 끝 Exon은 각각 앞, 뒤로 100bp씩 더 추가되어 선택됩니다.

.. _--exon_space : https://visbam.readthedocs.io/en/latest/input/optional.html#exon-space


not Exon_sliced
-------------

Span으로 따로 설정하여 나눌 시 `--draw_span`_ 값에 따라 나눕니다.

.. _--draw_span : https://visbam.readthedocs.io/en/latest/input/optional.html#draw-span
