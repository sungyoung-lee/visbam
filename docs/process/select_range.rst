Selecting Range to Display
==========================

그래프가 계산되고 보여질 범위를 설정합니다.
이 과정에서는 ``exon_sliced`` 옵션에 따라 크게 결정됩니다.

exon_sliced on
--------------

Exon으로 나눌 시, Exon 앞뒤 ``exon_space`` 값을 더한 범위가 선택됩니다.
첫 Exon과 끝 Exon은 각각 앞, 뒤로 100bp씩 더 추가되어 선택됩니다.


exon_sliced off
---------------

Span으로 따로 설정하여 나눌 시 ``draw_span`` 값에 따라 나눕니다.


