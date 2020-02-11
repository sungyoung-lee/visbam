Result of Clustering
====================

Clustering 과정에서 출력된 결과물에 대해서 설명합니다.
현재 Clustering 과정 중 Silhouette 과정만 별도의 출력을 하고 있습니다.


CI/Tau/Score Scatter Plot
-------------------------

CI(Confidence Interval), Tau에 따른 Score의 관계를 표시한 그래프 입니다.
Scatter Plot으로 출력되었고, X축은 Ci, Y축은 Tau입니다.

Plot의 크기가 Score의 크기 입니다.
우측 상단에 있는 Legend로 Score가 어느 정도 크기인지 표시하였습니다.
또, Plot의 색 또한 Score의 크기에 따라 다릅니다.
오른쪽 Color Bar를 통해 Score의 크기를 알 수 있습니다.
선택된 부분에 대한 Plot에는 빨간색 테두리로 강조 표시를 하였습니다.

이 그래프의 Width/Height는
``score_plot_width``, ``score_plot_height`` 로 지정할 수 있습니다.


Ratio/Score Scatter Plot
------------------------

Clustering 된 것과 그렇지 않은 것 사이의 비율(Ratio)를 X축으로 하고,
Score를 Y축으로 하여 표시한 그래프 입니다.


