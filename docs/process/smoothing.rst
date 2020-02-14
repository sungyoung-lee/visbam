Smoothing the Graph
===================

그래프를 더 간단히 볼 수 있도록 Smoothing 하는 과정입니다.
Smoothing은 두 가지 알고리즘을 선택하여 진행합니다.
``average`` 와 ``loess`` 입니다.


Average
-------

Average는 각 position 앞 뒤로 ``average`` 간격 만큼의 평균을
계산하여 그 값을 구하는 방식으로 Smoothing을 합니다.
이 때 해당 position 앞뒤로 ``average`` 만큼의 값이 없는
position의 경우에는 그 값이 0으로 처리됩니다.
`--fill`_ 옵션을 켜 Exon 앞뒤의 값을 더 불러온 뒤 처리할 수 있습니다.

.. _--fill : https://visbam.readthedocs.io/en/latest/input/optional.html#fill

Loess
------

Loess는 skmisc의 loess라이브러리를 이용하여 Smoothing하는 방식입니다.
자세한 방식은 skmisc_ 를 참조하십시오.

.. _skmisc: https://has2k1.github.io/scikit-misc/loess.html
