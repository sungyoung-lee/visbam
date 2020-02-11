Clustering Samples
==================

Sample을 두 그룹으로 Clustering 하는 Process입니다.
Clustering 알고리즘은 3가지가 있습니다.
``silhouette``, ``nmf``, ``splice_site`` 입니다.

Silhouette
----------

NMF
---

NMF를 이용하여 Exon의 시작과 끝 앞뒤로 ``interval`` 을 둔 범위의
Coverage 값을 모든 Sample에 대해 불러옵니다.
이 값을 NMF(Non-nefative matrix factorization)을 이용하여 W와 H 행렬로 나눕니다.
W와 H행렬의 각각 행과 열은 2로 나누어 줍니다.
이때 W행은 Sample 만큼의 행과 2열로 구성됩니다.
해당 Sample의 행에서 0번째 열보다 1번째 열이 더 클 경우 Clustering합니다.


Splice_Site
-----------

Generic Variants 데이터 중 Effect의 이름이 ``splice_site`` 가
있는 것들을 모아 줍니다.
Splice_Site가 있는 Sample을 모두 선택하여 Clustering 하여 줍니다.
