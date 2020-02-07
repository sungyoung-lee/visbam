Visbam
======

.. meta::
   :description lang=ko: Visualize coverages of multiple bam files

Visbam은 여러 Bam file들의 coverage를 계산하여 visualize하여 줍니다.
주요 기능은 다음과 같습니다.

Coverage Visualize
    여러 Bam파일들의 Coverage를 계산 후 하나의 Line plot으로 표현합니다.
    계산 후 결과를 Caching하여 한번 계산한 후에는
    더욱 빠르게 결과를 볼 수 있습니다.

Refseq Visualize
    Bam파일의 Coverage가 선택한 DNA의 어느 위치에 있는지 알기 쉽도록
    Coverage 계산 결과 하단에 Refseq Data를 토대로
    Exon정보 등을 표시합니다.

Sample Clustering
    Bam파일들을 Coverage의 계산 결과에 따라 두 그룹으로 나누어 줍니다.
    각 그룹은 붉은색과 초록색으로 표시됩니다.
    알고리즘을 선택하고 설정을 바꾸어 최적의 결과를 표시할 수 있습니다.
    현재 3가지 알고리즘을 지원하고 있습니다.
