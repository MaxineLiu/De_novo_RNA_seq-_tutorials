# [转录组组装软件：Trinity(转载)](https://links.jianshu.com/go?to=https%3A%2F%2Fbaijiahao.baidu.com%2Fs%3Fid%3D1597982350547011645%26wfr%3Dspider%26for%3Dpc)  

视频参考：[BroadE: Trinity - How it works](https://links.jianshu.com/go?to=https%3A%2F%2Fwww.youtube.com%2Fwatch%3Fv%3DGccnW_g-4nE%26t%3D176s)和基因课的转录组原理的trinity原理视频类似。

本文添加了一些其他内容便于理解，添加内容有 ★{添加的内容}做标记。

  

![](https://raw.githubusercontent.com/MaxineLiu/picBed/main/img/202111100243380.webp)

  

![](https://raw.githubusercontent.com/MaxineLiu/picBed/main/img/202111100243382.webp)

trinity拼接原理.jpg

Trinity，是由 the Broad Institute 开发的转录组denovo组装软件，由三个独立的软件模块组成：Inchworm, Chrysalis和Butterfly。三个软件依次来处理大规模的RNA-seq的reads数据。

## Trinity的简要工作流程

**Inchworm**  
将RNA-seq的原始reads数据组装成Unique序列；  
★｛unique sequences，也就是contigs｝

**Chrysalis**  
将上一步生成的contigs聚类，然后对每个components构建deBruijn图；

**Butterfly**  
处理这些deBruijn图，依据图中reads和成对的reads来寻找路径，从而得到具有可变剪接的全长转录子，同时将旁系同源基因的转录子分开。

#### Trinity原理详细的介绍

### 1）Inchworm生成contig

* * *

★｛kmer  
这里首先需要知道一个专有名词的概念，mer，其在分子生物学领域中意义为单体单元 （monomeric unit，mer）。通常用于核酸序列中的单位，代表nt或者bp，例如，100 mer DNA代表这段DNA序列单链长度100nt，或者双链长度100bp。

而k-mer则是指将核酸序列分成包含k个碱基的字符串，即从一段连续的核酸序列中迭代地选取长度为K个碱基的序列，若核酸序列长度为L，k-mer长度为K，那么可以得到L-K+1个k-mers。如下图所示，假设这里存在某序列长度为21，设定选取的k-mer长度为7，则得到（21-7+1=15）个7-mers。

![](https://raw.githubusercontent.com/MaxineLiu/picBed/main/img/202111100243383.webp)  

｝

---

假设kmer长度是k，将测序reads以k-1的overlap分割成长度为k的k-mer，去除可能错误的kmer，低复杂度和单一的kmers，构建成kmer库，以出现次数最高的kmer作为基序，基于k-1个overlap向两边贪婪延伸，直到不能延伸之后得到一个contig，去除已经使用过的kmer，对于剩余的kmer按照上述方法延伸得到contig。直到kmer用完，生成contig过程结束。具体过程如下图所示：

举个具体的例子，假如kmer是7，kmer库为：

![kmer_counts](https://raw.githubusercontent.com/MaxineLiu/picBed/main/img/202111100243384.webp)

其中出现次数最高的是GATTACA，以ATTACA为overlap向两端延伸，下一个碱基的可能是只有四种G、A、T、C。其中出现次数最高ATTACAG和ATTACAC（4次）。

![](https://raw.githubusercontent.com/MaxineLiu/picBed/main/img/202111100243385.webp)

![](https://raw.githubusercontent.com/MaxineLiu/picBed/main/img/202111100243386.webp)

分别以TTACAG和TTACAC为基序延伸，以此类推延伸（反向也是如此），直到不能再延伸为止。最终contig序列为：…AAGATTACAGA…

![](https://raw.githubusercontent.com/MaxineLiu/picBed/main/img/202111100243387.webp)

### 2）contig聚类成components

  ![](https://raw.githubusercontent.com/MaxineLiu/picBed/main/img/202111100243388.webp)

根据最小overlap聚类contig。每一个component 是由contigs组成的集合，这些contig可能是来自可变剪切体或者相近的旁系同源物。

contig聚类满足的条件：

（1） contig之间有k-1碱基的overlap

（2） 满足跨越两个contig的junction的最小reads数，且（k-1）mer的junction两端分别有（k-1）/2的碱基支持。

![](https://raw.githubusercontent.com/MaxineLiu/picBed/main/img/202111100243389.webp)

![](https://raw.githubusercontent.com/MaxineLiu/picBed/main/img/202111100243390.webp)

![](https://raw.githubusercontent.com/MaxineLiu/picBed/main/img/202111100243391.webp)

### 3）构建de Bruijn图

每一个component构建一个de Bruijn图，k-1个字节大小表示节点，k个字节大小表示连接节点的边，原始数据集中支持的（k-1）mer的数目作为边的权重。

#### a. de Bruijn图简化

合并de Bruijn图中的线性路径中的连续节点生成较长序列的节点，剔除可能由于测序错误（只有极少reads支持）的分叉，使边均匀。(多倍体多态性似乎比测序错误更常见，保留）

#### b. 寻找最佳路径

动态打分算法，鉴定被reads或双端reads支持的路径，剔除reads支持比较少的路径，将最佳路径上的碱基输出到fasta文件中。

![](https://raw.githubusercontent.com/MaxineLiu/picBed/main/img/202111100243393.webp)

## 结 果

按照 DNAXX，components, gene和 isoform分组的线性序列（见以下可左右滑动方框），Gene以下 isoform通常可能为同一基因的不同可变剪切。

> TRINITY\_DN6743\_c1\_g1\_i1 len:403_path:\[5739,5784,5857,5863,353\] TTGGGAGCCTGCCCAGGTTTTTGCTGGTACCAGGCTAAGTAGCTGCTAACACTCTGACTGGCCCGGCAGGTGATGGTGACTTTTTCCTCCTGAGACAAGGAGAGGGAGGCTGGAGACTGTGTCATCACGATTTCTCCGGTGATATCTGGGAGCCAGAGTAACAGAAGGCAGAGAAGGCGAGCTGGGGCTTCCATGGCTCACTCTGTGTCCTAACTGAGGCAGATCTCCCCCAGAGCACTGACCCAGCACTGATATGGGCTCTGGAGAGAAGAGTTTGCTAGGAGGAACATGCAAAGCAGCTGGGGAGGGGCATCTGGGCTTTCAGTTGCAGAGACCATTCACCTCCTCTTCTCTGCACTTGAGCAACCCATCCCCAGGTGGTCATGTCAGAAGACGCCTGGAG

