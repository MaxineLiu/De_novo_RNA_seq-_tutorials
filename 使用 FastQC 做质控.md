20160410 测序分析——使用 FastQC 做质控
=====================================

[孟浩巍](//www.zhihu.com/people/meng_howard)

[​](https://www.zhihu.com/question/48510028)

北京大学 生物学博士

**前言**

Hello 大家好！我是孟浩巍，一名生物信息方向的搬砖工！

之前我们已经介绍过了二代测序Illumina平台的测序原理，包括一些常用的专业术语，比如reads，flowcell，tail等等，如果大家有什么疑问请看之前的文章：

-   [Illumina测序原理介绍](http://zhuanlan.zhihu.com/p/20702684)

随后我们又介绍了，测序数据的储存格式FASTQ格式。今天要讲的内容会涉及到一些里面的知识。如果你忘记啦，请看看之前的文章：

-   [FASTA与FASTQ格式介绍](http://zhuanlan.zhihu.com/p/20714540)

Ok！那么今天我们要解决的问题是在从测序公司拿到原始数据以后，我们应该怎么评价这次的测序质量。是不是要做相应的一些后续处理。我们今天要使用的就是一个强大的工具——FastQC

**FastQC的基本介绍**

FastQC是一款基于Java的软件，一般都是在linux环境下使用命令行运行，它可以快速多线程地对测序数据进行质量评估（Quality Control），其官网地址为：[Babraham Bioinformatics](https://link.zhihu.com/?target=http%3A//www.bioinformatics.bbsrc.ac.uk/projects/fastqc/)

FastQC的下载和安装，和一般的Java软件没有什么区别，我们在这里就不做介绍了，在成功安装好以后，我们就在命令行模式下，输入fastqc就可以调用这个程序，界面如下：

![](https://pic1.zhimg.com/6cdffd5747cfc38375f481a141f96acc_b.png)

![](https://pic1.zhimg.com/80/6cdffd5747cfc38375f481a141f96acc_720w.png)

这时候我们可以选择 --help选项看一下帮助文档：

![](https://pic3.zhimg.com/315fc9215fc3242c2b4e1f44c4de0682_b.png)

![](https://pic3.zhimg.com/80/315fc9215fc3242c2b4e1f44c4de0682_720w.png)

    # 基本格式
    
    # fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam] [-c contaminant file] seqfile1 .. seqfileN
    
    # 主要是包括前面的各种选项和最后面的可以加入N个文件
    # -o --outdir FastQC生成的报告文件的储存路径，生成的报告的文件名是根据输入来定的
    # --extract 生成的报告默认会打包成1个压缩文件，使用这个参数是让程序不打包
    # -t --threads 选择程序运行的线程数，每个线程会占用250MB内存，越多越快咯
    # -c --contaminants 污染物选项，输入的是一个文件，格式是Name [Tab] Sequence，里面是可能的污染序列，如果有这个选项，FastQC会在计算时候评估污染的情况，并在统计的时候进行分析，一般用不到
    # -a --adapters 也是输入一个文件，文件的格式Name [Tab] Sequence，储存的是测序的adpater序列信息，如果不输入，目前版本的FastQC就按照通用引物来评估序列时候有adapter的残留
    # -q --quiet 安静运行模式，一般不选这个选项的时候，程序会实时报告运行的状况。

以我平时用的一个真实的例子：

    fastqc -o ./tmp.result/fastQC/ -t 6 ./tmp.data/fastq/H1EScell-dnase-2014-GSE56869_20151208_SRR1248176_1.fq 

使用的数据是2014年Dnase Hi-C的测序数据，数据下载地址：

    http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1370433

运行一段时间以后，就会出现报告：

    H1EScell-dnase-2014-GSE56869_20151208_SRR1248176_1.fq_fastqc.html
    H1EScell-dnase-2014-GSE56869_20151208_SRR1248176_1.fq_fastqc.zip

使用浏览器打开后缀是html的文件，就是图表化的fastqc报告：

![](https://pic1.zhimg.com/877967dad063658ee72cf567b4d69ab4_b.png)

![](https://pic1.zhimg.com/80/877967dad063658ee72cf567b4d69ab4_720w.png)

**FastQC的报告介绍**

-   总结信息

上图中Summary的部分就是整个报告的目录，整个报告分成若干个部分。合格会有个绿色的对勾，警告是个“!”，不合格是个红色的叉子。

-   基本信息

![](https://pic4.zhimg.com/b19f28a107a90f7c4f7b02510c7884ab_b.png)

![](https://pic4.zhimg.com/80/b19f28a107a90f7c4f7b02510c7884ab_720w.png)

    # Encoding指测序平台的版本和相应的编码版本号，这个在计算Phred反推error P的时候有用，如果不明白可以参考之前的文章。
    # Total Sequences记录了输入文本的reads的数量
    # Sequence length 是测序的长度
    # %GC 是我们需要重点关注的一个指标，这个值表示的是整体序列中的GC含量，这个数值一般是物种特意的，比如人类细胞就是42%左右。

-   序列测序质量统计

![](https://pic4.zhimg.com/38670ee6d5f373e326e3fe3e23ba4f9b_b.png)

![](https://pic4.zhimg.com/80/38670ee6d5f373e326e3fe3e23ba4f9b_720w.png)

    # 此图中的横轴是测序序列第1个碱基到第101个碱基
    # 纵轴是质量得分，Q = -10*log10（error P）即20表示1%的错误率，30表示0.1%
    # 图中每1个boxplot，都是该位置的所有序列的测序质量的一个统计，上面的bar是90%分位数，下面的bar是10%分位数，箱子的中间的横线是50%分位数，箱子的上边是75%分位数，下边是25%分位数
    # 图中蓝色的细线是各个位置的平均值的连线
    # 一般要求此图中，所有位置的10%分位数大于20,也就是我们常说的Q20过滤
    # 所以上面的这个测序结果，需要把后面的87bp以后的序列切除，从而保证后续分析的正确性
    # Warning 报警 如果任何碱基质量低于10,或者是任何中位数低于25
    # Failure 报错 如果任何碱基质量低于5,或者是任何中位数低于20

-   每个tail测序的情况

![](https://pic3.zhimg.com/ccf9bb69e20e2a3561581d08f57ef63e_b.png)

![](https://pic3.zhimg.com/80/ccf9bb69e20e2a3561581d08f57ef63e_720w.png)

    # 横轴和之前一样，代表101个碱基的每个不同位置
    # 纵轴是tail的Index编号
    # 这个图主要是为了防止，在测序过程中，某些tail受到不可控因素的影响而出现测序质量偏低
    # 蓝色代表测序质量很高，暖色代表测序质量不高，如果某些tail出现暖色，可以在后续分析中把该tail测序的结果全部都去除

-   每条序列的测序质量统计

![](https://pic4.zhimg.com/63086aca8162f21cb685d45c20e88b6f_b.png)

    # 假如我测的1条序列长度为101bp，那么这101个位置每个位置Q之的平均值就是这条reads的质量值
    # 该图横轴是0-40，表示Q值
    # 纵轴是每个值对应的reads数目
    # 我们的数据中，测序结果主要集中在高分中，证明测序质量良好！

-   GC 含量统计

![](https://pic1.zhimg.com/1b78609b8c4690eddf1eb1408c26f10c_b.png)

    # 横轴是1 - 101 bp；纵轴是百分比
    # 图中四条线代表A T C G在每个位置平均含量
    # 理论上来说，A和T应该相等，G和C应该相等，但是一般测序的时候，刚开始测序仪状态不稳定，很可能出现上图的情况。像这种情况，即使测序的得分很高，也需要cut开始部分的序列信息，一般像我碰到这种情况，会cut前面5bp

-   序列平均GC含量分布图

![](https://pic3.zhimg.com/5fead4df4b5c980b7519ddd536a7b196_b.png)



    # 横轴是0 - 100%； 纵轴是每条序列GC含量对应的数量
    # 蓝色的线是程序根据经验分布给出的理论值，红色是真实值，两个应该比较接近才比较好
    # 当红色的线出现双峰，基本肯定是混入了其他物种的DNA序列
    # 这张图中的信息良好

-   序列测序长度统计

![](https://pic1.zhimg.com/97d24fbdf9e42fc057072f7582e1532c_b.png)

![](https://pic1.zhimg.com/80/97d24fbdf9e42fc057072f7582e1532c_720w.png)

    # 每次测序仪测出来的长度在理论上应该是完全相等的，但是总会有一些偏差
    # 比如此图中，101bp是主要的，但是还是有少量的100和102bp的长度，不过数量比较少，不影响后续分析
    # 当测序的长度不同时，如果很严重，则表明测序仪在此次测序过程中产生的数据不可信 

-   序列Adapter

![](https://pic2.zhimg.com/7563b88d741b0f97e8546325f91e9969_b.png)

![](https://pic2.zhimg.com/80/7563b88d741b0f97e8546325f91e9969_720w.png)

    # 此图衡量的是序列中两端adapter的情况
    # 如果在当时fastqc分析的时候-a选项没有内容，则默认使用图例中的四种通用adapter序列进行统计
    # 本例中adapter都已经去除，如果有adapter序列没有去除干净的情况，在后续分析的时候需要先使用cutadapt软件进行去接头，这个软件以后我会介绍

-   重复短序列

![](https://pic4.zhimg.com/d36385723bfbe857c1d0516c622fc05f_b.png)

![](https://pic4.zhimg.com/80/d36385723bfbe857c1d0516c622fc05f_720w.png)

    # 这个图统计的是，在序列中某些特征的短序列重复出现的次数
    # 我们可以看到1-8bp的时候图例中的几种短序列都出现了非常多的次数，一般来说，出现这种情况，要么是adapter没有去除干净，而又没有使用-a参数；要么就是序列本身可能重复度比较高，如建库PCR的时候出现了bias
    # 对于这种情况，我的办法是可以cut掉前面的一些长度，可以试着cut 5~8bp

**本期总结与下期预告**

在本片文章中，我们已经介绍了测序质量评估最常用的FastQC软件，并详细解读了报告内容。希望能够对大家有所帮助。可是我们还是留有一些问题，比如我一直说cut序列，去adapter，却没有给出工具的用法。所以下期我们的主题主要是，围绕着FastQC中出现的问题使用不同的fastx\_trimmer，cutadapt等工具进行修正。

希望各位多多点赞支持！

孟浩巍

——————————————————————————

另外欢迎各位参加我们的知乎Live：

1. 知乎Live：如何快速入门生物信息学 （涉及内容：测序原理，生物信息学发展历史，软件的安装与调试，入门路线图，介绍了RNA-Seq的分析流程并给出实践代码）；

[如何快速入门生物信息学](https://www.zhihu.com/lives/851804180241846272)

2. 知乎Live: 生信进阶第1课-重复Nature文章 (涉及内容：肺癌相关研究现状，RNA-Seq单细胞测序，RNA-Seq的建库方法，RNA-Seq的分析流程细节，相关生信图的绘制）；

[生信进阶第1课-重复Nature文章](https://www.zhihu.com/lives/865204175334670336)

3. 知乎Live：生信进阶第2课-基因组序列

(涉及内容：介绍基因组的序列结构，hg19与hg38的区别，ENCODE计划，常用的表观组学实验原理ChIP-Seq，Hi-C等，ChIP-Seq的标准处理流程，绘图原理)

[生信进阶第2课-基因组序列](https://www.zhihu.com/lives/880463398117208064)

4. 知乎Live：不用编程怎么做生物信息学

(涉及内容：介绍生物信息学入门的几个层次，从命令行到图形界面再到命令行，绘制生物进化树，图形界面分析平台，使用图形界面处理RNA-Seq数据，使用图形界面分析ChIP-Seq数据，UCSC genome browser，WashU genome browser)

[不用编程怎么做生物信息学](https://www.zhihu.com/lives/893132277473746944)

编辑于 2017-09-21
