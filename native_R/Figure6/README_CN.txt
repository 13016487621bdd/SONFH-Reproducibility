Figure 6：原始计数到组织图的完整R复现包

全部解压，在RStudio打开RUN_FIGURE6.R，点击Source。
成功标志：PASS_COUNTS_TO_FIGURE6。
results/本次时间与进程编号/保存全部结果，figure/保存PNG、PDF、TIFF和英文图注。
软件：DESeq2 1.50.2、GSVA 2.4.9，以及readxl和SummarizedExperiment。
你之前组织审计已运行过这些包；本程序不自动安装或升级。

计算流程：
20345基因×9样本原始计数→ratio size factors→固定VST配方
(design ~1; blind=TRUE; nsub=1000; fitType=parametric)
→与归档全VST矩阵核对→冻结注释过滤为18572个一对一映射
→追加ENSG00000112096为SOD2，得到18573×9矩阵并逐值核对
→GSVA重新评分及九个样本内的组成分数z标准化
→原始及去除DUOX1/DUOX2两种情景的生成减分解差值
→全部126种5/9标签分配的双侧精确检验→重新画图。
另外输出两个情景的固定分数留一法均值差，原始情景与归档LOO核对。

两个面板均绘制4例骨折对照和5例SONFH的全部个体点；短横线是组均值，不是CI。
预期原始差值约1.2937，精确P=5/126=0.0397；
去除DUOX1/DUOX2后差值约1.3341，精确P=7/126=0.0556。
这里效应没有减弱，不能仅根据P越过0.05就写“效应衰减”。
两种情景各自计算组成分数的均值/标准差，差值不再次标准化。
组织对照是骨折对照，不是GC暴露匹配对照；分析属探索性方向评估。

重要输出：numerical_checks.csv、FINAL_STATUS.txt、sessionInfo.txt、输入和代码MD5、
upstream_RECONSTRUCTED.rds、Tissue_effects_RECALCULATED.csv、
Tissue_participant_scores_RECALCULATED.csv、All126_permutations_by_version.csv、
Tissue_fixed_score_LOO.csv，以及figure/内的图片和图注。

范围与状态：
从提供的原始基因计数开始，不包括FASTQ比对/定量，不实时查询新注释数据库。
冻结注释及基因集是输入，归档矩阵/分数仅用于核对。
VST配方依据你此前全矩阵零误差的恢复审计固定，不再搜索候选配方，不按P值选配方。
这是等价重建流程，不声称找回唯一的历史源码或唯一历史设置。
当前准备环境没有R，新入口和绘图尚未实际运行/目视检查。
请提供本机运行末尾输出及Figure6_R.png，以确认整合流程和排版。
