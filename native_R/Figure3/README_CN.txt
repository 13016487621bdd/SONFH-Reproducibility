Figure 3：表达矩阵到统计结果和图片的R运行包

这是完整计算版，与之前只读取统计表画图的Figure3_R包不同。

操作
1. 完整解压本包到桌面，保留 inputs、reference、code 文件夹。
2. 用RStudio打开最外层 RUN_FIGURE3.R。
3. 点击右上角 Source，等待运行完成。不要只运行选中的一行。
4. 成功时控制台最后显示 PASS_EXPRESSION_TO_FIGURE3。
5. 所有结果保存在 results/本次时间与进程编号/，图在其中
   figure_build/results_R/Figure3_R.png（同时输出PDF，并尝试输出TIFF）。

需要已安装的历史版本：GSVA 2.4.9、clubSandwich 0.7.0。
你之前的运行记录中已有这两个版本，通常无需再次安装。
本程序不自动安装、升级或降级软件；版本不符会明确停止。
解压到磁盘后再运行，不要直接从压缩包中双击R脚本。

计算顺序
1. 读取18833基因×40样本表达矩阵、GSM对应的分组和性别、85个冻结基因集。
2. 使用GSVA Gaussian/tau=1/maxDiff=TRUE/absRanking=FALSE/minSize=1/maxSize=Inf
   重新计算85×40个分数，并与归档分数逐项比较。
3. 各分数在40个样本中按样本标准差标准化。依冻结组内顺序选首个基因集为参照，
   调用原始Layer2函数重新拟合全部21组；参与者GSM聚类CR2，两个基因集时
   Satterthwaite检验，三个及以上时HTZ联合检验；21项内单独BH校正。
4. 依11对冻结A/B方向计算z(A)-z(B)，不再标准化这个差值。
   用组别+性别线性模型重新估计疾病系数，HC3标准误、残差自由度37的t检验
   及普通95%置信区间，11项内单独BH校正。
5. 与归档统计结果逐项核对后，将新计算的统计表传给R绘图函数，输出Figure3。
两个分析层为平行的预定义分析，不做逐层显著性筛选。

重要输出
FINAL_STATUS.txt：整次运行是否通过。只有完整计算、核对、绘图结束才写PASS。
numerical_checks.csv：每项核对的最大绝对误差及阈值（默认1e-8）。
sessionInfo.txt：本机R与软件版本。
input_and_code_MD5.csv：实际运行的输入及代码校验值。
GSVA_scores_RECALCULATED.csv：新计算的全部GSVA分数。
GSVA_z_scores_RECALCULATED.csv：标准化分数。
Layer2_members_and_references.csv：功能组成员及参照基因集。
Layer2_CR2_RECALCULATED.csv：21组新拟合结果。
Layer3_participant_contrasts_RECALCULATED.csv：每位参与者的11个差值。
Layer3_HC3_RECALCULATED.csv：11个对比的新拟合结果。
figure_build/inputs/：直接交给绘图函数的新统计表。
figure_build/results_R/：R生成的图、显示数据及绘图核对报告。

范围和验证状态
本包从归档的基因层面表达矩阵开始；不包含CEL-to-RMA、探针注释汇总或组织分析。
基因集定义和组内成员来自冻结记录；master文件中的既有beta/P/FDR不参与新模型拟合。
reference中的旧结果仅用于比较；图片由重新拟合后的结果生成。
原始Layer2函数按原样保留在code/layer2_original_function.R。
绘图代码来源于你已成功运行的基础R版，增加了接收新输入目录的参数。
本准备环境没有R，因此这个新整合入口尚未实际运行；不能把先前绘图PASS当作
本整合流程已通过。请运行一次，以实际FINAL_STATUS及核对报告为准。
如有报错，复制控制台完整文字；不需要不断重启R，也不要自行改动输入数据或阈值。
