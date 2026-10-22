Figure 2：完整R计算与绘图流程

操作：
1. 全部解压本包，保留 inputs、reference、code 文件夹。
2. 在RStudio打开最外层 RUN_FIGURE2.R，点击 Source。
3. 等待控制台出现 PASS_EXPRESSION_TO_FIGURE2。
4. results/本次运行时间与进程编号/figure/ 内为Figure2_R.png、PDF及TIFF。

需要：已安装GSVA 2.4.9。你刚才Figure3完整流程已成功使用这个版本，无需再安装。
其余步骤使用基础R，不自动安装或更新任何包。

计算链：
18833基因×40样本的归档表达矩阵 + 样本分组/性别 + 85个冻结基因集
→ 重新计算85×40个GSVA分数
→ 全部85个基因集各自拟合 score ~ group + sex，原始GSVA单位
→ HC3标准误、残差自由度37的t检验和普通95%置信区间
→ 在85项内进行BH校正
→ 按FDR<0.05选出8项，绘制森林图
→ 同8项在全40人中按行进行样本标准差标准化，绘制热图
→ 逐项核对归档统计结果和全部8×40个显示值。

森林图：调整性别的效应及普通95%置信区间；不是标准化系数。
热图：描述性、未调整的个体分数；对照10人在前，SONFH30人在后。
组内保留归档样本顺序，不聚类。蓝色低于该行全40人的平均值，橙色高于平均值。
该热图不是独立验证；颜色也不是不同基因集之间绝对数量的比较。
FDR在全部85项中计算，不是仅在展示的8项中计算。

主要输出：
FINAL_STATUS.txt：本次完整流程是否通过。
numerical_checks.csv：每项核对的误差与阈值。
sessionInfo.txt：本机软件版本。
input_and_code_MD5.csv：实际输入和代码校验值。
GSVA_scores_RECALCULATED.csv：全部新计算的原始分数。
Layer1_HC3_all85_RECALCULATED.csv：全部85项新拟合的统计结果。
figure/Figure2_estimates_RECALCULATED.csv：图中8项的新统计结果。
figure/Figure2_display_z_RECALCULATED.csv：热图显示值。
figure/Figure2_sample_order.csv：列顺序、分组和性别。

范围：从归档基因表达矩阵到Figure2；不包括CEL-to-RMA或探针注释汇总。
reference中旧统计结果仅用于核对，新模型不会读取旧效应/P值作为拟合输入。
R版图形字体与布局可能略有差异，数值及解释应与原图一致。
本准备环境无R，新的Figure2入口尚未实际执行。请以本机最终PASS和报告为准。
如果报错，请复制完整错误文字，不必反复重启R。
