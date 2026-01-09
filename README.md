# Power-System-Flow-Calculation-Program
利用牛顿-拉夫逊迭代法计算复杂电力系统潮流的C++程序
# 使用方法
在程序根目录下新建 input.txt 文件，按以下格式完成输入文件：
#节点数/支路数/计算精度
X X X.XXXX
#节点数据
#编号/类型/有功/无功/电压幅值/电压相角
X X X X X X

#支路数据
#编号/节点1/节点2/电阻/电抗
X X X X X
运行程序后会在程序根目录生成 output.txt ，其中即为潮流计算的结果
