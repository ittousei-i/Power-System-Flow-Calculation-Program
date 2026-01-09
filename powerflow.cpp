// powerflow.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

//程序功能：利用牛顿-拉夫逊迭代法，使用直角坐标对n节点n支路电力系统进行潮流计算
// ------------------------------------------------------------
// 功能特点：
// 1. 不采用稀疏矩阵技术，直接使用全矩阵
// 2. 支持三种节点类型：
//    - 平衡节点（Slack）
//    - PV节点
//    - PQ节点
// 3. 输入、输出均采用文本文件，便于调试和打印
// 4. 输出内容包括：
//    - 原始输入参数
//    - 节点导纳矩阵 Ybus
//    - 收敛后的各节点电压
//    - 各节点注入功率
//    - 各支路两端功率及功率损耗
// ------------------------------------------------------------

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
using namespace std;

//定义常量
const double PI = 3.14159265358979323846;

//封装一个复数类，便于进行运算
class Complex
{
private:
    double real, imag;
public:
    //构造函数
    //无参构造 复数设为 0 + j0
    Complex() :real(0), imag(0) {}
    //有参构造 复数设为 real + j*imag
    Complex(double r, double i) :real(r), imag(i) {}
    //使用重载运算符定义复数运算
    //加法
    Complex operator+(const Complex& other) const
    {
        return Complex(real + other.real, imag + other.imag);
    }
    //减法
    Complex operator-(const Complex& other) const
    {
        return Complex(real - other.real, imag - other.imag);
    }

    //乘法
    Complex operator*(const Complex& other) const
    {
        double new_real = real * other.real - imag * other.imag;
        double new_imag = real * other.imag + imag * other.real;
        return Complex(new_real, new_imag);
    }

    //除法
    Complex operator/(const Complex& other) const
    {
        //计算分母
        double denominator = other.real * other.real + other.imag * other.imag;
        //避免除以零
        if (denominator < 1e-9) //浮点型不能直接判0，用极小值判断
        {
            throw invalid_argument("除数不能为0（复数0+0i）");
        }
        double new_real = (real * other.real + imag * other.imag) / denominator;
        double new_imag = (imag * other.real - real * other.imag) / denominator;
        return Complex(new_real, new_imag);
    }

    //取模运算：计算复数的模 |z| = √(real² + imag²)
    double modulus() const
    {
        return sqrt(real * real + imag * imag);
    }

    //共轭运算：返回共轭复数（实部不变，虚部取反）
    Complex conjugate() const
    {
        return Complex(real, -imag);
    }

    //重载<<运算符，方便输出复数（友元函数，能访问私有成员）
    friend ostream& operator<<(ostream& os, const Complex& c)
    {
        //格式化输出，处理虚部正负的情况
        if (c.imag >= 0)
        {
            os << c.real << " + " << "j" << c.imag;
        }
        else
        {
            os << c.real << " - " << "j" << -c.imag;
        }
        return os;
    }

    //获取实部和虚部
    double getReal() const { return real; }
    double getImag() const { return imag; }
};

//节点结构体
struct Bus
{
    int id; //节点编号
    int type; //节点类型 0=Slack,1=PV,2=PQ
    double P, Q; //注入节点的有功和无功功率
    double e, f, delta, U; //电压的实部、虚部、相角和模值
    double dP, dQ, du, de, df; //不平衡量
};

//支路结构体
struct Branch
{
    int id; //支路编号
    int from, to; //构成支路的两节点的编号
    double R, X; //该支路的电阻和电抗 阻抗Z=R+jX
};

/*以下为函数部分*/

//读取文件用工具函数 功能为遇注释和空行则跳过
bool readDataLine(ifstream& fin, string& line)
{
    while (getline(fin, line))
    {
        if (line.empty()) continue;
        if (line[0] == '#') continue;
        return true;
    }
    return false;
}

//读取基础参数（节点数、支路数、计算精度）
bool readBasicParams(ifstream& fin, int& busNum, int& branchNum, double& epsilon)
{
    string line;
    if (!readDataLine(fin, line))
    {
        return false;
    }
    stringstream ss(line);
    ss >> busNum >> branchNum >> epsilon;
    return true;
}

//读取节点数据
void readBusData(ifstream& fin, int busNum, vector<Bus>& bus)
{
    bus.resize(busNum);
    string line;
    for (int i = 0; i < busNum; i++)
    {
        readDataLine(fin, line);
        stringstream sb(line);
        sb >> bus[i].id >> bus[i].type >> bus[i].P >> bus[i].Q
            >> bus[i].U >> bus[i].delta;

        if (bus[i].type == 2)
        {
            bus[i].e = bus[i].U;
            bus[i].f = 0.0;
        }
        else if (bus[i].type == 1)
        {
            bus[i].e = bus[i].U;
            bus[i].f = 0.0;
        }
        else
        {
            bus[i].e = 1.0;
            bus[i].f = 0.0;
        }
    }
}

//读取支路数据
void readBranchData(ifstream& fin, int branchNum, vector<Branch>& branch)
{
    branch.resize(branchNum);
    string line;
    for (int i = 0; i < branchNum; i++)
    {
        readDataLine(fin, line);
        stringstream sl(line);
        sl >> branch[i].id >> branch[i].from >> branch[i].to
            >> branch[i].R >> branch[i].X;
    }
}

//构建节点导纳矩阵Ybus及其实部G、虚部B矩阵
void buildYbus(int busNum, int branchNum, const vector<Branch>& branch,
    vector<vector<Complex>>& Ybus, vector<vector<double>>& G, vector<vector<double>>& B)
{
    //初始化Ybus
    Ybus.assign(busNum, vector<Complex>(busNum, Complex(0, 0)));

    //遍历每条支路，构建节点导纳矩阵
    for (int i = 0; i < branchNum; i++)
    {
        int from = branch[i].from - 1; //转换为0-based索引
        int to = branch[i].to - 1; //转换为0-based索引

        //计算支路导纳 y = 1/(R + jX)
        Complex Z(branch[i].R, branch[i].X);  //支路阻抗
        Complex y = Complex(1, 0) / Z;        //支路导纳

        //更新节点导纳矩阵
        //自导纳：对角元素加上支路导纳
        Ybus[from][from] = Ybus[from][from] + y;
        Ybus[to][to] = Ybus[to][to] + y;

        //互导纳：非对角元素减去支路导纳
        Ybus[from][to] = Ybus[from][to] - y;
        Ybus[to][from] = Ybus[to][from] - y;
    }

    //将Ybus分解为实部G和虚部B矩阵
    G.assign(busNum, vector<double>(busNum, 0));
    B.assign(busNum, vector<double>(busNum, 0));
    for (int i = 0; i < busNum; i++)
    {
        for (int j = 0; j < busNum; j++)
        {
            G[i][j] = Ybus[i][j].getReal();
            B[i][j] = Ybus[i][j].getImag();
        }
    }
}

//输出原始输入参数和节点导纳矩阵到文件
void outputOriginalParams(ofstream& fout, int busNum, int branchNum, double epsilon,
    const vector<Bus>& bus, const vector<Branch>& branch,
    const vector<vector<Complex>>& Ybus)
{
    //输出基础参数
    fout << "============= 原始输入参数 =============" << endl;
    fout << "节点数: " << busNum << " 支路数: " << branchNum << " 计算精度: " << epsilon << endl;
    fout << endl;

    //输出节点数据
    fout << "============= 节点数据 =============" << endl;
    fout << "#编号/类型/有功/无功/电压幅值/电压相角" << endl;
    for (int i = 0; i < busNum; i++)
    {
        fout << bus[i].id << " " << bus[i].type << " "
            << bus[i].P << " " << bus[i].Q << " "
            << bus[i].U << " " << bus[i].delta << endl;
    }
    fout << endl;

    //输出支路数据
    fout << "============= 支路数据 =============" << endl;
    fout << "#编号/节点1/节点2/电阻/电抗" << endl;
    for (int i = 0; i < branchNum; i++)
    {
        fout << branch[i].id << " " << branch[i].from << " " << branch[i].to << " "
            << branch[i].R << " " << branch[i].X << endl;
    }
    fout << endl;

    //输出节点导纳矩阵
    fout << "============= 节点导纳矩阵 Ybus =============" << endl;
    for (int i = 0; i < busNum; i++)
    {
        fout << "第" << i + 1 << "行: ";
        for (int j = 0; j < busNum; j++)
        {
            fout << Ybus[i][j];
            if (j < busNum - 1)
            {
                fout << "  ";
            }
        }
        fout << endl;
    }
    fout << endl;
}

//牛顿-拉夫逊迭代求解潮流方程
bool newtonRaphsonIteration(int busNum, double epsilon, int maxIter,
    vector<Bus>& bus, const vector<vector<double>>& G,
    const vector<vector<double>>& B, int& iter)
{
    //初始化：如果电压幅值为0（代表未知），赋予平直启动初值 1.0
    for (int i = 0; i < busNum; i++)
    {
        if (bus[i].type == 2 && bus[i].e == 0.0) //PQ节点平直启动
        {
            bus[i].e = 1.0;
            bus[i].f = 0.0;
        }
    }

    iter = 0;
    bool converged = false;
    int n = busNum - 1; //除去平衡节点的计算节点数
    int dim = 2 * n;    //不平衡量维数

    while (iter < maxIter && !converged)
    {
        iter++;
        vector<double> dM(dim, 0.0); //不平衡量向量 [dP1, dQ1, dP2, dU2^2...]
        vector<vector<double>> J(dim, vector<double>(dim, 0.0)); //雅可比矩阵
        double maxDelta = 0.0;

        //1.计算不平衡量和填充雅可比矩阵
        for (int i = 0; i < n; i++)
        {
            double ai = 0, bi = 0;
            for (int j = 0; j < busNum; j++)
            {
                ai += (G[i][j] * bus[j].e - B[i][j] * bus[j].f);
                bi += (G[i][j] * bus[j].f + B[i][j] * bus[j].e);
            }

            //计算 Pi, Qi
            double Pi = bus[i].e * ai + bus[i].f * bi;
            double Qi = bus[i].f * ai - bus[i].e * bi;

            //计算 dP
            bus[i].dP = bus[i].P - Pi;
            dM[2 * i] = bus[i].dP;
            if (abs(bus[i].dP) > maxDelta) maxDelta = abs(bus[i].dP);

            //PQ 节点：计算 dQ
            if (bus[i].type == 2)
            {
                bus[i].dQ = bus[i].Q - Qi;
                dM[2 * i + 1] = bus[i].dQ;
                if (abs(bus[i].dQ) > maxDelta) maxDelta = abs(bus[i].dQ);
            }
            //PV 节点：计算 dU^2
            else if (bus[i].type == 1)
            {
                double Ui2 = bus[i].e * bus[i].e + bus[i].f * bus[i].f;
                bus[i].du = bus[i].U * bus[i].U - Ui2;
                dM[2 * i + 1] = bus[i].du;
                if (abs(bus[i].du) > maxDelta) maxDelta = abs(bus[i].du);
            }

            //2.填充雅可比矩阵元素
            for (int j = 0; j < n; j++)
            {
                //对角块
                if (i == j)
                {
                    //dP/df, dP/de
                    J[2 * i][2 * j] = bi + G[i][i] * bus[i].f - B[i][i] * bus[i].e; //Hii = dPi/dfi
                    J[2 * i][2 * j + 1] = ai + G[i][i] * bus[i].e + B[i][i] * bus[i].f; //Nii = dPi/dei
                    //PQ: dQ/df, dQ/de
                    if (bus[i].type == 2)
                    {
                        J[2 * i + 1][2 * j] = ai - G[i][i] * bus[i].e - B[i][i] * bus[i].f; //Kii = dQi/dfi
                        J[2 * i + 1][2 * j + 1] = -bi + G[i][i] * bus[i].f - B[i][i] * bus[i].e; //Lii = dQi/dei
                    }
                    //PV: dU^2/df, dU^2/de
                    else
                    {
                        J[2 * i + 1][2 * j] = 2 * bus[i].f;
                        J[2 * i + 1][2 * j + 1] = 2 * bus[i].e;
                    }
                }
                else //非对角块
                {
                    J[2 * i][2 * j] = G[i][j] * bus[i].f - B[i][j] * bus[i].e; //Hij
                    J[2 * i][2 * j + 1] = G[i][j] * bus[i].e + B[i][j] * bus[i].f; //Nij

                    //PQ
                    if (bus[i].type == 2)
                    {
                        J[2 * i + 1][2 * j] = -J[2 * i][2 * j + 1]; //Kij = -Nij
                        J[2 * i + 1][2 * j + 1] = J[2 * i][2 * j];  //Lij = Hij
                    }
                    //PV
                    else
                    {
                        J[2 * i + 1][2 * j] = 0;
                        J[2 * i + 1][2 * j + 1] = 0;
                    }
                }
            }
        }

        //检查收敛
        if (maxDelta < epsilon)
        {
            converged = true;
            break;
        }

        //3.解线性方程组 J * dX = dM (高斯消元法)
        for (int i = 0; i < dim; i++)
        {
            int pivot = i;
            for (int j = i + 1; j < dim; j++)
            {
                if (abs(J[j][i]) > abs(J[pivot][i]))
                {
                    pivot = j;
                }
            }
            swap(J[i], J[pivot]);
            swap(dM[i], dM[pivot]);

            for (int j = i + 1; j < dim; j++)
            {
                double factor = J[j][i] / J[i][i];
                dM[j] -= factor * dM[i];
                for (int k = i; k < dim; k++)
                {
                    J[j][k] -= factor * J[i][k];
                }
            }
        }
        vector<double> dX(dim);
        for (int i = dim - 1; i >= 0; i--)
        {
            double sum = 0;
            for (int j = i + 1; j < dim; j++)
            {
                sum += J[i][j] * dX[j];
            }
            dX[i] = (dM[i] - sum) / J[i][i];
        }

        //4.更新电压 [df1, de1, df2, de2...]
        for (int i = 0; i < n; i++)
        {
            bus[i].f += dX[2 * i];
            bus[i].e += dX[2 * i + 1];
        }
    }

    return converged;
}

//计算平衡节点功率、节点电压幅值/相角、支路功率及损耗
void calculatePower(int busNum, int branchNum, vector<Bus>& bus,
    const vector<Branch>& branch, const vector<vector<double>>& G,
    const vector<vector<double>>& B)
{
    //更新所有节点的幅值和相角，计算节点注入功率
    for (int i = 0; i < busNum; i++)
    {
        bus[i].U = sqrt(bus[i].e * bus[i].e + bus[i].f * bus[i].f);
        bus[i].delta = atan2(bus[i].f, bus[i].e) * 180.0 / PI;

        //计算节点的注入功率（包括平衡节点和PV节点的无功）
        double ai = 0, bi = 0;
        for (int j = 0; j < busNum; j++)
        {
            ai += (G[i][j] * bus[j].e - B[i][j] * bus[j].f);
            bi += (G[i][j] * bus[j].f + B[i][j] * bus[j].e);
        }
        bus[i].P = bus[i].e * ai + bus[i].f * bi;
        bus[i].Q = bus[i].f * ai - bus[i].e * bi;
    }
}

//输出潮流计算结果到文件
void outputResults(ofstream& fout, int busNum, int branchNum, int iter,
    const vector<Bus>& bus, const vector<Branch>& branch)
{
    //输出迭代结果
    fout << "============= 计算结果 (迭代次数: " << iter << ") =============" << endl;
    fout << "节点编号  电压幅值U    电压相角(deg)  注入有功P    注入无功Q" << endl;
    for (int i = 0; i < busNum; i++)
    {
        fout << bus[i].id << "\t" << bus[i].U << "\t" << bus[i].delta << "\t"
            << bus[i].P << "\t" << bus[i].Q << endl;
    }
    fout << endl;

    //输出支路功率
    fout << "============= 支路功率流动 =============" << endl;
    fout << "支路  从节点->到节点    始端功率(Sij)          末端功率(Sji)          功率损耗" << endl;
    for (int i = 0; i < branchNum; i++)
    {
        int u = branch[i].from - 1;
        int v = branch[i].to - 1;
        Complex Vij(bus[u].e, bus[u].f);
        Complex Vji(bus[v].e, bus[v].f);
        Complex yij = Complex(1, 0) / Complex(branch[i].R, branch[i].X);

        Complex Sij = Vij * ((Vij - Vji) * yij).conjugate();
        Complex Sji = Vji * ((Vji - Vij) * yij).conjugate();
        Complex Loss = Sij + Sji;

        fout << branch[i].id << "     " << branch[i].from << " -> " << branch[i].to << "    "
            << Sij << "    " << Sji << "    " << Loss << endl;
    }
}

//主函数
int main()
{
    //1.文件输入输出初始化
    ifstream fin("input.txt");
    ofstream fout("output.txt");
    if (!fin)
    {
        cout << "无法打开 input.txt" << endl;
        return 1;
    }

    //2.读取基础参数
    int busNum, branchNum;
    double epsilon;
    if (!readBasicParams(fin, busNum, branchNum, epsilon))
    {
        cout << "读取基础参数失败" << endl;
        fin.close();
        fout.close();
        return 1;
    }

    //3.读取节点和支路数据
    vector<Bus> bus;
    vector<Branch> branch;
    readBusData(fin, busNum, bus);
    readBranchData(fin, branchNum, branch);

    //4.构建节点导纳矩阵
    vector<vector<Complex>> Ybus;
    vector<vector<double>> G, B;
    buildYbus(busNum, branchNum, branch, Ybus, G, B);

    //5.输出原始输入参数
    outputOriginalParams(fout, busNum, branchNum, epsilon, bus, branch, Ybus);

    //6.牛顿-拉夫逊迭代求解
    int iter = 0;
    const int maxIter = 20;
    newtonRaphsonIteration(busNum, epsilon, maxIter, bus, G, B, iter);

    //7.计算节点功率和支路功率
    calculatePower(busNum, branchNum, bus, branch, G, B);

    //8.输出计算结果
    outputResults(fout, busNum, branchNum, iter, bus, branch);

    //资源释放
    cout << "计算完成，请打开 output.txt 查看潮流计算结果" << endl;
    fin.close();
    fout.close();

    return 0;
}
