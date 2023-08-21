#ifndef _KMC_INIT_H_
#define _KMC_INIT_H_
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "kmc_input.h"
#include "blitz/array.h"


//晶格基类
class Base {
public:
    int n1nbr;                          // 第一近邻数
    int n2nbr;                          // 第二近邻数
    int Length;                         // 长度
    int Bnum;                           // 晶胞中B原子数
    int Vnum;                           // 晶胞中V原子数
    int Anum;                           // 晶胞中A原子数
    int Cnum;                           // 晶胞中C原子数
    int Dnum;                           // 晶胞中D原子数
    int vacx;                           // 空位坐标x
    int vacy;                           // 空位坐标y
    int vacz;                           // 空位坐标z
    int vacnm;                          // 空位编号
    blitz::Array<int, 4> Site;          // 网格点坐标
    blitz::Array<double, 2> BList;      // B能量列表
    blitz::Array<double, 2> AList;      // A能量列表
    std::vector<int>modla;              // 模量向量la
    std::vector<int>modlb;              // 模量向量lb
    std::vector<int>modlc;              // 模量向量lc
    int(*modlx);                        // 模量向量lx
    int(*modly);                        // 模量向量ly
    int(*modlz);                        // 模量向量lz
    typedef struct {                    // 定义空位结构体
        int x;                          // 空位坐标x
        int y;                          // 空位坐标y
        int z;                          // 空位坐标z
        int nm;                         // 空位编号
    } VAC;

    blitz::Array<std::vector<VAC>, 3> Vaclists;    // 空位列表
    //SimulationParameters p;//指向 SimulationParameters 类型的指针

    virtual ~Base() = default;           // 析构函数
    void Energy(SimulationParameters& parameter); // 计算能量
    void InitStatesArray(SimulationParameters& parameter); // 初始化晶胞
    void Periodic(SimulationParameters& parameter); // 处理周期性边界条件
    void RecycleFNeighbors(SimulationParameters& parameter); // 更新FCC近邻
    void RecycleBNeighbors(SimulationParameters& parameter); // 更新BCC近邻

    virtual void CountNeighbors(int nm, int i, int  j, int  k, std::vector<int> count1, std::vector<int> count2) = 0; // 计算近邻数

private:
    std::vector<int> CalculateN1Neighbors(SimulationParameters& parameter, int nm, int nmn, int i, int  j, int  k);
    std::vector<int> CalculateN1FNeighbors(SimulationParameters& parameter, int nm, int nmn, int i, int  j, int  k);
    std::vector<int> CalculateN2Neighbors(SimulationParameters& parameter, int nm, int nmn, int i, int  j, int  k);
    void CalculateAtomNumbers(SimulationParameters& parameter); // 计算原子数量
    void GenerateLatticePoints(SimulationParameters& parameter); // 生成晶格点
    void AddAtom(SimulationParameters& parameter, int type, int ii, int jj, int zz, int nm); // 添加原子
    void TofileStatesArray(SimulationParameters& parameter); // 将晶胞状态写入文件
    void AbInitStatesArrays(SimulationParameters& parameter); // 初始化晶胞状态
    std::vector<std::string> Split(const std::string& s, const char& delim);

};

// 声明一个单空位BCC晶格类SingleBcc，继承自基类Base
class SingleBcc : public Base {
public:
    blitz::Array<int, 4> NN1, NN2; // 4维整数数组，存储每个位置的一级和二级邻居数
    SimulationParameters* p;
    // SingleBcc类的构造函数，输入参数为SimulationParameters的引用
    SingleBcc(SimulationParameters& parameter){
        p = &parameter;
        Periodic(parameter); // 执行周期性边界条件的初始化
        n1nbr = 8, n2nbr = 6, Length = 2;// 一级邻居数，二级邻居数，晶格长度的赋值
        // 根据晶格长度和模拟参数的xyz尺寸初始化Site数组，并赋初值为0
        Site.resize(Length, parameter.nx, parameter.ny, parameter.nz); Site = 0;
        NN1.resize(Length, parameter.nx, parameter.ny, parameter.nz); NN1 = 0;
        NN2.resize(Length, parameter.nx, parameter.ny, parameter.nz); NN2 = 0;
        // 根据给定范围初始化BList数组，并赋初值为0.0
        BList.resize(blitz::Range(-8, 8), blitz::Range(-6, 6)); BList = 0.0;
        AList.resize(blitz::Range(-8, 8), blitz::Range(-6, 6));AList = 0.0;

        Energy(parameter); // 计算能量
        InitStatesArray(parameter); // 初始化晶格状态数组
        RecycleBNeighbors(parameter); // 计算每个位置的一级和二级邻居数
        
    }
    // 重写CountNeighbors函数，计算每个位置的一级和二级邻居数
    void CountNeighbors(int nm, int i, int  j, int  k, std::vector<int> count1, std::vector<int> count2) override;
 
};
// 声明一个单空位FCC晶格类SingleFcc，继承自基类Base  同BCC
class SingleFcc : public Base {
public:
    blitz::Array<int, 4>NN1, NN2;
    SimulationParameters* p;
    SingleFcc(SimulationParameters& parameter) {
        p = &parameter;
        Periodic(parameter);
        n1nbr = 12, n2nbr = 6, Length = 4;
        Site.resize(Length, parameter.nx, parameter.ny, parameter.nz); Site = 0;
        NN1.resize(Length, parameter.nx, parameter.ny, parameter.nz); NN1 = 0;
        NN2.resize(Length, parameter.nx, parameter.ny, parameter.nz); NN2 = 0;
        BList.resize(blitz::Range(-12, 12), blitz::Range(-6, 6)); BList = 0.0;
        AList.resize(blitz::Range(-12, 12), blitz::Range(-6, 6)); AList = 0.0;
        
        Energy(parameter);
        InitStatesArray(parameter);
        RecycleFNeighbors(parameter);

       
    };

    void CountNeighbors(int nm,  int i, int  j, int  k, std::vector<int> count1, std::vector<int> count2) override;
    //void Neighbor(SimulationParameters& parameter) override;

};
// 声明一个多空位BCC晶格类MultiBcc，继承自基类Base
class MultiBcc : public Base {
public:
   
    blitz::Array<int, 5>NN1, NN2; //分别是 5 维数组，用于存储邻居信息
    SimulationParameters *p; //指向 SimulationParameters 类型的指针。

    MultiBcc(SimulationParameters& parameter) {
        p = &parameter;
        Periodic(parameter);
        n1nbr = 8, n2nbr = 6, Length = 2;
        Site.resize(Length, parameter.nx, parameter.ny, parameter.nz); Site = 0;
        NN1.resize(3,Length, parameter.nx, parameter.ny, parameter.nz); NN1 = 0;
        NN2.resize(3,Length, parameter.nx, parameter.ny, parameter.nz); NN2 = 0;
        InitStatesArray(parameter);
        RecycleBNeighbors(parameter);
        
    };
    //覆盖了基类中的 CountNeighbors() 函数，用于计算当前空位的邻居信息。
    void CountNeighbors(int nm,  int i, int  j, int  k, std::vector<int> count1, std::vector<int> count2) override;

};
// 声明一个多空位BCC晶格类MultiFcc，继承自基类Base
class MultiFcc : public Base {
public:
    blitz::Array<int, 5>NN1, NN2;//分别是 5 维数组，用于存储邻居信息
    SimulationParameters* p;//指向 SimulationParameters 类型的指针
    MultiFcc(SimulationParameters& parameter) {
        p = &parameter;
        Periodic(parameter);
        n1nbr = 12, n2nbr = 6, Length = 4;
        Site.resize(Length, parameter.nx, parameter.ny, parameter.nz); Site = 0;
        NN1.resize(3,Length, parameter.nx, parameter.ny, parameter.nz); NN1 = 0;
        NN2.resize(3,Length, parameter.nx, parameter.ny, parameter.nz); NN2 = 0;
        InitStatesArray(parameter);
        RecycleFNeighbors(parameter);
       
        
    }
    //覆盖了基类中的 CountNeighbors() 函数，用于计算当前空位的邻居信息。
    void CountNeighbors(int nm,int i, int  j, int  k, std::vector<int> count1, std::vector<int> count2) override;
};
// 声明一个多类型BCC晶格类（ABVCD 5种类型）MultiSizeBcc，继承自基类Base
class MultiSizeBcc : public Base {
public:
    blitz::Array<int, 5>NN1, NN2;
    SimulationParameters* p;

    MultiSizeBcc(SimulationParameters& parameter) {
        p = &parameter;
        Periodic(parameter);
        n1nbr = 8, n2nbr = 6, Length = 2;
        Site.resize(Length, parameter.nx, parameter.ny, parameter.nz); Site = 0;
        NN1.resize(5, Length, parameter.nx, parameter.ny, parameter.nz); NN1 = 0;
        NN2.resize(5, Length, parameter.nx, parameter.ny, parameter.nz); NN2 = 0;
        InitStatesArray(parameter);
        RecycleBNeighbors(parameter);
    
    };

    void CountNeighbors(int nm, int i, int  j, int  k, std::vector<int> count1, std::vector<int> count2) override;
};
// 声明一个多类型FCC晶格类（ABVCD 5种类型）MultiSizeFcc，继承自基类Base
class MultiSizeFcc : public Base {
public:
    blitz::Array<int, 5>NN1, NN2;
    SimulationParameters* p;

    MultiSizeFcc(SimulationParameters& parameter) {
        p = &parameter;
        Periodic(parameter);
        n1nbr = 12, n2nbr = 6, Length = 4;
        Site.resize(Length, parameter.nx, parameter.ny, parameter.nz); Site = 0;
        NN1.resize(5, Length, parameter.nx, parameter.ny, parameter.nz); NN1 = 0;
        NN2.resize(5, Length, parameter.nx, parameter.ny, parameter.nz); NN2 = 0;
        InitStatesArray(parameter);
        RecycleFNeighbors(parameter);
        
    };

    void CountNeighbors(int nm,  int i, int  j, int  k, std::vector<int> count1, std::vector<int> count2) override;
};
// 声明一个并行多类型BCC晶格类（ABVCD 5种类型）ParaBcc ，继承自基类Base
class ParaBcc : public Base {
public:
   
    blitz::Array<int, 5>NN1, NN2;
    SimulationParameters* p;
    
    ParaBcc(SimulationParameters& parameter) {
        p = &parameter;
        Periodic(parameter);
        n1nbr = 8, n2nbr = 6, Length = 2;
        Site.resize(Length, parameter.nx, parameter.ny, parameter.nz); Site = 0;
        NN1.resize(5, Length, parameter.nx, parameter.ny, parameter.nz); NN1 = 0;
        NN2.resize(5, Length, parameter.nx, parameter.ny, parameter.nz); NN2 = 0;
        Vaclists.resize(parameter.parallelx, parameter.parallely, parameter.parallelz); 
        for (int i = 0; i < parameter.parallelx; i++) {
            for (int j = 0; j < parameter.parallely; j++) {
                for (int k = 0; k < parameter.parallelz; k++) {
                    Vaclists(i, j, k) = std::vector<VAC>(); // 创建一个空的 VAC vector，并将其赋值给 Vaclists 中的元素
                }
            }
        }
        InitStatesArray(parameter);
        RecycleBNeighbors(parameter);

        for (int i = 0; i < Vaclists.extent(0); i++) {
            for (int j = 0; j < Vaclists.extent(1); j++) {
                for (int k = 0; k < Vaclists.extent(2); k++) {
                    std::vector<VAC>& current_vaclist = Vaclists(i, j, k);
                    for (const auto& vac : current_vaclist) {
                        // 在此处对每个 vacancy 做操作
                        std::cout << vac.nm<< "  " << vac.x << "  " << vac.y << "  " << vac.z<<std::endl;
                    }
                }
            }
        }
    };

    void CountNeighbors(int nm, int i, int  j, int  k, std::vector<int> count1, std::vector<int> count2) override;
};


class TypeBuilder {
public:
    // 定义一个 CreatorFunc 类型，指向一个接受 SimulationParameters 引用参数并返回 unique_ptr<Base> 类型的函数
    using CreatorFunc = std::unique_ptr<Base>(*)(SimulationParameters&);

    // 创建指定参数类型的对象
    static std::unique_ptr<Base> create(SimulationParameters& parameter) {
        // 定义一个类型到函数指针的映射表
        static const std::map<std::string, CreatorFunc> typeMap{
            {"BCC", &createBcc},
            {"FCC", &createFcc}
        };

        // 如果参数指定为并行模式，则创建 ParaBcc 对象并返回
        if (parameter.Paralle) {
            return std::make_unique<ParaBcc>(parameter);
        }

        // 根据参数类型查找对应的函数指针
        auto iter = typeMap.find(parameter.par_ltc);
        if (iter == typeMap.end()) {
            return nullptr;
        }

        // 调用对应的函数创建对象并返回
        return iter->second(parameter);
    }

private:
    // 创建 BCC 类型的对象
    static std::unique_ptr<Base> createBcc(SimulationParameters& parameter) {
        if (parameter.MutipleSize) {
            return std::make_unique<MultiSizeBcc>(parameter);
        }
        else if (parameter.Mutiple) {
            return std::make_unique<MultiBcc>(parameter);
        }
        else {
            return std::make_unique<SingleBcc>(parameter);
        }
    }

    // 创建 FCC 类型的对象
    static std::unique_ptr<Base> createFcc(SimulationParameters& parameter) {
        if (parameter.MutipleSize) {
            return std::make_unique<MultiSizeFcc>(parameter);
        }
        else if (parameter.Mutiple) {
            return std::make_unique<MultiFcc>(parameter);
        }
        else {
            return std::make_unique<SingleFcc>(parameter);
        }
    }
};



#endif