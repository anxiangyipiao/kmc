#ifndef _KMC_INIT_H_
#define _KMC_INIT_H_
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "kmc_input.h"
#include "blitz/array.h"


//�������
class Base {
public:
    int n1nbr;                          // ��һ������
    int n2nbr;                          // �ڶ�������
    int Length;                         // ����
    int Bnum;                           // ������Bԭ����
    int Vnum;                           // ������Vԭ����
    int Anum;                           // ������Aԭ����
    int Cnum;                           // ������Cԭ����
    int Dnum;                           // ������Dԭ����
    int vacx;                           // ��λ����x
    int vacy;                           // ��λ����y
    int vacz;                           // ��λ����z
    int vacnm;                          // ��λ���
    blitz::Array<int, 4> Site;          // ���������
    blitz::Array<double, 2> BList;      // B�����б�
    blitz::Array<double, 2> AList;      // A�����б�
    std::vector<int>modla;              // ģ������la
    std::vector<int>modlb;              // ģ������lb
    std::vector<int>modlc;              // ģ������lc
    int(*modlx);                        // ģ������lx
    int(*modly);                        // ģ������ly
    int(*modlz);                        // ģ������lz
    typedef struct {                    // �����λ�ṹ��
        int x;                          // ��λ����x
        int y;                          // ��λ����y
        int z;                          // ��λ����z
        int nm;                         // ��λ���
    } VAC;

    blitz::Array<std::vector<VAC>, 3> Vaclists;    // ��λ�б�
    //SimulationParameters p;//ָ�� SimulationParameters ���͵�ָ��

    virtual ~Base() = default;           // ��������
    void Energy(SimulationParameters& parameter); // ��������
    void InitStatesArray(SimulationParameters& parameter); // ��ʼ������
    void Periodic(SimulationParameters& parameter); // ���������Ա߽�����
    void RecycleFNeighbors(SimulationParameters& parameter); // ����FCC����
    void RecycleBNeighbors(SimulationParameters& parameter); // ����BCC����

    virtual void CountNeighbors(int nm, int i, int  j, int  k, std::vector<int> count1, std::vector<int> count2) = 0; // ���������

private:
    std::vector<int> CalculateN1Neighbors(SimulationParameters& parameter, int nm, int nmn, int i, int  j, int  k);
    std::vector<int> CalculateN1FNeighbors(SimulationParameters& parameter, int nm, int nmn, int i, int  j, int  k);
    std::vector<int> CalculateN2Neighbors(SimulationParameters& parameter, int nm, int nmn, int i, int  j, int  k);
    void CalculateAtomNumbers(SimulationParameters& parameter); // ����ԭ������
    void GenerateLatticePoints(SimulationParameters& parameter); // ���ɾ����
    void AddAtom(SimulationParameters& parameter, int type, int ii, int jj, int zz, int nm); // ���ԭ��
    void TofileStatesArray(SimulationParameters& parameter); // ������״̬д���ļ�
    void AbInitStatesArrays(SimulationParameters& parameter); // ��ʼ������״̬
    std::vector<std::string> Split(const std::string& s, const char& delim);

};

// ����һ������λBCC������SingleBcc���̳��Ի���Base
class SingleBcc : public Base {
public:
    blitz::Array<int, 4> NN1, NN2; // 4ά�������飬�洢ÿ��λ�õ�һ���Ͷ����ھ���
    SimulationParameters* p;
    // SingleBcc��Ĺ��캯�����������ΪSimulationParameters������
    SingleBcc(SimulationParameters& parameter){
        p = &parameter;
        Periodic(parameter); // ִ�������Ա߽������ĳ�ʼ��
        n1nbr = 8, n2nbr = 6, Length = 2;// һ���ھ����������ھ��������񳤶ȵĸ�ֵ
        // ���ݾ��񳤶Ⱥ�ģ�������xyz�ߴ��ʼ��Site���飬������ֵΪ0
        Site.resize(Length, parameter.nx, parameter.ny, parameter.nz); Site = 0;
        NN1.resize(Length, parameter.nx, parameter.ny, parameter.nz); NN1 = 0;
        NN2.resize(Length, parameter.nx, parameter.ny, parameter.nz); NN2 = 0;
        // ���ݸ�����Χ��ʼ��BList���飬������ֵΪ0.0
        BList.resize(blitz::Range(-8, 8), blitz::Range(-6, 6)); BList = 0.0;
        AList.resize(blitz::Range(-8, 8), blitz::Range(-6, 6));AList = 0.0;

        Energy(parameter); // ��������
        InitStatesArray(parameter); // ��ʼ������״̬����
        RecycleBNeighbors(parameter); // ����ÿ��λ�õ�һ���Ͷ����ھ���
        
    }
    // ��дCountNeighbors����������ÿ��λ�õ�һ���Ͷ����ھ���
    void CountNeighbors(int nm, int i, int  j, int  k, std::vector<int> count1, std::vector<int> count2) override;
 
};
// ����һ������λFCC������SingleFcc���̳��Ի���Base  ͬBCC
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
// ����һ�����λBCC������MultiBcc���̳��Ի���Base
class MultiBcc : public Base {
public:
   
    blitz::Array<int, 5>NN1, NN2; //�ֱ��� 5 ά���飬���ڴ洢�ھ���Ϣ
    SimulationParameters *p; //ָ�� SimulationParameters ���͵�ָ�롣

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
    //�����˻����е� CountNeighbors() ���������ڼ��㵱ǰ��λ���ھ���Ϣ��
    void CountNeighbors(int nm,  int i, int  j, int  k, std::vector<int> count1, std::vector<int> count2) override;

};
// ����һ�����λBCC������MultiFcc���̳��Ի���Base
class MultiFcc : public Base {
public:
    blitz::Array<int, 5>NN1, NN2;//�ֱ��� 5 ά���飬���ڴ洢�ھ���Ϣ
    SimulationParameters* p;//ָ�� SimulationParameters ���͵�ָ��
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
    //�����˻����е� CountNeighbors() ���������ڼ��㵱ǰ��λ���ھ���Ϣ��
    void CountNeighbors(int nm,int i, int  j, int  k, std::vector<int> count1, std::vector<int> count2) override;
};
// ����һ��������BCC�����ࣨABVCD 5�����ͣ�MultiSizeBcc���̳��Ի���Base
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
// ����һ��������FCC�����ࣨABVCD 5�����ͣ�MultiSizeFcc���̳��Ի���Base
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
// ����һ�����ж�����BCC�����ࣨABVCD 5�����ͣ�ParaBcc ���̳��Ի���Base
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
                    Vaclists(i, j, k) = std::vector<VAC>(); // ����һ���յ� VAC vector�������丳ֵ�� Vaclists �е�Ԫ��
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
                        // �ڴ˴���ÿ�� vacancy ������
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
    // ����һ�� CreatorFunc ���ͣ�ָ��һ������ SimulationParameters ���ò��������� unique_ptr<Base> ���͵ĺ���
    using CreatorFunc = std::unique_ptr<Base>(*)(SimulationParameters&);

    // ����ָ���������͵Ķ���
    static std::unique_ptr<Base> create(SimulationParameters& parameter) {
        // ����һ�����͵�����ָ���ӳ���
        static const std::map<std::string, CreatorFunc> typeMap{
            {"BCC", &createBcc},
            {"FCC", &createFcc}
        };

        // �������ָ��Ϊ����ģʽ���򴴽� ParaBcc ���󲢷���
        if (parameter.Paralle) {
            return std::make_unique<ParaBcc>(parameter);
        }

        // ���ݲ������Ͳ��Ҷ�Ӧ�ĺ���ָ��
        auto iter = typeMap.find(parameter.par_ltc);
        if (iter == typeMap.end()) {
            return nullptr;
        }

        // ���ö�Ӧ�ĺ����������󲢷���
        return iter->second(parameter);
    }

private:
    // ���� BCC ���͵Ķ���
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

    // ���� FCC ���͵Ķ���
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