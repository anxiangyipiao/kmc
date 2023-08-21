#pragma once
#ifndef _KMC_JUMP_H_
#define _KMC_JUMP_H_
#include "kmc_init.h"
#include <algorithm>


class JumpBase {
public:

    int path;
    blitz::Array<int, 1> x, y, z; blitz::Array<double, 1> w;
    virtual ~JumpBase() = default;
    // 计算能量的虚函数，需要在派生类中进行实现
    //CalculatedEnergy
    // BCC 类型
    template<typename T>
    void CalculatedBEnergy(int& nf, T* s);
    // FCC 类型
    template<typename T>
    void CalculatedFEnergy(blitz::Array<int, 1>& fnf,T* s);

    virtual void CalculatedEnergy(int inv, int nm, int inv0) = 0;

    std::vector<int> ResponseMultilist1(std::vector<int> listBase) {

        std::vector<int> List;
        int i = listBase[3] - 1 - listBase[0]; int j = listBase[4] - listBase[1];
        int k = listBase[5] + listBase[0] + 2 - listBase[2] - listBase[3]; int n = listBase[1] - listBase[4]; int m = listBase[2] - listBase[5] - 1;
        int ii = listBase[9] - listBase[6]; int jj = listBase[10] - listBase[7]; int kk = listBase[11] + listBase[6] - listBase[8] - listBase[9];
        int nn = listBase[7] - listBase[10]; int mm = listBase[8] - listBase[11];

        List = {i,j,k,n,m,ii,jj,kk,nn,mm};

        return List;
    };

    std::vector<int> ResponseMultilist2(std::vector<int> listBase) {
       
        std::vector<int> List;

        int i = listBase[4] - 1 - listBase[1]; int j = listBase[3] - listBase[0];
        int k = listBase[5] + listBase[1] + 2 - listBase[4] - listBase[2]; int n = listBase[0] - listBase[3]; int m = listBase[2] - listBase[5] - 1;
        int ii = listBase[10] - listBase[7]; int jj = listBase[9] - listBase[6];
        int kk = listBase[11] + listBase[7]- listBase[8] - listBase[10]; int nn = listBase[6] - listBase[9]; int mm = listBase[8] - listBase[11];
       
        List = { i,j,k,n,m,ii,jj,kk,nn,mm };

        return List;
    };


    //SumAndChoice
    template<typename T>
    void SumAndChoice(T* s);
    template<typename T>
    void SumAndChoice2(T* s);
    
    //Neighbor
    template<typename T>
    void UpdateMutipleN1BeforeNeighbor(T* s, int x, int y, int z, int nf, int vacnm, int type);
    template<typename T>
    void UpdateMutipleN1AfterNeighbor(T* s, int x, int y, int z, int nf,int vacnm,int type);
    
    template<typename T>
    void UpdateMutipleFN1BeforeNeighbor(T* s, int x, int y, int z, int fnf, int type);
    template<typename T>
    void UpdateMutipleFN1AfterNeighbor(T* s, int x, int y, int z, int vacnm, int type);
    
    template<typename T>
    void UpdateMutipleN2BeforeNeighbor(T* s, int x, int y, int z, int nf, int type);
    template<typename T>
    void UpdateMutipleN2AfterNeighbor(T* s, int x, int y, int z, int nf, int type);
    
    template<typename T>
    void UpdateSingleN2BeforeNeighbor(T* s, int x, int y, int z, int nf);
    template<typename T>
    void UpdateSingleN2AafterNeighbor(T* s, int x, int y, int z, int nf);
    

    virtual void UpdateNeighbor() = 0;
    //ExchangeSite
    template<typename T>
    void ExchangeSite(T* s,int nf);
   
    
    //run type(BCC FCC）
    template<typename T>
    void runBSimulationLoop(int& nf, T* s);

    template<typename T>
    void runFSimulationLoop(blitz::Array<int, 1>& fnf, T* s); 
   
    //run 
    virtual void runSimulationLoop() = 0;
};


class BccJump :public JumpBase {
public:
    BccJump(SingleBcc* s) {
        x.resize(8); y.resize(8); z.resize(8); w.resize(8);
        x = 0; y = 0; z = 0; w = 0; nf = 0;
        this->s = s;
    }
   
    void CalculatedEnergy(int inv, int nm, int inv0) override;
    void UpdateNeighbor() override;
    void runSimulationLoop() override;
private:
    SingleBcc* s;
    int nf;
};


class FccJump :public JumpBase {
public:
    FccJump(SingleFcc* s) {
        x.resize(12); y.resize(12); z.resize(12); fnf.resize(12); w.resize(12);
        x = 0; y = 0; z = 0; w = 0; fnf = 0;
        this->s = s;
    }
    void CalculatedEnergy(int inv, int nm, int inv0) override;
    void UpdateNeighbor() override;
    void runSimulationLoop() override;
private:
    SingleFcc* s;
    blitz::Array<int, 1> fnf;

};


class MultiBccJump :public JumpBase {
public:
    MultiBccJump(MultiBcc* s) {
        x.resize(8); y.resize(8); z.resize(8); w.resize(8);
        x = 0; y = 0; z = 0; w = 0; nf = 0;
        this->s = s;
    }
    std::vector<int> ReturnList(int inv);
    void CalculatedEnergy(int inv, int nm, int inv0) override;
    void UpdateNeighbor() override;
    void runSimulationLoop() override ;
private:
    MultiBcc* s;
    int nf;

};


class MultiFccJump :public JumpBase {
public:
    MultiFccJump(MultiFcc* s) {
        x.resize(12); y.resize(12); z.resize(12); fnf.resize(12);
        x = 0; y = 0; z = 0; w.resize(12); w = 0; fnf = 0;
        this->s = s;
    }
    std::vector<int> ReturnList(int inv,int nm);
    void CalculatedEnergy(int inv, int nm, int inv0) override;
    void UpdateNeighbor() override;
    void runSimulationLoop() override;

private:
    MultiFcc* s;
    blitz::Array<int, 1> fnf;

};


class MultiSizeBccJump :public JumpBase {
public:
    MultiSizeBccJump(MultiSizeBcc* s) {
        x.resize(8); y.resize(8); z.resize(8); w.resize(8);
        x = 0; y = 0; z = 0; w = 0; nf = 0;
        this->s = s;
    }

    void CalculatedEnergy(int inv, int nm, int inv0) override;
    void UpdateNeighbor() override;
    void runSimulationLoop() override;
private:
    MultiSizeBcc* s;
    int nf;
};


class MultiSizeFccJump :public JumpBase {
public:
    MultiSizeFccJump(MultiSizeFcc* s) {
        x.resize(12); y.resize(12); z.resize(12); fnf.resize(12);
        x = 0; y = 0; z = 0; w.resize(12); w = 0; fnf = 0;
        this->s = s;
    }

    void CalculatedEnergy(int inv, int nm, int inv0) override;
    void UpdateNeighbor() override;
    void runSimulationLoop() override;

private:
    MultiSizeFcc* s;
    blitz::Array<int, 1> fnf;

};

// 并行策略  分核运行 如果一个区域内没有空位，那么就跳过，每个区域运行固定时间片
class ParaBccJump :public JumpBase {
public:
    ParaBccJump(ParaBcc* s) {
        this->s = s;
        Divide();
    }
    void runSimulationLoop() override {};
    void UpdateNeighbor() override {};
    void CalculatedEnergy(int inv, int nm, int inv0) override {};

    void Divide();
    void Sector(int i,int j,int k,int idx,int idy,int idz);
    void CalculatedEnergys(int vacx, int vacy, int vacz, int vacnm, int& nf,
        blitz::Array<int, 1> x, blitz::Array<int, 1> y,
        blitz::Array<int, 1> z, blitz::Array<double, 1> w);
    void Choice(int& path, blitz::Array<double, 1> w);
    void UpdateNeighbors(int vacx, int vacy, int vacz, int vacnm, int& nf,
        blitz::Array<int, 1> x, blitz::Array<int, 1> y,
        blitz::Array<int, 1> z, blitz::Array<double, 1> w,int& path);
    void ExchangeSites(int& vacx, int& vacy, int& vacz, int& vacnm, int& nf,
        blitz::Array<int, 1> x, blitz::Array<int, 1> y,
        blitz::Array<int, 1> z, blitz::Array<double, 1> w, int& path);

    bool Border(int& vacx, int& vacy, int& vacz, int& vacnm,
        int i,int j,int k, int idx, int idy, int idz);
private:
    ParaBcc* s;
};




class JumpBuilder {
public:
    static std::unique_ptr<JumpBase> create(std::unique_ptr<Base>& obj) {
        if (SingleBcc* sb = dynamic_cast<SingleBcc*>(obj.get())) {
            // obj 指向 SingleBcc 实例
            return std::make_unique<BccJump>(sb);
        }
        else if (SingleFcc* sf = dynamic_cast<SingleFcc*>(obj.get())) {
            // obj 指向 SingleFcc 实例
            return std::make_unique<FccJump>(sf);
        }
        else if (MultiBcc* mb = dynamic_cast<MultiBcc*>(obj.get())) {
            // obj 指向 MultiBcc 实例
            return std::make_unique<MultiBccJump>(mb);
        }
        else if (MultiFcc* mf = dynamic_cast<MultiFcc*>(obj.get())) {
            // obj 指向 MultiFcc 实例
            return std::make_unique<MultiFccJump>(mf);
        }
        else if (MultiSizeBcc* mf = dynamic_cast<MultiSizeBcc*>(obj.get())) {
            // obj 指向 MultiFcc 实例
            return std::make_unique<MultiSizeBccJump>(mf);
        }
        else if (MultiSizeFcc* mf = dynamic_cast<MultiSizeFcc*>(obj.get())) {
            // obj 指向 MultiFcc 实例
            return std::make_unique<MultiSizeFccJump>(mf);
        }
        else if (ParaBcc* mf = dynamic_cast<ParaBcc*>(obj.get())) {
            // obj 指向 MultiFcc 实例
            return std::make_unique<ParaBccJump>(mf);
        }
        else {
            // obj 指向一个未知类型的实例，抛出异常或采取其他适当的措施
            throw std::invalid_argument("Unknown object type");
        }
    }
};





#endif


