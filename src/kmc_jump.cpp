#include "kmc_jump.h"


//JumpBase


/**
 * @brief 计算 BCC 结构下的能量
 *
 * 该函数用于计算 BCC 结构下的能量，并更新传入的 `nf` 参数以反映空位邻居的位置。
 *
 * @tparam T 泛型类型，通常为类的指针，包含用于能量计算的相关数据。
 * @param nf 用于存储是反映空位邻居的位置。
 * @param s 指向包含相关数据的对象的指针。
 */
template<typename T>
void JumpBase::CalculatedBEnergy(int& nf, T* s)
{
	nf = (s->vacnm == 0) ? 1 : 0;
	for (int inv = 0; inv < 8; inv++) {
	
		x(inv) = (s->modlx[s->vacx + v1nbr_bcc[s->vacnm][0][inv]]);
		y(inv) = (s->modly[s->vacy + v1nbr_bcc[s->vacnm][1][inv]]);
		z(inv) = (s->modlz[s->vacz + v1nbr_bcc[s->vacnm][2][inv]]);
		CalculatedEnergy(inv, 0, 0);
	
	}
}

template<typename T>
void JumpBase::CalculatedFEnergy(blitz::Array<int, 1>& fnf,T* s)
{
	int inv0 = 0;
	for (int nm = 0; nm < 4; nm++) {
		if (nm != s->vacnm) {
			for (int inv = 0; inv < 4; inv++) {

				x(inv0) = s->modlx[s->vacx + v1nbr_fcc[0][s->vacnm][nm][inv]];
				y(inv0) = s->modly[s->vacy + v1nbr_fcc[1][s->vacnm][nm][inv]];
				z(inv0) = s->modlz[s->vacz + v1nbr_fcc[2][s->vacnm][nm][inv]];
				fnf(inv0) = (nm);
				CalculatedEnergy(inv,nm,inv0);
				inv0++;
			}
		}
	}
}

/**
 * @brief 计算权重总和并根据随机值选择路径       direct KMC
 *
 * 该函数用于计算权重数组中所有元素的总和，并根据生成的随机数选择路径。
 *
 * @tparam T 泛型类型，通常为类的指针，包含用于权重计算的相关数据。
 * @param s 指向包含相关数据的对象的指针。
 */
template<typename T>
void JumpBase::SumAndChoice(T* s)
{
	double sum = 0;
	// 计算权重数组中所有元素的总和
	for (int i = 0; i < s->n1nbr; i++) {
		sum = sum + w(i);
	}

	// 生成一个 [0, sum) 范围内的随机数
	double ran = u(e) * sum;

	double boundry = 0;
	for (int i = 0; i < s->n1nbr; i++) {
		boundry = boundry + w(i);
		if (ran < boundry) {
			path = i;  // 根据随机数选择路径
			break;
		}
	}
}

/**
 * @brief 计算权重数组中的最大权重，并根据随机值选择路径      reject KMC
 *
 * 该函数用于计算权重数组中的最大权重，并根据生成的随机数选择路径。
 *
 * @tparam T 泛型类型，通常为类的指针，包含用于权重计算的相关数据。
 * @param s 指向包含相关数据的对象的指针。
 */
template<typename T>
void JumpBase::SumAndChoice2(T* s)
{
	// 找到权重数组中的最大权重
	double max_w = w(0);
	for (int i = 1; i < s->n1nbr; i++) {
		if (w(i) > max_w) {
			max_w = w(i);
		}
	}

	double r = u(e) * s->n1nbr;

	int j = r + 1;

	for (int i = 0; i < s->n1nbr; i++) {
		// 根据随机数和最大权重选择路径
		if (j - r < (w(i) / max_w)) {
			path = j - 1;
			break;
		}
	}

}


/**
 * @brief 更新单空位跳跃之前的第二近邻(NN2)数量（邻居数量）。
 *
 * 该函数用于在进行跳跃之前，更新单一晶格点的第二近邻(NN2)数量，以反映跳跃前的状态。
 *
 * @tparam T 泛型类型，通常为类的指针，包含晶格信息和邻居数量的相关数据。
 * @param s 指向包含晶格信息和邻居数量的对象的指针。
 * @param x 晶格点的 x 坐标。
 * @param y 晶格点的 y 坐标。
 * @param z 晶格点的 z 坐标。
 * @param nf 跳跃之前的晶格点类型或状态。
 */
template<typename T>
void JumpBase::UpdateSingleN2BeforeNeighbor(T* s, int x, int y, int z, int nf)
{
	for (int jj = 0; jj < 6; jj++) {
		int xx = s->modlx[x + v2nbr_bcc[0][jj]];
		int	yy = s->modly[y + v2nbr_bcc[1][jj]];
		int	zz = s->modlz[z + v2nbr_bcc[2][jj]];
		s->NN2(nf, xx, yy, zz) = s->NN2(nf, xx, yy, zz) - 1;
	}

}

/**
 * @brief 更新单一晶格点的第二近邻(NN2)数量（邻居数量）。
 *
 * 该函数用于在进行跳跃之后，更新单一晶格点的第二近邻(NN2)数量，以反映跳跃后的状态。
 */
template<typename T>
void JumpBase::UpdateSingleN2AafterNeighbor(T* s, int x, int y, int z, int nf)
{
	for (int jj = 0; jj < 6; jj++) {
		int xx = s->modlx[x + v2nbr_bcc[0][jj]];
		int	yy = s->modly[y + v2nbr_bcc[1][jj]];
		int	zz = s->modlz[z + v2nbr_bcc[2][jj]];
		s->NN2(nf, xx, yy, zz) = s->NN2(nf, xx, yy, zz) + 1;
	}
}

/**
 * @brief 更新多空位晶格点的第一近邻(NN1)数量（邻居数量）。
 *
 * 该函数用于在进行跳跃之前，更新多类型晶格点的第一近邻(NN1)数量，以反映跳跃前的状态。
 *
 * @tparam T 泛型类型，通常为类的指针，包含晶格信息和邻居数量的相关数据。
 * @param s 指向包含晶格信息和邻居数量的对象的指针。
 * @param x 晶格点的 x 坐标。
 * @param y 晶格点的 y 坐标。
 * @param z 晶格点的 z 坐标。
 * @param nf 跳跃前的晶格点类型或状态。
 * @param vacnm 空位的标识号。
 * @param type 要更新的晶格点类型。
 */
template<typename T>
void JumpBase::UpdateMutipleN1BeforeNeighbor(T* s,int x, int y, int z, int nf, int vacnm, int type)
{
	for (int ii = 0; ii < 8; ii++) {
		int xx = s->modlx[x + v1nbr_bcc[nf][0][ii]];
		int yy = s->modly[y + v1nbr_bcc[nf][1][ii]];
		int zz = s->modlz[z + v1nbr_bcc[nf][2][ii]];

		s->NN1(type, vacnm, xx, yy, zz) = s->NN1(type, vacnm, xx, yy, zz) - 1;
		s->NN1(2, vacnm, xx, yy, zz) = s->NN1(2, vacnm, xx, yy, zz) + 1;
	}
}

template<typename T>
void JumpBase::UpdateMutipleN1AfterNeighbor(T* s, int x, int y, int z, int nf, int vacnm, int type)
{
	for (int ii = 0; ii < 8; ii++) {
		int xx = s->modlx[x + v1nbr_bcc[vacnm][0][ii]];
		int yy = s->modly[y + v1nbr_bcc[vacnm][1][ii]];
		int zz = s->modlz[z + v1nbr_bcc[vacnm][2][ii]];

		s->NN1(type, nf, xx, yy, zz) = s->NN1(type, nf, xx, yy, zz) + 1;
		s->NN1(2, nf, xx, yy, zz) = s->NN1(2, nf, xx, yy, zz) - 1;
	}
}

//多空位FCC 第一
template<typename T>
void JumpBase::UpdateMutipleFN1BeforeNeighbor(T* s, int x, int y, int z, int fnf, int type)
{
	for (int ii = 0; ii < 4; ii++) {
		if (ii != fnf) {
			for (int kk = 0; kk < 4; kk++) {
				int xx = s->modlx[x + v1nbr_fcc[0][fnf][ii][kk]];
				int yy = s->modly[y + v1nbr_fcc[1][fnf][ii][kk]];
				int zz = s->modlz[z + v1nbr_fcc[2][fnf][ii][kk]];

				s->NN1(type, ii, xx, yy, zz) = s->NN1(type, ii, xx, yy, zz) - 1;
				s->NN1(2, ii, xx, yy, zz) = s->NN1(2, ii, xx, yy, zz) + 1;
			}
		}
	}
}

template<typename T>
void JumpBase::UpdateMutipleFN1AfterNeighbor(T* s, int x, int y, int z, int vacnm, int type)
{
	for (int ii = 0; ii < 4; ii++) {
		if (ii != vacnm) {
			for (int kk = 0; kk < 4; kk++) {
				int xx = s->modlx[x + v1nbr_fcc[0][vacnm][ii][kk]];
				int yy = s->modly[y + v1nbr_fcc[1][vacnm][ii][kk]];
				int zz = s->modlz[z + v1nbr_fcc[2][vacnm][ii][kk]];

				s->NN1(type, ii, xx, yy, zz) = s->NN1(type, ii, xx, yy, zz) + 1;
				s->NN1(2, ii, xx, yy, zz) = s->NN1(2, ii, xx, yy, zz) - 1;

			}
		}
	}
	
}

//多空位 第二
template<typename T>
void JumpBase::UpdateMutipleN2BeforeNeighbor(T* s, int x, int y, int z, int nf,int type)
{
	for (int jj = 0; jj < 6; jj++) {
		int xx = s->modlx[x + v2nbr_bcc[0][jj]];
		int	yy = s->modly[y + v2nbr_bcc[1][jj]];
		int	zz = s->modlz[z + v2nbr_bcc[2][jj]];

		s->NN2(type, nf, xx, yy, zz) = s->NN2(type, nf, xx, yy, zz) - 1;
		s->NN2(2, nf, xx, yy, zz) = s->NN2(2, nf, xx, yy, zz) + 1;
	}
}

template<typename T>
void JumpBase::UpdateMutipleN2AfterNeighbor(T* s, int x, int y, int z, int vacnm,int type)
{
	for (int jj = 0; jj < 6; jj++) {
		int xx = s->modlx[x + v2nbr_bcc[0][jj]];
		int	yy = s->modly[y + v2nbr_bcc[1][jj]];
		int	zz = s->modlz[z + v2nbr_bcc[2][jj]];

		s->NN2(type, vacnm, xx, yy, zz) = s->NN2(type, vacnm, xx, yy, zz) + 1;
		s->NN2(2, vacnm, xx, yy, zz) = s->NN2(2, vacnm, xx, yy, zz) - 1;

	}
}

/**
 * @brief 交换晶格点上的原子状态。
 *
 * 该函数用于交换两个晶格点上的原子状态，实现跳跃操作。交换后，空位的位置和标识将更新为新的位置和类型。
 *
 * @tparam T 泛型类型，通常为类的指针，包含晶格信息和邻居数量的相关数据。
 * @param s 指向包含晶格信息和邻居数量的对象的指针。
 * @param nf 跳跃前的晶格点类型或状态。
 */
template<typename T>
void JumpBase::ExchangeSite(T* s,int nf)
{
	
	s->Site(s->vacnm, s->vacx, s->vacy, s->vacz) = s->Site(nf, x(path), y(path), z(path));
	s->vacx = x(path); s->vacy = y(path); s->vacz = z(path); s->vacnm = nf;
	s->Site(s->vacnm, s->vacx, s->vacy, s->vacz) = 2;
}

/**
 * @brief 运行BCC结构的晶格跳跃模拟循环。
 *
 * 该函数执行BCC结构晶格的跳跃模拟循环，包括计算能量、选择跳跃路径、更新邻居和交换原子位置。
 *
 * @tparam T 泛型类型，通常为类的指针，包含晶格信息和邻居数量的相关数据。
 * @param nf 用于存储跳跃前的晶格点类型或状态的变量。
 * @param s 指向包含晶格信息和邻居数量的对象的指针。
 */
template<typename T>
void JumpBase::runBSimulationLoop(int& nf, T* s)
{
	// 计算BCC结构的能量，同时更新跳跃前的晶格点类型或状态
		CalculatedBEnergy(nf, s);
		// 选择跳跃路径
		SumAndChoice(s);
		// 更新邻居信息
		UpdateNeighbor();
		// 执行原子位置交换，完成跳跃操作
		ExchangeSite(s,nf);

}

/**
 * @brief 运行FCC结构的晶格跳跃模拟循环。
 *
 * 该函数执行FCC结构晶格的跳跃模拟循环，包括计算能量、选择跳跃路径、更新邻居和交换原子位置。
 *
 * @tparam T 泛型类型，通常为类的指针，包含晶格信息和邻居数量的相关数据。
 * @param fnf 用于存储跳跃前的晶格点类型或状态的数组。
 * @param s 指向包含晶格信息和邻居数量的对象的指针。
 */
template<typename T>
void JumpBase::runFSimulationLoop(blitz::Array<int, 1>& fnf, T* s)
{

		CalculatedFEnergy(fnf, s);
		SumAndChoice(s);
		UpdateNeighbor();
		ExchangeSite(s, fnf(path));
	
}


/**
 * @brief 计算BCC结构晶格点的能量变化。
 *
 * 该函数用于计算BCC结构晶格点跳跃前后的能量变化。
 *
 * @param inv 跳跃路径的索引，用于选择邻居位置。
 * @param nm 跳跃前晶格点的类型或状态。
 * @param inv0 跳跃前晶格点的索引。
 */
void BccJump::CalculatedEnergy(int inv, int nm, int inv0)
{
	int dn1 = s->NN1(s->vacnm, s->vacx, s->vacy, s->vacz) - s->NN1(nf, x(inv), y(inv), z(inv));
	int dn2 = s->NN2(s->vacnm, s->vacx, s->vacy, s->vacz) - s->NN2(nf, x(inv), y(inv), z(inv));

	if (s->Site(nf, x(inv), y(inv), z(inv)) == 1) {

		w(inv) = (s->BList(dn1, dn2));
		
	}
	else {

		w(inv) = (s->AList(dn1, dn2));

	}
	
}

/**
 * @brief 更新BCC结构晶格点的邻居状态。
 *
 * 该函数用于更新BCC结构晶格点跳跃前后的邻居状态。
 *
 * 如果跳跃前晶格点类型为1，即B原子，则执行以下操作：
 * 1. 更新跳跃前晶格点的NN1邻居状态。
 * 2. 调用 UpdateSingleN2BeforeNeighbor 更新跳跃前晶格点的NN2邻居状态。
 * 3. 更新跳跃后晶格点的NN1邻居状态。
 * 4. 调用 UpdateSingleN2AfterNeighbor 更新跳跃后晶格点的NN2邻居状态。
 *
 * @note 该函数假设 BCC 结构中原子类型为 1 表示 B 原子，其他类型表示 A 原子。
 */
void BccJump::UpdateNeighbor()
{
	if (s->Site(nf, x(path), y(path), z(path)) == 1) {
		
		for (int ii = 0; ii < 8; ii++) {
			int xx = s->modlx[x(path) + v1nbr_bcc[nf][0][ii]];
			int yy = s->modly[y(path) + v1nbr_bcc[nf][1][ii]];
			int zz = s->modlz[z(path) + v1nbr_bcc[nf][2][ii]];
			s->NN1(s->vacnm, xx, yy, zz) = s->NN1(s->vacnm, xx, yy, zz) - 1;
		}

		UpdateSingleN2BeforeNeighbor(s, x(path), y(path), z(path), nf);
		
		for (int ii = 0; ii < 8; ii++) {
			int xx = s->modlx[s->vacx + v1nbr_bcc[s->vacnm][0][ii]];
			int yy = s->modly[s->vacy + v1nbr_bcc[s->vacnm][1][ii]];
			int zz = s->modlz[s->vacz + v1nbr_bcc[s->vacnm][2][ii]];
			s->NN1(nf, xx, yy, zz) = s->NN1(nf, xx, yy, zz) + 1;
		}

		UpdateSingleN2AafterNeighbor(s, s->vacx, s->vacy, s->vacz, s->vacnm);
		
	}

}

/**
 * @brief 运行BCC结构跳跃模拟循环。
 *
 * 该函数用于运行BCC结构的跳跃模拟循环，调用了基类中的 runBSimulationLoop 函数。
 */
void BccJump::runSimulationLoop()
{
	BccJump::runBSimulationLoop(nf, s);
}


//FCC
void FccJump::CalculatedEnergy(int inv, int nm, int inv0)
{
	
	int dn1 = s->NN1(s->vacnm, s->vacx, s->vacy, s->vacz) - s->NN1(nm, x(inv), y(inv), z(inv));
	int dn2 = s->NN2(s->vacnm, s->vacx, s->vacy, s->vacz) - s->NN2(nm, x(inv), y(inv), z(inv));

	if (s->Site(nm, x(inv), y(inv), z(inv)) == 1) {
			w(inv0) = s->BList(dn1, dn2);
	}
	else {
			w(inv0) = s->AList(dn1, dn2);
	}

}

void FccJump::UpdateNeighbor()
{
	if (s->Site(fnf(path), x(path), y(path), z(path)) == 1) {

		for (int ii = 0; ii < 4; ii++) {
			if (ii != fnf(path)) {
				for (int kk = 0; kk < 4; kk++) {
					int xx = s->modlx[x(path) + v1nbr_fcc[0][fnf(path)][ii][kk]];
					int yy = s->modly[y(path) + v1nbr_fcc[1][fnf(path)][ii][kk]];
					int zz = s->modlz[z(path) + v1nbr_fcc[2][fnf(path)][ii][kk]];
					s->NN1(ii, xx, yy, zz) = s->NN1(ii, xx, yy, zz) - 1;
				}
			}
		}

		UpdateSingleN2BeforeNeighbor(s, x(path), y(path), z(path), fnf(path));
	
		for (int ii = 0; ii < 4; ii++) {
			if (ii != s->vacnm) {
				for (int kk = 0; kk < 4; kk++) {
					int xx = s->modlx[s->vacx + v1nbr_fcc[0][s->vacnm][ii][kk]];
					int yy = s->modly[s->vacy + v1nbr_fcc[1][s->vacnm][ii][kk]];
					int zz = s->modlz[s->vacz + v1nbr_fcc[2][s->vacnm][ii][kk]];
					s->NN1(ii, xx, yy, zz) = s->NN1(ii, xx, yy, zz) + 1;
				}
			}
		}

		UpdateSingleN2BeforeNeighbor(s, s->vacx, s->vacy, s->vacz, s->vacnm);
	
	}
}

void FccJump::runSimulationLoop()
{
	FccJump::runFSimulationLoop(fnf, s);
}

/**
 * @brief 计算并返回一个包含多个整数值的列表。
 *
 * 这个函数计算了一组整数值，表示了特定位置和邻居的状态信息，包括两种不同的原子类型（1和0）以及相邻位置的NN1和NN2的值。
 * 返回的列表依次包括：NN1(1), NN1(0), NN1(2), NN1(1)（空位的NN1(1)）, NN1(0)（空位的NN1(0)）, NN1(2)（空位的NN1(2)）,
 * NN2(1), NN2(0), NN2(2), NN2(1)（空位的NN2(1)）, NN2(0)（空位的NN2(0)）, NN2(2)（空位的NN2(2)）。
 *
 * @param inv 用于索引位置和邻居状态的索引。
 * @return 包含计算的整数值的列表。
 */
std::vector<int> MultiBccJump::ReturnList(int inv)
{
	int a = s->NN1(1, nf, x(inv), y(inv), z(inv)); int b = s->NN1(0, nf, x(inv), y(inv), z(inv)); int c = s->NN1(2, nf, x(inv), y(inv), z(inv));
	int d = s->NN1(1, s->vacnm, s->vacx, s->vacy, s->vacz); int e = s->NN1(0, s->vacnm, s->vacx, s->vacy, s->vacz); int f = s->NN1(2, s->vacnm, s->vacx, s->vacy, s->vacz);
	
	int aa = s->NN2(1, nf, x(inv), y(inv), z(inv)); int bb = s->NN2(0, nf, x(inv), y(inv), z(inv)); int cc = s->NN2(2, nf, x(inv), y(inv), z(inv));
	int dd = s->NN2(1, s->vacnm, s->vacx, s->vacy, s->vacz); int ee = s->NN2(0, s->vacnm, s->vacx, s->vacy, s->vacz); int ff = s->NN2(2, s->vacnm, s->vacx, s->vacy, s->vacz);

	std::vector<int> list = { a,b,c,d,e,f,aa,bb,cc,dd,ee,ff };
	return list;
}

/**
 * @brief 计算并设置跳跃能量和权重。
 *
 * 这个函数根据不同的情况计算跳跃能量，并为跳跃设置权重。
 * 如果当前位置为1，计算 tde_vb 并计算权重。
 * 如果当前位置为0，计算 tde_va 并计算权重。
 * 如果当前位置不是1或0，则权重为0。
 *
 * @param inv 用于索引位置和邻居状态的索引。
 * @param nm 未使用的参数。
 * @param inv0 未使用的参数。
 */
void MultiBccJump::CalculatedEnergy(int inv, int nm, int inv0)
{
	
	if (s->Site(nf, x(inv), y(inv), z(inv)) == 1)
		{
			
			std::vector<int> listBase = ReturnList(inv);

			std::vector<int> list = ResponseMultilist1(listBase);

			double tde_vb = list[0] * s->p->par_eSPB1B + list[1] * s->p->par_eSPA1B + list[2] * s->p->par_eSPB1V + list[3] * s->p->par_eSPA1V + list[4] * s->p->par_eSPV1V +
				list[5] * s->p->par_eSPB2B + list[6] * s->p->par_eSPA2B + list[7] * s->p->par_eSPB2V + list[8] * s->p->par_eSPA2V + list[9] * s->p->par_eSPV2V;

			w(inv) = exp(-(tde_vb / 2 + s->p->par_eSPB) / (s->p->par_beta * s->p->par_temp));

		}
	else if (s->Site(nf, x(inv), y(inv), z(inv)) == 0)
		{
			
			std::vector<int> listBase = ReturnList(inv);

			std::vector<int> list = ResponseMultilist2(listBase);

			double tde_va = list[0] * s->p->par_eSPA1A + list[1] * s->p->par_eSPA1B + list[2] * s->p->par_eSPA1V + list[3] * s->p->par_eSPB1V + list[4] * s->p->par_eSPV1V +
				list[5] * s->p->par_eSPA2A + list[6] * s->p->par_eSPA2B + list[7] * s->p->par_eSPA2V + list[8] * s->p->par_eSPB2V + list[9] * s->p->par_eSPV2V;

			w(inv) = exp(-(tde_va / 2 + s->p->par_eSPA) / (s->p->par_beta * s->p->par_temp));

		}
	else {
		w(inv) = 0.0;
	}

}

/**
 * @brief 根据当前位置的状态更新邻居信息。
 *
 * 如果当前位置为0（Fe），则使用0来更新邻居信息。
 * 如果当前位置为1（Cu），则使用1来更新邻居信息。
 * 更新过程中，会先更新N1邻居，然后更新N2邻居。
 *
 * @note 在这个函数中，通过调用不同的更新函数来更新N1和N2邻居。
 *
 */
void MultiBccJump::UpdateNeighbor()
{
	if (s->Site(nf, x(path), y(path), z(path)) == 0) {
		
		// 0 - Fe    NN1 里是 0 

		UpdateMutipleN1BeforeNeighbor(s, x(path), y(path), z(path), nf, s->vacnm, 0);
	
		UpdateMutipleN2BeforeNeighbor(s, x(path), y(path), z(path), nf, 0);
	
		UpdateMutipleN1AfterNeighbor(s, s->vacx, s->vacy, s->vacz, nf, s->vacnm, 0);

		UpdateMutipleN2AfterNeighbor(s, s->vacx, s->vacy, s->vacz, s->vacnm, 0);

	}
	else if (s->Site(nf, x(path), y(path), z(path)) == 1) {
		//1 - Cu    NN储存在1位置

		UpdateMutipleN1BeforeNeighbor(s, x(path), y(path), z(path), nf, s->vacnm, 1);

		UpdateMutipleN2BeforeNeighbor(s, x(path), y(path), z(path), nf, 1);

		UpdateMutipleN1AfterNeighbor(s, s->vacx, s->vacy, s->vacz, nf, s->vacnm, 1);
	
		UpdateMutipleN2AfterNeighbor(s, s->vacx, s->vacy, s->vacz, s->vacnm, 1);

	}
}

void MultiBccJump::runSimulationLoop()
{
	MultiBccJump::runBSimulationLoop(nf, s);
}

//MultiFcc
std::vector<int> MultiFccJump::ReturnList(int inv,int nm)
{
	int a = s->NN1(1, nm, x(inv), y(inv), z(inv)); int b = s->NN1(0, nm, x(inv), y(inv), z(inv)); int c = s->NN1(2, nm, x(inv), y(inv), z(inv));
	int d = s->NN1(1, s->vacnm, s->vacx, s->vacy, s->vacz); int e = s->NN1(0, s->vacnm, s->vacx, s->vacy, s->vacz); int f = s->NN1(2, s->vacnm, s->vacx, s->vacy, s->vacz);

	int aa = s->NN2(1, nm, x(inv), y(inv), z(inv)); int bb = s->NN2(0, nm, x(inv), y(inv), z(inv)); int cc = s->NN2(2, nm, x(inv), y(inv), z(inv));
	int dd = s->NN2(1, s->vacnm, s->vacx, s->vacy, s->vacz); int ee = s->NN2(0, s->vacnm, s->vacx, s->vacy, s->vacz); int ff = s->NN2(2, s->vacnm, s->vacx, s->vacy, s->vacz);

	std::vector<int> list = { a,b,c,d,e,f,aa,bb,cc,dd,ee,ff };
	return list;
}

void MultiFccJump::CalculatedEnergy(int inv, int nm, int inv0)
{
		if (s->Site(nm, x(inv), y(inv), z(inv)) == 1)
		{
					
				std::vector<int> listBase = ReturnList(inv,nm);

				std::vector<int> list = ResponseMultilist1(listBase);

				double tde_vb = list[0] * s->p->par_eSPB1B + list[1] * s->p->par_eSPA1B + list[2] * s->p->par_eSPB1V + list[3] * s->p->par_eSPA1V + list[4] * s->p->par_eSPV1V +
						list[5] * s->p->par_eSPB2B + list[6] * s->p->par_eSPA2B + list[7] * s->p->par_eSPB2V + list[8] * s->p->par_eSPA2V + list[9] * s->p->par_eSPV2V;

				w(inv0) = exp(-(tde_vb / 2 + s->p->par_eSPB) / (s->p->par_beta * s->p->par_temp));

		}
		else if (s->Site(nm, x(inv), y(inv), z(inv)) == 0)
		{
					

				std::vector<int> listBase = ReturnList(inv, nm);

				std::vector<int> list = ResponseMultilist2(listBase);

				double tde_va = list[0] * s->p->par_eSPA1A + list[1] * s->p->par_eSPA1B + list[2] * s->p->par_eSPA1V + list[3] * s->p->par_eSPB1V + list[4] * s->p->par_eSPV1V +
						list[5] * s->p->par_eSPA2A + list[6] * s->p->par_eSPA2B + list[7] * s->p->par_eSPA2V + list[8] * s->p->par_eSPB2V + list[9] * s->p->par_eSPV2V;

				w(inv0) = exp(-(tde_va / 2 + s->p->par_eSPA) / (s->p->par_beta * s->p->par_temp));

		}
		else {
				w(inv0) = 0.0;
		}
}

void MultiFccJump::UpdateNeighbor()
{
	if (s->Site(fnf(path), x(path), y(path), z(path)) == 0)
	{	
		UpdateMutipleFN1BeforeNeighbor(s, x(path), y(path), z(path), fnf(path),0);
	
		UpdateMutipleN2BeforeNeighbor(s, x(path), y(path), z(path), fnf(path), 0);
		
		UpdateMutipleFN1AfterNeighbor(s, s->vacx, s->vacy, s->vacz,s->vacnm, 0);
	
		UpdateMutipleN2AfterNeighbor(s, s->vacx, s->vacy, s->vacz, s->vacnm, 0);
	
	}
	else if (s->Site(fnf(path), x(path), y(path), z(path)) == 1)
	{
		UpdateMutipleFN1BeforeNeighbor(s, x(path), y(path), z(path), fnf(path), 1);

		UpdateMutipleN2BeforeNeighbor(s, x(path), y(path), z(path), fnf(path), 1);

		UpdateMutipleFN1AfterNeighbor(s, s->vacx, s->vacy, s->vacz, s->vacnm, 1);

		UpdateMutipleN2AfterNeighbor(s, s->vacx, s->vacy, s->vacz, s->vacnm, 1);
		
	}
}

void MultiFccJump::runSimulationLoop()
{
	MultiFccJump::runFSimulationLoop(fnf, s);
}


//MultiSizeBcc
void MultiSizeBccJump::CalculatedEnergy(int inv, int nm, int inv0)
{
		if (s->Site(nf, x(inv), y(inv), z(inv)) == 1)
		{
			int a = s->NN1(1, nf, x(inv), y(inv), z(inv)); int b = s->NN1(0, nf, x(inv), y(inv), z(inv)); int c = s->NN1(2, nf, x(inv), y(inv), z(inv));
			int h = s->NN1(3, nf, x(inv), y(inv), z(inv)); int i = s->NN1(4, nf, x(inv), y(inv), z(inv));
			int d = s->NN1(1, s->vacnm, s->vacx, s->vacy, s->vacz); int e = s->NN1(0, s->vacnm, s->vacx, s->vacy, s->vacz); int f = s->NN1(2, s->vacnm, s->vacx, s->vacy, s->vacz);
			int j = s->NN1(3, s->vacnm, s->vacx, s->vacy, s->vacz); int k = s->NN1(4, s->vacnm, s->vacx, s->vacy, s->vacz);

			int l = d - 1 - a; int m = e - b; int n = f + a + 2 - c - d; int o = b - e; int p = c - f - 1;
			int q = h - j; int r = i - k; int t = j - h; int u = k - i;

			int aa = s->NN2(1, nf, x(inv), y(inv), z(inv)); int bb = s->NN2(0, nf, x(inv), y(inv), z(inv)); int cc = s->NN2(2, nf, x(inv), y(inv), z(inv));
			int hh = s->NN2(3, nf, x(inv), y(inv), z(inv)); int ii = s->NN2(4, nf, x(inv), y(inv), z(inv));
			int dd = s->NN2(1, s->vacnm, s->vacx, s->vacy, s->vacz); int ee = s->NN2(0, s->vacnm, s->vacx, s->vacy, s->vacz); int ff = s->NN2(2, s->vacnm, s->vacx, s->vacy, s->vacz);
			int jj = s->NN2(3, s->vacnm, s->vacx, s->vacy, s->vacz); int kk = s->NN2(4, s->vacnm, s->vacx, s->vacy, s->vacz);

			int ll = dd - aa; int mm = ee - bb; int nn = ff + aa - cc - dd; int oo = bb - ee; int pp = cc - ff;
			int qq = hh - jj; int rr = ii - kk; int tt = jj - hh; int uu = kk - ii;


			double tde_vb = l * s->p->par_eSPB1B + m * s->p->par_eSPA1B + n * s->p->par_eSPB1V + o * s->p->par_eSPA1V + p * s->p->par_eSPV1V +
				q * s->p->par_eSPC1V + r * s->p->par_eSPD1V + t * s->p->par_eSPA1C + u * s->p->par_eSPA1D +
				ll * s->p->par_eSPB2B + mm * s->p->par_eSPA2B + nn * s->p->par_eSPB2V + oo * s->p->par_eSPA2V + pp * s->p->par_eSPV2V +
				qq * s->p->par_eSPC2V + rr * s->p->par_eSPD2V + tt * s->p->par_eSPA2C + uu * s->p->par_eSPA2D;

			w(inv) = exp(-(tde_vb / 2 + s->p->par_eSPB) / (s->p->par_beta * s->p->par_temp));


		}
		else if (s->Site(nf, x(inv), y(inv), z(inv)) == 0)
		{
			int a = s->NN1(1, nf, x(inv), y(inv), z(inv)); int b = s->NN1(0, nf, x(inv), y(inv), z(inv)); int c = s->NN1(2, nf, x(inv), y(inv), z(inv));
			int h = s->NN1(3, nf, x(inv), y(inv), z(inv)); int i = s->NN1(4, nf, x(inv), y(inv), z(inv));
			int d = s->NN1(1, s->vacnm, s->vacx, s->vacy, s->vacz); int e = s->NN1(0, s->vacnm, s->vacx, s->vacy, s->vacz); int f = s->NN1(2, s->vacnm, s->vacx, s->vacy, s->vacz);
			int j = s->NN1(3, s->vacnm, s->vacx, s->vacy, s->vacz); int k = s->NN1(4, s->vacnm, s->vacx, s->vacy, s->vacz);

			int l = e - 1 - b; int m = d - a; int n = f + b + 2 - e - c; int o = a - d; int p = c - f - 1;
			int q = h - j; int r = i - k; int t = j - h; int u = k - i;

			int aa = s->NN2(1, nf, x(inv), y(inv), z(inv)); int bb = s->NN2(0, nf, x(inv), y(inv), z(inv)); int cc = s->NN2(2, nf, x(inv), y(inv), z(inv));
			int hh = s->NN2(3, nf, x(inv), y(inv), z(inv)); int ii = s->NN2(4, nf, x(inv), y(inv), z(inv));
			int dd = s->NN2(1, s->vacnm, s->vacx, s->vacy, s->vacz); int ee = s->NN2(0, s->vacnm, s->vacx, s->vacy, s->vacz); int ff = s->NN2(2, s->vacnm, s->vacx, s->vacy, s->vacz);
			int jj = s->NN2(3, s->vacnm, s->vacx, s->vacy, s->vacz); int kk = s->NN2(4, s->vacnm, s->vacx, s->vacy, s->vacz);

			int ll = ee - bb; int mm = dd - aa; int nn = ff + bb - cc - ee; int oo = aa - dd; int pp = cc - ff;
			int qq = hh - jj; int rr = ii - kk; int tt = jj - hh; int uu = kk - ii;


			double tde_va = l * s->p->par_eSPA1A + m * s->p->par_eSPA1B + n * s->p->par_eSPA1V + o * s->p->par_eSPB1V + p * s->p->par_eSPV1V +
				q * s->p->par_eSPC1V + r * s->p->par_eSPD1V + t * s->p->par_eSPA1C + u * s->p->par_eSPA1D +
				ll * s->p->par_eSPA2A + mm * s->p->par_eSPA2B + nn * s->p->par_eSPA2V + oo * s->p->par_eSPB2V + pp * s->p->par_eSPV2V +
				qq * s->p->par_eSPC2V + rr * s->p->par_eSPD2V + tt * s->p->par_eSPA2C + uu * s->p->par_eSPA2D;

			w(inv) = exp(-(tde_va / 2 + s->p->par_eSPA) / (s->p->par_beta * s->p->par_temp));

		}
		else if (s->Site(nf, x(inv), y(inv), z(inv)) == 3)
		{
			int a = s->NN1(1, nf, x(inv), y(inv), z(inv)); int b = s->NN1(0, nf, x(inv), y(inv), z(inv)); int c = s->NN1(2, nf, x(inv), y(inv), z(inv));
			int h = s->NN1(3, nf, x(inv), y(inv), z(inv)); int i = s->NN1(4, nf, x(inv), y(inv), z(inv));
			int d = s->NN1(1, s->vacnm, s->vacx, s->vacy, s->vacz); int e = s->NN1(0, s->vacnm, s->vacx, s->vacy, s->vacz); int f = s->NN1(2, s->vacnm, s->vacx, s->vacy, s->vacz);
			int j = s->NN1(3, s->vacnm, s->vacx, s->vacy, s->vacz); int k = s->NN1(4, s->vacnm, s->vacx, s->vacy, s->vacz);

			int l = j - 1 - h; int m = e - b; int n = f + h + 2 - c - j; int o = b - e; int p = c - f - 1;
			int q = a - d; int r = i - k; int t = d - a; int u = k - i;

			int aa = s->NN2(1, nf, x(inv), y(inv), z(inv)); int bb = s->NN2(0, nf, x(inv), y(inv), z(inv)); int cc = s->NN2(2, nf, x(inv), y(inv), z(inv));
			int hh = s->NN2(3, nf, x(inv), y(inv), z(inv)); int ii = s->NN2(4, nf, x(inv), y(inv), z(inv));
			int dd = s->NN2(1, s->vacnm, s->vacx, s->vacy, s->vacz); int ee = s->NN2(0, s->vacnm, s->vacx, s->vacy, s->vacz); int ff = s->NN2(2, s->vacnm, s->vacx, s->vacy, s->vacz);
			int jj = s->NN2(3, s->vacnm, s->vacx, s->vacy, s->vacz); int kk = s->NN2(4, s->vacnm, s->vacx, s->vacy, s->vacz);

			int ll = jj - hh; int mm = ee - bb; int nn = ff + hh - cc - jj; int oo = bb - ee; int pp = cc - ff;
			int qq = aa - dd; int rr = ii - kk; int tt = dd - aa; int uu = kk - ii;


			double tde_va = l * s->p->par_eSPC1C + m * s->p->par_eSPA1C + n * s->p->par_eSPC1V + o * s->p->par_eSPA1V + p * s->p->par_eSPV1V +
				q * s->p->par_eSPB1V + r * s->p->par_eSPD1V + t * s->p->par_eSPB1C + u * s->p->par_eSPC1D +
				ll * s->p->par_eSPC2C + mm * s->p->par_eSPA2C + nn * s->p->par_eSPC2V + oo * s->p->par_eSPA2V + pp * s->p->par_eSPV2V +
				qq * s->p->par_eSPB2V + rr * s->p->par_eSPD2V + tt * s->p->par_eSPB2C + uu * s->p->par_eSPC2D;

			w(inv) = exp(-(tde_va / 2 + s->p->par_eSPC) / (s->p->par_beta * s->p->par_temp));

		}
		else if (s->Site(nf, x(inv), y(inv), z(inv)) == 4)
		{
			int a = s->NN1(1, nf, x(inv), y(inv), z(inv)); int b = s->NN1(0, nf, x(inv), y(inv), z(inv)); int c = s->NN1(2, nf, x(inv), y(inv), z(inv));
			int h = s->NN1(3, nf, x(inv), y(inv), z(inv)); int i = s->NN1(4, nf, x(inv), y(inv), z(inv));
			int d = s->NN1(1, s->vacnm, s->vacx, s->vacy, s->vacz); int e = s->NN1(0, s->vacnm, s->vacx, s->vacy, s->vacz); int f = s->NN1(2, s->vacnm, s->vacx, s->vacy, s->vacz);
			int j = s->NN1(3, s->vacnm, s->vacx, s->vacy, s->vacz); int k = s->NN1(4, s->vacnm, s->vacx, s->vacy, s->vacz);

			int l = k - 1 - i; int m = e - b; int n = f + i + 2 - c - k; int o = b - e; int p = c - f - 1;
			int q = a - d; int r = h - j; int t = d - a; int u = j - h;

			int aa = s->NN2(1, nf, x(inv), y(inv), z(inv)); int bb = s->NN2(0, nf, x(inv), y(inv), z(inv)); int cc = s->NN2(2, nf, x(inv), y(inv), z(inv));
			int hh = s->NN2(3, nf, x(inv), y(inv), z(inv)); int ii = s->NN2(4, nf, x(inv), y(inv), z(inv));
			int dd = s->NN2(1, s->vacnm, s->vacx, s->vacy, s->vacz); int ee = s->NN2(0, s->vacnm, s->vacx, s->vacy, s->vacz); int ff = s->NN2(2, s->vacnm, s->vacx, s->vacy, s->vacz);
			int jj = s->NN2(3, s->vacnm, s->vacx, s->vacy, s->vacz); int kk = s->NN2(4, s->vacnm, s->vacx, s->vacy, s->vacz);

			int ll = kk - ii; int mm = ee - bb; int nn = ff + ii - cc - kk; int oo = bb - ee; int pp = cc - ff;
			int qq = aa - dd; int rr = hh - jj; int tt = dd - aa; int uu = jj - hh;


			double tde_va = l * s->p->par_eSPD1D + m * s->p->par_eSPA1D + n * s->p->par_eSPA1V + o * s->p->par_eSPD1V + p * s->p->par_eSPV1V +
				q * s->p->par_eSPB1V + r * s->p->par_eSPC1V + t * s->p->par_eSPB1D + u * s->p->par_eSPC1D +
				ll * s->p->par_eSPD2D + mm * s->p->par_eSPA2D + nn * s->p->par_eSPA2V + oo * s->p->par_eSPD2V + pp * s->p->par_eSPV2V +
				qq * s->p->par_eSPB2V + rr * s->p->par_eSPC2V + tt * s->p->par_eSPB2D + uu * s->p->par_eSPC2D;

			w(inv) = exp(-(tde_va / 2 + s->p->par_eSPD) / (s->p->par_beta * s->p->par_temp));

		}
		else {
			w(inv) = 0.0;
		}
	
}

void MultiSizeBccJump::UpdateNeighbor()
{
	if (s->Site(nf, x(path), y(path), z(path)) == 0) {

		UpdateMutipleN1BeforeNeighbor(s, x(path), y(path), z(path), nf, s->vacnm, 0);
	
		UpdateMutipleN2BeforeNeighbor(s, x(path), y(path), z(path), nf, 0);
	
		UpdateMutipleN1AfterNeighbor(s, s->vacx, s->vacy, s->vacz, nf, s->vacnm, 0);

		UpdateMutipleN2AfterNeighbor(s, s->vacx, s->vacy, s->vacz, s->vacnm, 0);
		
	}
	else if (s->Site(nf, x(path), y(path), z(path)) == 1) {

		UpdateMutipleN1BeforeNeighbor(s, x(path), y(path), z(path), nf, s->vacnm, 1);
	
		UpdateMutipleN2BeforeNeighbor(s, x(path), y(path), z(path), nf, 1);
	
		UpdateMutipleN1AfterNeighbor(s, s->vacx, s->vacy, s->vacz, nf, s->vacnm, 1);
	
		UpdateMutipleN2AfterNeighbor(s, s->vacx, s->vacy, s->vacz, s->vacnm, 1);
	
	}
	else if (s->Site(nf, x(path), y(path), z(path)) == 3) {

		UpdateMutipleN1BeforeNeighbor(s, x(path), y(path), z(path), nf, s->vacnm, 3);

		UpdateMutipleN2BeforeNeighbor(s, x(path), y(path), z(path), nf, 3);

		UpdateMutipleN1AfterNeighbor(s, s->vacx, s->vacy, s->vacz, nf, s->vacnm, 3);

		UpdateMutipleN2AfterNeighbor(s, s->vacx, s->vacy, s->vacz, s->vacnm, 3);

	}
	else if (s->Site(nf, x(path), y(path), z(path)) == 4) {

		UpdateMutipleN1BeforeNeighbor(s, x(path), y(path), z(path), nf, s->vacnm, 4);

		UpdateMutipleN2BeforeNeighbor(s, x(path), y(path), z(path), nf, 4);

		UpdateMutipleN1AfterNeighbor(s, s->vacx, s->vacy, s->vacz, nf, s->vacnm, 4);

		UpdateMutipleN2AfterNeighbor(s, s->vacx, s->vacy, s->vacz, s->vacnm, 4);

	}
}

void MultiSizeBccJump::runSimulationLoop()
{
	MultiSizeBccJump::runBSimulationLoop(nf, s);
}


//MultiSizeFcc
void MultiSizeFccJump::CalculatedEnergy(int inv, int nf, int inv0)
{
	if (s->Site(nf, x(inv), y(inv), z(inv)) == 1)
	{
		int a = s->NN1(1, nf, x(inv), y(inv), z(inv)); int b = s->NN1(0, nf, x(inv), y(inv), z(inv)); int c = s->NN1(2, nf, x(inv), y(inv), z(inv));
		int h = s->NN1(3, nf, x(inv), y(inv), z(inv)); int i = s->NN1(4, nf, x(inv), y(inv), z(inv));
		int d = s->NN1(1, s->vacnm, s->vacx, s->vacy, s->vacz); int e = s->NN1(0, s->vacnm, s->vacx, s->vacy, s->vacz); int f = s->NN1(2, s->vacnm, s->vacx, s->vacy, s->vacz);
		int j = s->NN1(3, s->vacnm, s->vacx, s->vacy, s->vacz); int k = s->NN1(4, s->vacnm, s->vacx, s->vacy, s->vacz);

		int l = d - 1 - a; int m = e - b; int n = f + a + 2 - c - d; int o = b - e; int p = c - f - 1;
		int q = h - j; int r = i - k; int t = j - h; int u = k - i;

		int aa = s->NN2(1, nf, x(inv), y(inv), z(inv)); int bb = s->NN2(0, nf, x(inv), y(inv), z(inv)); int cc = s->NN2(2, nf, x(inv), y(inv), z(inv));
		int hh = s->NN2(3, nf, x(inv), y(inv), z(inv)); int ii = s->NN2(4, nf, x(inv), y(inv), z(inv));
		int dd = s->NN2(1, s->vacnm, s->vacx, s->vacy, s->vacz); int ee = s->NN2(0, s->vacnm, s->vacx, s->vacy, s->vacz); int ff = s->NN2(2, s->vacnm, s->vacx, s->vacy, s->vacz);
		int jj = s->NN2(3, s->vacnm, s->vacx, s->vacy, s->vacz); int kk = s->NN2(4, s->vacnm, s->vacx, s->vacy, s->vacz);

		int ll = dd - aa; int mm = ee - bb; int nn = ff + aa - cc - dd; int oo = bb - ee; int pp = cc - ff;
		int qq = hh - jj; int rr = ii - kk; int tt = jj - hh; int uu = kk - ii;


		double tde_vb = l * s->p->par_eSPB1B + m * s->p->par_eSPA1B + n * s->p->par_eSPB1V + o * s->p->par_eSPA1V + p * s->p->par_eSPV1V +
			q * s->p->par_eSPC1V + r * s->p->par_eSPD1V + t * s->p->par_eSPA1C + u * s->p->par_eSPA1D +
			ll * s->p->par_eSPB2B + mm * s->p->par_eSPA2B + nn * s->p->par_eSPB2V + oo * s->p->par_eSPA2V + pp * s->p->par_eSPV2V +
			qq * s->p->par_eSPC2V + rr * s->p->par_eSPD2V + tt * s->p->par_eSPA2C + uu * s->p->par_eSPA2D;

		w(inv0) = exp(-(tde_vb / 2 + s->p->par_eSPB) / (s->p->par_beta * s->p->par_temp));


	}
	else if (s->Site(nf, x(inv), y(inv), z(inv)) == 0)
	{
		int a = s->NN1(1, nf, x(inv), y(inv), z(inv)); int b = s->NN1(0, nf, x(inv), y(inv), z(inv)); int c = s->NN1(2, nf, x(inv), y(inv), z(inv));
		int h = s->NN1(3, nf, x(inv), y(inv), z(inv)); int i = s->NN1(4, nf, x(inv), y(inv), z(inv));
		int d = s->NN1(1, s->vacnm, s->vacx, s->vacy, s->vacz); int e = s->NN1(0, s->vacnm, s->vacx, s->vacy, s->vacz); int f = s->NN1(2, s->vacnm, s->vacx, s->vacy, s->vacz);
		int j = s->NN1(3, s->vacnm, s->vacx, s->vacy, s->vacz); int k = s->NN1(4, s->vacnm, s->vacx, s->vacy, s->vacz);

		int l = e - 1 - b; int m = d - a; int n = f + b + 2 - e - c; int o = a - d; int p = c - f - 1;
		int q = h - j; int r = i - k; int t = j - h; int u = k - i;

		int aa = s->NN2(1, nf, x(inv), y(inv), z(inv)); int bb = s->NN2(0, nf, x(inv), y(inv), z(inv)); int cc = s->NN2(2, nf, x(inv), y(inv), z(inv));
		int hh = s->NN2(3, nf, x(inv), y(inv), z(inv)); int ii = s->NN2(4, nf, x(inv), y(inv), z(inv));
		int dd = s->NN2(1, s->vacnm, s->vacx, s->vacy, s->vacz); int ee = s->NN2(0, s->vacnm, s->vacx, s->vacy, s->vacz); int ff = s->NN2(2, s->vacnm, s->vacx, s->vacy, s->vacz);
		int jj = s->NN2(3, s->vacnm, s->vacx, s->vacy, s->vacz); int kk = s->NN2(4, s->vacnm, s->vacx, s->vacy, s->vacz);

		int ll = ee - bb; int mm = dd - aa; int nn = ff + bb - cc - ee; int oo = aa - dd; int pp = cc - ff;
		int qq = hh - jj; int rr = ii - kk; int tt = jj - hh; int uu = kk - ii;


		double tde_va = l * s->p->par_eSPA1A + m * s->p->par_eSPA1B + n * s->p->par_eSPA1V + o * s->p->par_eSPB1V + p * s->p->par_eSPV1V +
			q * s->p->par_eSPC1V + r * s->p->par_eSPD1V + t * s->p->par_eSPA1C + u * s->p->par_eSPA1D +
			ll * s->p->par_eSPA2A + mm * s->p->par_eSPA2B + nn * s->p->par_eSPA2V + oo * s->p->par_eSPB2V + pp * s->p->par_eSPV2V +
			qq * s->p->par_eSPC2V + rr * s->p->par_eSPD2V + tt * s->p->par_eSPA2C + uu * s->p->par_eSPA2D;

		w(inv0) = exp(-(tde_va / 2 + s->p->par_eSPA) / (s->p->par_beta * s->p->par_temp));

	}
	else if (s->Site(nf, x(inv), y(inv), z(inv)) == 3)
	{
		int a = s->NN1(1, nf, x(inv), y(inv), z(inv)); int b = s->NN1(0, nf, x(inv), y(inv), z(inv)); int c = s->NN1(2, nf, x(inv), y(inv), z(inv));
		int h = s->NN1(3, nf, x(inv), y(inv), z(inv)); int i = s->NN1(4, nf, x(inv), y(inv), z(inv));
		int d = s->NN1(1, s->vacnm, s->vacx, s->vacy, s->vacz); int e = s->NN1(0, s->vacnm, s->vacx, s->vacy, s->vacz); int f = s->NN1(2, s->vacnm, s->vacx, s->vacy, s->vacz);
		int j = s->NN1(3, s->vacnm, s->vacx, s->vacy, s->vacz); int k = s->NN1(4, s->vacnm, s->vacx, s->vacy, s->vacz);

		int l = j - 1 - h; int m = e - b; int n = f + h + 2 - c - j; int o = b - e; int p = c - f - 1;
		int q = a - d; int r = i - k; int t = d - a; int u = k - i;

		int aa = s->NN2(1, nf, x(inv), y(inv), z(inv)); int bb = s->NN2(0, nf, x(inv), y(inv), z(inv)); int cc = s->NN2(2, nf, x(inv), y(inv), z(inv));
		int hh = s->NN2(3, nf, x(inv), y(inv), z(inv)); int ii = s->NN2(4, nf, x(inv), y(inv), z(inv));
		int dd = s->NN2(1, s->vacnm, s->vacx, s->vacy, s->vacz); int ee = s->NN2(0, s->vacnm, s->vacx, s->vacy, s->vacz); int ff = s->NN2(2, s->vacnm, s->vacx, s->vacy, s->vacz);
		int jj = s->NN2(3, s->vacnm, s->vacx, s->vacy, s->vacz); int kk = s->NN2(4, s->vacnm, s->vacx, s->vacy, s->vacz);

		int ll = jj - hh; int mm = ee - bb; int nn = ff + hh - cc - jj; int oo = bb - ee; int pp = cc - ff;
		int qq = aa - dd; int rr = ii - kk; int tt = dd - aa; int uu = kk - ii;


		double tde_va = l * s->p->par_eSPC1C + m * s->p->par_eSPA1C + n * s->p->par_eSPC1V + o * s->p->par_eSPA1V + p * s->p->par_eSPV1V +
			q * s->p->par_eSPB1V + r * s->p->par_eSPD1V + t * s->p->par_eSPB1C + u * s->p->par_eSPC1D +
			ll * s->p->par_eSPC2C + mm * s->p->par_eSPA2C + nn * s->p->par_eSPC2V + oo * s->p->par_eSPA2V + pp * s->p->par_eSPV2V +
			qq * s->p->par_eSPB2V + rr * s->p->par_eSPD2V + tt * s->p->par_eSPB2C + uu * s->p->par_eSPC2D;

		w(inv0) = exp(-(tde_va / 2 + s->p->par_eSPC) / (s->p->par_beta * s->p->par_temp));

	}
	else if (s->Site(nf, x(inv), y(inv), z(inv)) == 4)
	{
		int a = s->NN1(1, nf, x(inv), y(inv), z(inv)); int b = s->NN1(0, nf, x(inv), y(inv), z(inv)); int c = s->NN1(2, nf, x(inv), y(inv), z(inv));
		int h = s->NN1(3, nf, x(inv), y(inv), z(inv)); int i = s->NN1(4, nf, x(inv), y(inv), z(inv));
		int d = s->NN1(1, s->vacnm, s->vacx, s->vacy, s->vacz); int e = s->NN1(0, s->vacnm, s->vacx, s->vacy, s->vacz); int f = s->NN1(2, s->vacnm, s->vacx, s->vacy, s->vacz);
		int j = s->NN1(3, s->vacnm, s->vacx, s->vacy, s->vacz); int k = s->NN1(4, s->vacnm, s->vacx, s->vacy, s->vacz);

		int l = k - 1 - i; int m = e - b; int n = f + i + 2 - c - k; int o = b - e; int p = c - f - 1;
		int q = a - d; int r = h - j; int t = d - a; int u = j - h;

		int aa = s->NN2(1, nf, x(inv), y(inv), z(inv)); int bb = s->NN2(0, nf, x(inv), y(inv), z(inv)); int cc = s->NN2(2, nf, x(inv), y(inv), z(inv));
		int hh = s->NN2(3, nf, x(inv), y(inv), z(inv)); int ii = s->NN2(4, nf, x(inv), y(inv), z(inv));
		int dd = s->NN2(1, s->vacnm, s->vacx, s->vacy, s->vacz); int ee = s->NN2(0, s->vacnm, s->vacx, s->vacy, s->vacz); int ff = s->NN2(2, s->vacnm, s->vacx, s->vacy, s->vacz);
		int jj = s->NN2(3, s->vacnm, s->vacx, s->vacy, s->vacz); int kk = s->NN2(4, s->vacnm, s->vacx, s->vacy, s->vacz);

		int ll = kk - ii; int mm = ee - bb; int nn = ff + ii - cc - kk; int oo = bb - ee; int pp = cc - ff;
		int qq = aa - dd; int rr = hh - jj; int tt = dd - aa; int uu = jj - hh;


		double tde_va = l * s->p->par_eSPD1D + m * s->p->par_eSPA1D + n * s->p->par_eSPA1V + o * s->p->par_eSPD1V + p * s->p->par_eSPV1V +
			q * s->p->par_eSPB1V + r * s->p->par_eSPC1V + t * s->p->par_eSPB1D + u * s->p->par_eSPC1D +
			ll * s->p->par_eSPD2D + mm * s->p->par_eSPA2D + nn * s->p->par_eSPA2V + oo * s->p->par_eSPD2V + pp * s->p->par_eSPV2V +
			qq * s->p->par_eSPB2V + rr * s->p->par_eSPC2V + tt * s->p->par_eSPB2D + uu * s->p->par_eSPC2D;

		w(inv0) = exp(-(tde_va / 2 + s->p->par_eSPD) / (s->p->par_beta * s->p->par_temp));

	}
	else {
		w(inv0) = 0.0;
	}
			
				
}

void MultiSizeFccJump::UpdateNeighbor()
{
	if (s->Site(fnf(path), x(path), y(path), z(path)) == 0)
	{
		UpdateMutipleFN1BeforeNeighbor(s, x(path), y(path), z(path), fnf(path), 0);
		
		UpdateMutipleN2BeforeNeighbor(s, x(path), y(path), z(path), fnf(path), 0);
		
		UpdateMutipleFN1AfterNeighbor(s, s->vacx, s->vacy, s->vacz, s->vacnm, 0);
	
		UpdateMutipleN2AfterNeighbor(s, s->vacx, s->vacy, s->vacz, s->vacnm, 0);

	}
	else if (s->Site(fnf(path), x(path), y(path), z(path)) == 1)
	{
		UpdateMutipleFN1BeforeNeighbor(s, x(path), y(path), z(path), fnf(path), 1);
		
		UpdateMutipleN2BeforeNeighbor(s, x(path), y(path), z(path), fnf(path), 1);
	
		UpdateMutipleFN1AfterNeighbor(s, s->vacx, s->vacy, s->vacz, s->vacnm, 1);
		
		UpdateMutipleN2AfterNeighbor(s, s->vacx, s->vacy, s->vacz, s->vacnm, 1);

	}
	else if (s->Site(fnf(path), x(path), y(path), z(path)) == 3)
	{
		UpdateMutipleFN1BeforeNeighbor(s, x(path), y(path), z(path), fnf(path), 3);
	
		UpdateMutipleN2BeforeNeighbor(s, x(path), y(path), z(path), fnf(path), 3);
		
		UpdateMutipleFN1AfterNeighbor(s, s->vacx, s->vacy, s->vacz, s->vacnm, 3);

		UpdateMutipleN2AfterNeighbor(s, s->vacx, s->vacy, s->vacz, s->vacnm, 3);

	}
	else if (s->Site(fnf(path), x(path), y(path), z(path)) == 4)
	{
		UpdateMutipleFN1BeforeNeighbor(s, x(path), y(path), z(path), fnf(path), 4);

		UpdateMutipleN2BeforeNeighbor(s, x(path), y(path), z(path), fnf(path), 4);
		
		UpdateMutipleFN1AfterNeighbor(s, s->vacx, s->vacy, s->vacz, s->vacnm, 4);
		
		UpdateMutipleN2AfterNeighbor(s, s->vacx, s->vacy, s->vacz, s->vacnm, 4);
	
	}
}

void MultiSizeFccJump::runSimulationLoop()
{
	MultiSizeFccJump::runFSimulationLoop(fnf, s);
}



/**
 * @brief 执行蒙特卡洛模拟 (MC) 的一个步骤，根据规则交换空位的位置。
 *
 * @param i 当前处理区域的 x 索引
 * @param j 当前处理区域的 y 索引
 * @param k 当前处理区域的 z 索引
 * @param idx x 方向的索引步长
 * @param idy y 方向的索引步长
 * @param idz z 方向的索引步长
 *
 * 该函数执行蒙特卡洛模拟 (MC) 的一个步骤。它随机选择当前处理区域内的一个空位进行操作，
 * 然后根据一定的规则交换空位的位置。如果空位越出当前处理区域的边界，它将被加入到越界区域的 vaclists 中，
 * 并从当前区域的 vaclists 中删除。在一系列 MC 步骤中，该函数不断选择并处理不同的空位，
 * 直到满足某个退出条件为止。
 */
void ParaBccJump::Sector(int i, int j, int k, int idx, int idy, int idz) {

		// 初始化随机数生成器，通常使用当前时间作为种子
		srand(static_cast<unsigned int>(time(nullptr)));

		// 随机选择s->Vaclists(i, j, k)中的一个进行mc
		int randomNumber1 = rand() % (s->Vaclists(i, j, k).size());
		int randomNumber2 = -1; //使用-1表示 randomNumber2 没有被赋值
		
		auto* vac = &s->Vaclists(i, j, k)[randomNumber1];

		int nf = 0; int path1 = 0;//局部变量优先级更高
		blitz::Array<int, 1> x1, y1, z1; blitz::Array<double, 1> w1;
		x1.resize(8); y1.resize(8); z1.resize(8); w1.resize(8);
		x1 = 0; y1 = 0; z1 = 0; w1 = 0;

		for (int step = 0; step < 1000; step++) {


			nf = (vac->nm == 0) ? 1 : 0;

			CalculatedEnergys(vac->x, vac->y, vac->z, vac->nm, nf, x1, y1, z1, w1);

			Choice(path1, w1);

			if (s->Site(nf, x1(path1), y1(path1), z1(path1)) == 2) {
				// 如果交换后的位置为2，表示无法交换，重新选择一个空位
				randomNumber2 = rand() % (s->Vaclists(i, j, k).size());
				vac = &(s->Vaclists(i, j, k)[randomNumber2]);
				continue;
			}

			UpdateNeighbors(vac->x, vac->y, vac->z, vac->nm, nf, x1, y1, z1, w1, path1);
			ExchangeSites(vac->x, vac->y, vac->z, vac->nm, nf, x1, y1, z1, w1, path1);

			//如果空位跳出区域
			if (Border(vac->x, vac->y, vac->z, vac->nm, i, j, k, idx, idy, idz)) {

				int nexti = vac->x / idx; int nextj = vac->y / idy; int nextk = vac->z / idz;
				// 将空位加入到越界区域的vaclists
				s->Vaclists(nexti, nextj, nextk).emplace_back(*vac);

				// 从当前区域的vaclists中删除空位
				if (randomNumber2 != -1) {
					s->Vaclists(i, j, k).erase(s->Vaclists(i, j, k).begin() + randomNumber2);
				}
				else {
					s->Vaclists(i, j, k).erase(s->Vaclists(i, j, k).begin() + randomNumber1);
				}

				// 选择一个新的vac
				if (s->Vaclists(i, j, k).empty()) {
					break;
				}
				else {

					randomNumber2 = rand() % (s->Vaclists(i, j, k).size());
					vac = &(s->Vaclists(i, j, k)[randomNumber2]);
				}

			}
			else {
				// 如果空位没有跳出区域，选择一个新的空位进行下一步操作
				randomNumber2 = rand() % (s->Vaclists(i, j, k).size());
				vac = &(s->Vaclists(i, j, k)[randomNumber2]);
			}
		}
}


/**
 * @brief 并行分割任务以在多个线程上运行 Monte Carlo 模拟。
 *
 * 使用 OpenMP 进行并行化处理，根据并行块的索引和循环遍历，以便在不相邻的区域上运行 Monte Carlo 模拟。
 * 每个线程负责处理一个不相邻区域的任务。
 */
void ParaBccJump::Divide()
{
	//for (int n = 0; n < 8; n++) 是单线程执行，循环迭代不会并行化。
	//内部的 #pragma omp parallel for 是多线程执行，
	//用于并行化内部的三个嵌套循环，以充分利用多核处理器的性能。

	omp_set_num_threads(s->p->numThreads);


	int idx = s->p->nx / s->p->parallelx;
	int idy = s->p->ny / s->p->parallely;
	int idz = s->p->nz / s->p->parallelz;


	int nums[8][3] = { {0,0,0},{1,0,0},{1,1,0},{1,1,1},{1,0,1},{0,0,1},{0,1,1},{0,1,0} };
	int a, b, c;

	for (int n = 0; n < 8; n++) {

		a = nums[n][0]; b = nums[n][1]; c = nums[n][2];

		#pragma omp parallel for
		for (int i = a; i < s->p->parallelx; i += 2) {
			for (int j = b; j < s->p->parallely; j += 2) {
				for (int k = c; k < s->p->parallelz; k += 2) {

					// 遍历当前不相邻的区域
					auto& vac_list = s->Vaclists(i, j, k);
					if (vac_list.empty()) {

						continue;
					}
					else {

						//s->Vaclists(i, j, k)内有空位，那么就运行mc
						Sector(i, j, k, idx, idy, idz);

					}
				}
			}
		}
		//#pragma omp barrier

	}

}

void ParaBccJump::CalculatedEnergys(int vacx, int vacy, int vacz, int vacnm, int nf,
									blitz::Array<int, 1> x, blitz::Array<int, 1> y,
									blitz::Array<int, 1> z, blitz::Array<double, 1> w)
{
	for (int inv = 0; inv < 8; inv++) {

		x(inv) = s->modlx[vacx + v1nbr_bcc[vacnm][0][inv]];
		y(inv) = s->modly[vacy + v1nbr_bcc[vacnm][1][inv]];
		z(inv) = s->modlz[vacz + v1nbr_bcc[vacnm][2][inv]];

		if (s->Site(nf, x(inv), y(inv), z(inv)) == 1)
		{
			
			int a = s->NN1(1, nf, x(inv), y(inv), z(inv)); int b = s->NN1(0, nf, x(inv), y(inv), z(inv)); int c = s->NN1(2, nf, x(inv), y(inv), z(inv));
			int h = s->NN1(3, nf, x(inv), y(inv), z(inv)); int i = s->NN1(4, nf, x(inv), y(inv), z(inv));
			int d = s->NN1(1, vacnm, vacx, vacy, vacz); int e = s->NN1(0, vacnm, vacx, vacy, vacz); int f = s->NN1(2, vacnm, vacx, vacy, vacz);
			int j = s->NN1(3, vacnm, vacx, vacy, vacz); int k = s->NN1(4, vacnm,vacx, vacy, vacz);

			int l = d - 1 - a; int m = e - b; int n = f + a + 2 - c - d; int o = b - e; int p = c - f - 1;
			int q = h - j; int r = i - k; int t = j - h; int u = k - i;

			int aa = s->NN2(1, nf, x(inv), y(inv), z(inv)); int bb = s->NN2(0, nf, x(inv), y(inv), z(inv)); int cc = s->NN2(2, nf, x(inv), y(inv), z(inv));
			int hh = s->NN2(3, nf, x(inv), y(inv), z(inv)); int ii = s->NN2(4, nf, x(inv), y(inv), z(inv));
			int dd = s->NN2(1,vacnm, vacx, vacy, vacz); int ee = s->NN2(0, vacnm, vacx, vacy, vacz); int ff = s->NN2(2, vacnm, vacx, vacy, vacz);
			int jj = s->NN2(3, vacnm, vacx, vacy, vacz); int kk = s->NN2(4, vacnm, vacx, vacy, vacz);

			int ll = dd - aa; int mm = ee - bb; int nn = ff + aa - cc - dd; int oo = bb - ee; int pp = cc - ff;
			int qq = hh - jj; int rr = ii - kk; int tt = jj - hh; int uu = kk - ii;


			double tde_vb = l * s->p->par_eSPB1B + m * s->p->par_eSPA1B + n * s->p->par_eSPB1V + o * s->p->par_eSPA1V + p * s->p->par_eSPV1V +
				q * s->p->par_eSPC1V + r * s->p->par_eSPD1V + t * s->p->par_eSPA1C + u * s->p->par_eSPA1D +
				ll * s->p->par_eSPB2B + mm * s->p->par_eSPA2B + nn * s->p->par_eSPB2V + oo * s->p->par_eSPA2V + pp * s->p->par_eSPV2V +
				qq * s->p->par_eSPC2V + rr * s->p->par_eSPD2V + tt * s->p->par_eSPA2C + uu * s->p->par_eSPA2D;

			w(inv) = exp(-(tde_vb / 2 + s->p->par_eSPB) / (s->p->par_beta * s->p->par_temp));


		}
		else if (s->Site(nf, x(inv), y(inv), z(inv)) == 0)
		{
			int a = s->NN1(1, nf, x(inv), y(inv), z(inv)); int b = s->NN1(0, nf, x(inv), y(inv), z(inv)); int c = s->NN1(2, nf, x(inv), y(inv), z(inv));
			int h = s->NN1(3, nf, x(inv), y(inv), z(inv)); int i = s->NN1(4, nf, x(inv), y(inv), z(inv));
			int d = s->NN1(1, vacnm, vacx, vacy, vacz); int e = s->NN1(0, vacnm, vacx, vacy, vacz); int f = s->NN1(2, vacnm, vacx, vacy, vacz);
			int j = s->NN1(3, vacnm, vacx, vacy, vacz); int k = s->NN1(4, vacnm, vacx, vacy, vacz);

			int l = e - 1 - b; int m = d - a; int n = f + b + 2 - e - c; int o = a - d; int p = c - f - 1;
			int q = h - j; int r = i - k; int t = j - h; int u = k - i;

			int aa = s->NN2(1, nf, x(inv), y(inv), z(inv)); int bb = s->NN2(0, nf, x(inv), y(inv), z(inv)); int cc = s->NN2(2, nf, x(inv), y(inv), z(inv));
			int hh = s->NN2(3, nf, x(inv), y(inv), z(inv)); int ii = s->NN2(4, nf, x(inv), y(inv), z(inv));
			int dd = s->NN2(1, vacnm, vacx, vacy, vacz); int ee = s->NN2(0, vacnm, vacx, vacy, vacz); int ff = s->NN2(2, vacnm, vacx, vacy, vacz);
			int jj = s->NN2(3, vacnm, vacx, vacy, vacz); int kk = s->NN2(4, vacnm, vacx, vacy, vacz);

			int ll = ee - bb; int mm = dd - aa; int nn = ff + bb - cc - ee; int oo = aa - dd; int pp = cc - ff;
			int qq = hh - jj; int rr = ii - kk; int tt = jj - hh; int uu = kk - ii;


			double tde_va = l * s->p->par_eSPA1A + m * s->p->par_eSPA1B + n * s->p->par_eSPA1V + o * s->p->par_eSPB1V + p * s->p->par_eSPV1V +
				q * s->p->par_eSPC1V + r * s->p->par_eSPD1V + t * s->p->par_eSPA1C + u * s->p->par_eSPA1D +
				ll * s->p->par_eSPA2A + mm * s->p->par_eSPA2B + nn * s->p->par_eSPA2V + oo * s->p->par_eSPB2V + pp * s->p->par_eSPV2V +
				qq * s->p->par_eSPC2V + rr * s->p->par_eSPD2V + tt * s->p->par_eSPA2C + uu * s->p->par_eSPA2D;

			w(inv) = exp(-(tde_va / 2 + s->p->par_eSPA) / (s->p->par_beta * s->p->par_temp));

		}
		else if (s->Site(nf, x(inv), y(inv), z(inv)) == 3)
		{
			int a = s->NN1(1, nf, x(inv), y(inv), z(inv)); int b = s->NN1(0, nf, x(inv), y(inv), z(inv)); int c = s->NN1(2, nf, x(inv), y(inv), z(inv));
			int h = s->NN1(3, nf, x(inv), y(inv), z(inv)); int i = s->NN1(4, nf, x(inv), y(inv), z(inv));
			int d = s->NN1(1, vacnm, vacx, vacy, vacz); int e = s->NN1(0, vacnm, vacx, vacy, vacz); int f = s->NN1(2, vacnm, vacx, vacy, vacz);
			int j = s->NN1(3, vacnm, vacx, vacy, vacz); int k = s->NN1(4, vacnm, vacx, vacy, vacz);

			int l = j - 1 - h; int m = e - b; int n = f + h + 2 - c - j; int o = b - e; int p = c - f - 1;
			int q = a - d; int r = i - k; int t = d - a; int u = k - i;

			int aa = s->NN2(1, nf, x(inv), y(inv), z(inv)); int bb = s->NN2(0, nf, x(inv), y(inv), z(inv)); int cc = s->NN2(2, nf, x(inv), y(inv), z(inv));
			int hh = s->NN2(3, nf, x(inv), y(inv), z(inv)); int ii = s->NN2(4, nf, x(inv), y(inv), z(inv));
			int dd = s->NN2(1, vacnm, vacx, vacy, vacz); int ee = s->NN2(0, vacnm, vacx, vacy, vacz); int ff = s->NN2(2, vacnm, vacx, vacy, vacz);
			int jj = s->NN2(3, vacnm, vacx, vacy, vacz); int kk = s->NN2(4, vacnm, vacx, vacy, vacz);

			int ll = jj - hh; int mm = ee - bb; int nn = ff + hh - cc - jj; int oo = bb - ee; int pp = cc - ff;
			int qq = aa - dd; int rr = ii - kk; int tt = dd - aa; int uu = kk - ii;


			double tde_va = l * s->p->par_eSPC1C + m * s->p->par_eSPA1C + n * s->p->par_eSPC1V + o * s->p->par_eSPA1V + p * s->p->par_eSPV1V +
				q * s->p->par_eSPB1V + r * s->p->par_eSPD1V + t * s->p->par_eSPB1C + u * s->p->par_eSPC1D +
				ll * s->p->par_eSPC2C + mm * s->p->par_eSPA2C + nn * s->p->par_eSPC2V + oo * s->p->par_eSPA2V + pp * s->p->par_eSPV2V +
				qq * s->p->par_eSPB2V + rr * s->p->par_eSPD2V + tt * s->p->par_eSPB2C + uu * s->p->par_eSPC2D;

			w(inv) = exp(-(tde_va / 2 + s->p->par_eSPC) / (s->p->par_beta * s->p->par_temp));

		}
		else if (s->Site(nf, x(inv), y(inv), z(inv)) == 4)
		{
			int a = s->NN1(1, nf, x(inv), y(inv), z(inv)); int b = s->NN1(0, nf, x(inv), y(inv), z(inv)); int c = s->NN1(2, nf, x(inv), y(inv), z(inv));
			int h = s->NN1(3, nf, x(inv), y(inv), z(inv)); int i = s->NN1(4, nf, x(inv), y(inv), z(inv));
			int d = s->NN1(1, vacnm, vacx, vacy, vacz); int e = s->NN1(0, vacnm, vacx, vacy, vacz); int f = s->NN1(2, vacnm, vacx, vacy, vacz);
			int j = s->NN1(3, vacnm, vacx, vacy, vacz); int k = s->NN1(4, vacnm, vacx, vacy, vacz);

			int l = k - 1 - i; int m = e - b; int n = f + i + 2 - c - k; int o = b - e; int p = c - f - 1;
			int q = a - d; int r = h - j; int t = d - a; int u = j - h;

			int aa = s->NN2(1, nf, x(inv), y(inv), z(inv)); int bb = s->NN2(0, nf, x(inv), y(inv), z(inv)); int cc = s->NN2(2, nf, x(inv), y(inv), z(inv));
			int hh = s->NN2(3, nf, x(inv), y(inv), z(inv)); int ii = s->NN2(4, nf, x(inv), y(inv), z(inv));
			int dd = s->NN2(1, vacnm, vacx, vacy, vacz); int ee = s->NN2(0, vacnm, vacx, vacy, vacz); int ff = s->NN2(2, vacnm, vacx, vacy, vacz);
			int jj = s->NN2(3, vacnm, vacx, vacy, vacz); int kk = s->NN2(4, vacnm, vacx, vacy, vacz);

			int ll = kk - ii; int mm = ee - bb; int nn = ff + ii - cc - kk; int oo = bb - ee; int pp = cc - ff;
			int qq = aa - dd; int rr = hh - jj; int tt = dd - aa; int uu = jj - hh;


			double tde_va = l * s->p->par_eSPD1D + m * s->p->par_eSPA1D + n * s->p->par_eSPA1V + o * s->p->par_eSPD1V + p * s->p->par_eSPV1V +
				q * s->p->par_eSPB1V + r * s->p->par_eSPC1V + t * s->p->par_eSPB1D + u * s->p->par_eSPC1D +
				ll * s->p->par_eSPD2D + mm * s->p->par_eSPA2D + nn * s->p->par_eSPA2V + oo * s->p->par_eSPD2V + pp * s->p->par_eSPV2V +
				qq * s->p->par_eSPB2V + rr * s->p->par_eSPC2V + tt * s->p->par_eSPB2D + uu * s->p->par_eSPC2D;

			w(inv) = exp(-(tde_va / 2 + s->p->par_eSPD) / (s->p->par_beta * s->p->par_temp));

		}
		else {
			w(inv) = 0.0;
		}
	}
}


void ParaBccJump::Choice(int& path, blitz::Array<double, 1> w)
{
	double sum = 0;

	for (int i = 0; i < s->n1nbr; i++) {
		sum = sum + w(i);
	}

	double ran = u(e) * sum;

	double boundry = 0;
	for (int i = 0; i < s->n1nbr; i++) {
		boundry = boundry + w(i);
		if (ran < boundry) {
			path = i;
			break;
		}
	}
}


void ParaBccJump::UpdateNeighbors(int vacx, int vacy, int vacz, int vacnm, int nf,
	blitz::Array<int, 1> x, blitz::Array<int, 1> y,
	blitz::Array<int, 1> z, blitz::Array<double, 1> w,int path)
{
	if (s->Site(nf, x(path), y(path), z(path)) == 0) {

		
		UpdateMutipleN1BeforeNeighbor(s, x(path), y(path), z(path), nf, vacnm, 0);

		UpdateMutipleN2BeforeNeighbor(s, x(path), y(path), z(path), nf, 0);

		UpdateMutipleN1AfterNeighbor(s, vacx, vacy, vacz, nf, vacnm, 0);

		UpdateMutipleN2AfterNeighbor(s, vacx, vacy, vacz, vacnm, 0);

	}
	else if (s->Site(nf, x(path), y(path), z(path)) == 1) {

		UpdateMutipleN1BeforeNeighbor(s, x(path), y(path), z(path), nf, vacnm, 1);

		UpdateMutipleN2BeforeNeighbor(s, x(path), y(path), z(path), nf, 1);

		UpdateMutipleN1AfterNeighbor(s, vacx, vacy, vacz, nf, vacnm, 1);

		UpdateMutipleN2AfterNeighbor(s, vacx, vacy, vacz, vacnm, 1);
		
	}
	else if (s->Site(nf, x(path), y(path), z(path)) == 3) {

		
		UpdateMutipleN1BeforeNeighbor(s, x(path), y(path), z(path), nf, vacnm, 3);

		UpdateMutipleN2BeforeNeighbor(s, x(path), y(path), z(path), nf, 3);

		UpdateMutipleN1AfterNeighbor(s, vacx, vacy, vacz, nf, vacnm, 3);

		UpdateMutipleN2AfterNeighbor(s, vacx, vacy, vacz, vacnm, 3);
	}
	else if (s->Site(nf, x(path), y(path), z(path)) == 4) {

		
		UpdateMutipleN1BeforeNeighbor(s, x(path), y(path), z(path), nf, vacnm, 4);

		UpdateMutipleN2BeforeNeighbor(s, x(path), y(path), z(path), nf, 4);

		UpdateMutipleN1AfterNeighbor(s, vacx, vacy, vacz, nf, vacnm, 4);

		UpdateMutipleN2AfterNeighbor(s, vacx, vacy, vacz, vacnm, 4);
	}
}



void ParaBccJump::ExchangeSites(int& vacx, int& vacy, int& vacz, int& vacnm, int nf,
	blitz::Array<int, 1> x, blitz::Array<int, 1> y,
	blitz::Array<int, 1> z, blitz::Array<double, 1> w, int path)
{

	int temp = s->Site(vacnm, vacx, vacy, vacz);
	s->Site(vacnm, vacx, vacy, vacz) = s->Site(nf, x(path), y(path), z(path));
	s->Site(nf, x(path), y(path), z(path)) = temp;
	 
	vacx = x(path); vacy = y(path); vacz = z(path); vacnm = nf;

}


/**
 * @brief 检查点是否越界并处理越界情况。
 *
 * @param vacx 空位的 x 坐标
 * @param vacy 空位的 y 坐标
 * @param vacz 空位的 z 坐标
 * @param vacnm 空位的类型
 * @param i 当前处理区域的 x 索引
 * @param j 当前处理区域的 y 索引
 * @param k 当前处理区域的 z 索引
 * @param idx x 方向的索引步长
 * @param idy y 方向的索引步长
 * @param idz z 方向的索引步长
 * @return 如果点位于当前处理区域内，则返回 false；如果点位于当前处理区域之外（越界），则返回 true。
 *
 * 该函数用于检查一个点是否越界。如果点位于当前处理区域内，它将返回 false，
 * 否则，它将返回 true，并表示该点越界。如果点越界，通常需要将它加入到越界区域的 vaclists 中，
 * 然后重新选择一个不越界的空位，如果没有可用的不越界空位，可能需要结束模拟进程。
 */
bool ParaBccJump::Border(int& vacx, int& vacy, int& vacz, int& vacnm,
	int i, int j, int k,int idx,int idy,int idz) {
	//检查是否越界
	//判断要交换的点是否越界，不越界，交换 
	//越界 把vac接入到越界的vaclists，然后重新选择一个空位，如果没有空位，此进程结束模拟。

	if(vacx >= i * idx && vacx < (i + 1) * idx && vacy >= j * idy && vacy < (j + 1) * idy && vacz >= k * idz && vacz < (k + 1) * idz){
		return false;
	}
	else {
		return true;
	}
}




