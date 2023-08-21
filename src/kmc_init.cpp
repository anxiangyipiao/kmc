
#include "kmc_init.h"


// 计算不同原子种类的数量
void Base::CalculateAtomNumbers(SimulationParameters& parameter) {
	// 计算不同原子种类的数量
	Bnum = int(parameter.nx * parameter.ny * parameter.nz * Length * parameter.par_compB);
	Cnum = int(parameter.nx * parameter.ny * parameter.nz * Length * parameter.par_compC);
	Dnum = int(parameter.nx * parameter.ny * parameter.nz * Length * parameter.par_compD);
	Vnum = parameter.Mutiple ? int(parameter.nx * parameter.ny * parameter.nz * Length * parameter.par_compV) : 1;
	Anum = parameter.nx * parameter.ny * parameter.nz * Length - Bnum - Vnum - Cnum - Dnum;

}
/*
函数名：GenerateLatticePoints
参数：SimulationParameters& parameter，包含模拟参数的结构体
功能：根据模拟参数，在晶格中生成指定数量的B、C、D、V原子，并将其添加到晶格点上
返回值：无
*/
void Base::GenerateLatticePoints(SimulationParameters& parameter) {
	int total_num = Bnum + Cnum + Dnum + Vnum;
	for (int i = 0; i < total_num;) {
		// 随机生成晶格点坐标和原子种类
		int ii = int(u(e) * parameter.nx);
		int jj = int(u(e) * parameter.ny);
		int zz = int(u(e) * parameter.nz);
		int nm = int(u(e) * Length);
		// 如果晶格点为空，则添加指定原子种类到晶格点上
		if (Site(nm, ii, jj, zz) == 0) {
			if (i < Bnum) {
				AddAtom(parameter,1, ii, jj, zz, nm);  // 添加B原子
			}
			else if (i < Bnum + Cnum) {
				AddAtom(parameter, 3, ii, jj, zz, nm);  // Add C atom
			}
			else if (i < Bnum + Cnum + Dnum) {
				AddAtom(parameter, 4, ii, jj, zz, nm);  // Add D atom
			}
			else {
				AddAtom(parameter, 2, ii, jj, zz, nm);  // Add V atom
			}
			i++;
		}
	}
}

void Base::AddAtom(SimulationParameters& parameter,int type, int ii, int jj, int zz, int nm) {
	// 设置晶格对应位置的原子类型
	Site(nm, ii, jj, zz) = type;
	// 如果添加的是V原子
	if (type == 2) { 
		// 如果开启了并行计算，将V原子添加到对应的子区域中
		if (parameter.Paralle) {
			// V atom
			int subregion_size_x = parameter.nx / parameter.parallelx;
			int subregion_size_y = parameter.ny / parameter.parallely;
			int subregion_size_z = parameter.nz / parameter.parallelz;

			int subregion_x = ii / subregion_size_x;
			int subregion_y = jj / subregion_size_y;
			int subregion_z = zz / subregion_size_z;
			VAC Vac = { ii, jj, zz, nm };
			Vaclists(subregion_x, subregion_y, subregion_z).push_back(Vac);
			vacx = ii;vacy = jj;vacz = zz;vacnm = nm;
		}
		else {
			// 如果未开启并行计算，直接记录V原子的位置
			vacx = ii;vacy = jj;vacz = zz;vacnm = nm;
		}		
	}
}

/**
 * Calculate double precision energy values Benergy and Aenergy, and compute the values of two 2D arrays BList and AList.
 *
 * @param parameter A structure containing simulation parameters.
 */
void Base::Energy(SimulationParameters& parameter)
{	
	double par_KBT = parameter.par_beta * parameter.par_temp;

	double  Benergy, Aenergy;

		for (int j = -n2nbr; j <= n2nbr; j++) {
			for (int i = -n1nbr; i <= n1nbr; i++) {
	
				Benergy = (i - 1) * parameter.par_eSPB1B - (i - 1) * parameter.par_eSPA1B - (i - 1) * parameter.par_eSPB1V + (i - 1) * parameter.par_eSPA1V +
					(j - 1) * parameter.par_eSPB2B - (j - 1) * parameter.par_eSPA2B - (j - 1) * parameter.par_eSPB2V + (j - 1) * parameter.par_eSPA2V;
				Aenergy = -i * parameter.par_eSPA1A + i * parameter.par_eSPA1V + i * parameter.par_eSPA1B - i * parameter.par_eSPB1V -
					j * parameter.par_eSPA2A + j * parameter.par_eSPA2B - j * parameter.par_eSPB2V + j * parameter.par_eSPA2V;
	
				double Del_B = exp(-(Benergy / 2 + parameter.par_eSPB) / par_KBT);
				double Del_A = exp(-(Aenergy / 2 + parameter.par_eSPA) / par_KBT);
	
				BList(i, j) = Del_B;
				AList(i, j) = Del_A;
			}
		}
}


 /**

 Initializes the states array based on the given simulation parameters.

 @param parameter The simulation parameters to be used for initializing the states array.
 */
void Base::InitStatesArray(SimulationParameters& parameter)
{
	// If read_file flag is set, read the states array from file and return
	if (parameter.read_file) {
		TofileStatesArray(parameter);
		return;
	}

	// If none of the above flags are set, create the states array using the AbStatesArray function
	AbInitStatesArrays(parameter);
}


//Reads the states array from a file specified by the filepath parameter.
void Base::TofileStatesArray(SimulationParameters& parameter)
{
		std::ifstream fin(parameter.filepath, std::ios::in);
		if (!fin.is_open())
		{
			std::cout << "无法找到这个文件！" << std::endl;
			return;
		}

		std::string line;
		while (getline(fin, line))
		{
			std::stringstream ss(line);
			int a, b, c, d, e;
			ss >> a >> b >> c >> d >> e;
			Site(d, a, b, c) = e;
			if (e == 2) {
				vacx = a;
				vacy = b;
				vacz = c;
				vacnm = d;
			}
		}
		fin.close();
	

}

// Creates a new state array by spreading B and V atoms randomly
//void Base::AbStatesArray(SimulationParameters& parameter)
//{
//	// STATE 2: vacancy, 0: A atom, 1: B atom 
//	// Calculate the number of B, V and A atoms
//	Bnum = int(parameter.nx * parameter.ny * parameter.nz * Length * parameter.par_compB);
//	Vnum = parameter.Mutiple ? int(parameter.nx * parameter.ny * parameter.nz * Length * parameter.par_compV) : 1;
//	Anum = parameter.nx * parameter.ny * parameter.nz * Length - Bnum - Vnum;
//
//	// spread B and V
//	for (int i = 0; i < Bnum + Vnum;) {
//		int ii = int(u(e) * parameter.nx);
//		int jj = int(u(e) * parameter.ny);
//		int zz = int(u(e) * parameter.nz);
//		int nm = int(u(e) * Length);
//
//		if (Site(nm, ii, jj, zz) == 0) {
//			if (i < Bnum) {
//				Site(nm, ii, jj, zz) = 1; // Add B atom
//			}
//			else {
//				Site(nm, ii, jj, zz) = 2; // Add V atom
//				vacx = ii;
//				vacy = jj;
//				vacz = zz;
//				vacnm = nm;
//			}
//			i++;
//		}
//	}
//
//
//}

void Base::AbInitStatesArrays(SimulationParameters& parameter)
{
	// STATE 2: vacancy, 0: A atom, 1: B atom, 3: C atom, 4: D atom
	// Calculate the number of atoms for each type
	CalculateAtomNumbers(parameter);

	GenerateLatticePoints(parameter);

}


std::vector<std::string> Base::Split(const std::string& s, const char& delim)
{
		std::vector<std::string> tokens;
		size_t lastPos = s.find_first_not_of(delim, 0);
		size_t pos = s.find(delim, lastPos);
		while (lastPos != std::string::npos) {
		    tokens.emplace_back(s.substr(lastPos, pos - lastPos));
		    lastPos = s.find_first_not_of(delim, pos);
		    pos = s.find(delim, lastPos);

		}
		return tokens;
		
}

void Base::Periodic(SimulationParameters& parameter)
{
	for (int i = -1; i <= parameter.nx; i++) {

		if (i == -1) {
			modla.push_back(-1 + parameter.nx);
		}
		else {
			modla.push_back(i % parameter.nx);
		}
	}
	for (int i = -1; i <= parameter.ny; i++) {

		if (i == -1) {
			modlb.push_back(-1 + parameter.ny);
		}
		else {
			modlb.push_back(i % parameter.ny);
		}
	}
	for (int i = -1; i <= parameter.nz; i++) {

		if (i == -1) {
			modlc.push_back(-1 + parameter.nz);
		}
		else {
			modlc.push_back(i % parameter.nz);
		}
	}

	modlx = &modla[1]; modly = &modlb[1]; modlz = &modlc[1];
}


// This function calculates the nearest neighbor atoms of each atom in the BCC lattice
// It first loops through each site in the lattice and checks its neighbors using the v1nbr_bcc and v2nbr_bcc arrays
// It adds 1 to the NN1 and NN2 arrays if a neighbor is found at that site

void Base::RecycleBNeighbors(SimulationParameters& parameter) {

	int nmn;
	#pragma omp parallel for private(nmn)
	for (int nm = 0; nm < Length; nm++) {
		// 根据nm值设置nmn
		nmn = (nm == 0) ? 1 : 0;

		for (int k = 0; k < parameter.nz; k++) {
			for (int j = 0; j < parameter.ny; j++) {
				for (int i = 0; i < parameter.nx; i++) {

					//NN1 数量
					std::vector<int> count1 = CalculateN1Neighbors(parameter,nm, nmn, i, j, k);
					//NN2 数量
					std::vector<int> count2 = CalculateN2Neighbors(parameter, nm, nmn, i, j, k);
					CountNeighbors(nm, i, j, k, count1, count2);

				}
			}
		}
	}
}
//
// fcc NN1 and NN2
void Base::RecycleFNeighbors(SimulationParameters& parameter) {

	int nmn = 0;
	for (int nm = 0; nm < 4; nm++) {
		for (int k = 0; k < parameter.nz; k++) {
			for (int j = 0; j < parameter.ny; j++) {
				for (int i = 0; i < parameter.nx; i++) {
					std::vector<int> count1 = CalculateN1FNeighbors(parameter, nm, nmn, i, j, k);
					std::vector<int> count2 = CalculateN2Neighbors(parameter, nm, nmn, i, j, k);
					CountNeighbors(nm, i, j, k, count1, count2);
				}
			}
		}
	}
}

// 此函数计算bcc晶格中给定点 (i, j, k,nm) 的最近邻点。
// 它接收一个 SimulationParameters 对象，三个整数作为坐标 (i, j, k,nm)，以及两个整数作为最近邻点的索引（nm 和 nmn）。
// 函数返回一个整型向量，其中第 i 个元素表示最近邻点中值等于 i 的点的数量。
std::vector<int> Base::CalculateN1Neighbors(SimulationParameters& parameter, int nm, int nmn, int i, int  j, int  k)
{
	std::vector<int> count(5, 0);
	for (int inv = 0; inv < 8; inv++) {
		int ii = modlx[i + v1nbr_bcc[nm][0][inv]];
		int jj = modly[j + v1nbr_bcc[nm][1][inv]];
		int kk = modlz[k + v1nbr_bcc[nm][2][inv]];

		int site_value = Site(nmn, ii, jj, kk);
		if (site_value >= 0 && site_value <= 4) {
			count[site_value]++;
		}
	}

	return count;
}
//计算fcc晶格中给定点 (i, j, k,nm) 的最近邻点。
std::vector<int> Base::CalculateN1FNeighbors(SimulationParameters& parameter, int nm, int nmn, int i, int  j, int  k)
{
	std::vector<int> count(5, 0);
	for (nmn = 0; nmn < 4; nmn++) {
		if (nmn != nm) {
			for (int inv = 0; inv < 4; inv++) {

				int ii = modlx[i + v1nbr_fcc[0][nm][nmn][inv]];
				int jj = modly[j + v1nbr_fcc[1][nm][nmn][inv]];
				int kk = modlz[k + v1nbr_fcc[2][nm][nmn][inv]];

				int site_value = Site(nmn, ii, jj, kk);
				if (site_value >= 0 && site_value <= 4) {
					count[site_value]++;
				}
			}
		}
	}

	return count;
}
//计算晶格中给定点 (i, j, k,nm) 的第二近邻点。
std::vector<int> Base::CalculateN2Neighbors(SimulationParameters& parameter, int nm, int nmn, int i, int  j, int  k)
{
		std::vector<int> count(5, 0);
		for (int inv = 0; inv < 6; inv++) {
			int ii = modlx[i + v2nbr_bcc[0][inv]];
			int jj = modly[j + v2nbr_bcc[1][inv]];
			int kk = modlz[k + v2nbr_bcc[2][inv]];

			int site_value = Site(nm, ii, jj, kk);
			if (site_value >= 0 && site_value <= 4) {
				count[site_value]++;
			}
		}
		return count;
}


void SingleBcc::CountNeighbors(int nm,int i, int  j, int  k, std::vector<int> count1, std::vector<int> count2)
{
	
	NN1(nm, i, j, k) = count1[0];

	NN2(nm, i, j, k) = count2[0];
	
}


void SingleFcc::CountNeighbors(int nm,int i, int j, int k, std::vector<int> count1, std::vector<int> count2)
{
	
	NN1(nm, i, j, k) = count1[0];

	NN2(nm, i, j, k) = count2[0];
							
}


void MultiBcc::CountNeighbors(int nm, int i, int j, int k, std::vector<int> count1, std::vector<int> count2)
{
	
	for (int l = 0; l < 3; l++) {
		NN1(l, nm, i, j, k) = count1[l];
	}
	for (int l = 0; l < 3; l++) {
		NN2(l, nm, i, j, k) = count2[l];
	}
}

void MultiFcc::CountNeighbors(int nm, int i, int j, int k, std::vector<int> count1, std::vector<int> count2)
{	
	for (int l = 0; l < 3; l++) {
		NN1(l, nm, i, j, k) = count1[l];
	}
	for (int l = 0; l < 3; l++) {
		NN2(l, nm, i, j, k) = count2[l];
	}
}


void MultiSizeBcc::CountNeighbors(int nm, int i, int j, int k, std::vector<int> count1, std::vector<int> count2)
{
	for (int l = 0; l < 5; l++) {
		NN1(l, nm, i, j, k) = count1[l];
	}
	for (int l = 0; l < 5; l++) {
		NN2(l, nm, i, j, k) = count2[l];
	}

}

void MultiSizeFcc::CountNeighbors(int nm, int i, int j, int k, std::vector<int> count1, std::vector<int> count2)
{

	for (int l = 0; l < 5; l++) {
		NN1(l, nm, i, j, k) = count1[l];
	}
	for (int l = 0; l < 5; l++) {
		NN2(l, nm, i, j, k) = count2[l];
	}

}

void ParaBcc::CountNeighbors(int nm, int i, int  j, int  k, std::vector<int> count1, std::vector<int> count2) {

	for (int l = 0; l < 5; l++) {
		NN1(l, nm, i, j, k) = count1[l];
	}
	for (int l = 0; l < 5; l++) {
		NN2(l, nm, i, j, k) = count2[l];
	}
}














