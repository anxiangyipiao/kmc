#include "kmc_init.h"


/**
 * @brief 计算不同原子种类的数量
 *
 * 此函数根据给定的模拟参数，计算不同原子种类（A、B、C、D、V）的数量。
 *
 * @param parameter 模拟参数对象，包含了模拟的各种配置信息。
 *
 * @details
 * 计算不同原子种类的数量使用以下公式：
 * - Bnum = int(parameter.nx * parameter.ny * parameter.nz * Length * parameter.par_compB)
 * - Cnum = int(parameter.nx * parameter.ny * parameter.nz * Length * parameter.par_compC)
 * - Dnum = int(parameter.nx * parameter.ny * parameter.nz * Length * parameter.par_compD)
 * - 如果 parameter.Mutiple 为 true，则 Vnum = int(parameter.nx * parameter.ny * parameter.nz * Length * parameter.par_compV)，
 *   否则 Vnum 为 1。
 * - Anum = parameter.nx * parameter.ny * parameter.nz * Length - Bnum - Vnum - Cnum - Dnum。
 *
 * 函数将计算结果存储在类成员变量 Bnum、Cnum、Dnum、Vnum 和 Anum 中。
 */
void Base::CalculateAtomNumbers(SimulationParameters& parameter) {

	Bnum = int(parameter.nx * parameter.ny * parameter.nz * Length * parameter.par_compB);
	Cnum = int(parameter.nx * parameter.ny * parameter.nz * Length * parameter.par_compC);
	Dnum = int(parameter.nx * parameter.ny * parameter.nz * Length * parameter.par_compD);
	Vnum = parameter.Mutiple ? int(parameter.nx * parameter.ny * parameter.nz * Length * parameter.par_compV) : 1;
	Anum = parameter.nx * parameter.ny * parameter.nz * Length - Bnum - Vnum - Cnum - Dnum;

}

/**
 * @brief 生成晶格点和不同原子种类
 *
 * 此函数根据给定的模拟参数，随机生成晶格点坐标和不同原子种类，并将它们添加到晶格中。
 *
 * @param parameter 模拟参数对象，包含了模拟的各种配置信息。
 *
 * @details
 * 函数的主要逻辑如下：
 * - 根据参数计算总原子数 total_num = Bnum + Cnum + Dnum + Vnum。
 * - 使用随机数生成晶格点坐标 (ii, jj, zz) 和原子种类 nm。
 * - 如果晶格点为空（Site(nm, ii, jj, zz) == 0），则将指定的原子种类添加到晶格点上。
 * - 根据计数器 i 的值，依次添加不同原子种类（B、C、D、V）到晶格中，直到达到总原子数 total_num。
 *
 * @param Site(nm, ii, jj, zz) 用于检查晶格点是否为空。
 * @param AddAtom(parameter, type, ii, jj, zz, nm) 用于将指定原子种类添加到晶格点上。
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


/**
 * @brief 添加原子到晶格中
 *
 * 此函数将指定类型的原子添加到晶格的指定位置，并根据模拟参数中的配置将V原子添加到子区域（如果开启了并行计算）。
 *
 * @param parameter 模拟参数对象，包含了模拟的各种配置信息。
 * @param type 要添加的原子类型，1 表示B原子，2 表示V原子，3 表示C原子，4 表示D原子。
 * @param ii 晶格点的 x 坐标。
 * @param jj 晶格点的 y 坐标。
 * @param zz 晶格点的 z 坐标。
 * @param nm 原子的编号或标识。
 *
 * @details
 * 函数执行以下操作：
 * - 设置晶格对应位置的原子类型，将原子类型 type 存储在 Site(nm, ii, jj, zz) 中。
 * - 如果添加的是 V 原子（type == 2）：
 *   - 如果开启了并行计算（parameter.Paralle == true），将 V 原子添加到对应的子区域中。
 *     - 计算子区域的大小 subregion_size_x、subregion_size_y 和 subregion_size_z。
 *     - 根据晶格点坐标 ii、jj 和 zz 计算子区域的索引 subregion_x、subregion_y 和 subregion_z。
 *     - 将 V 原子的位置信息记录在 Vaclists(subregion_x, subregion_y, subregion_z) 中。
 *   - 如果未开启并行计算，直接记录 V 原子的位置信息（vacx、vacy、vacz 和 vacnm）。
 */
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
 * @brief 计算单空位能量项 
 *
 * 此函数计算并存储晶格点之间的能量项，用于模拟中的物理过程。
 *
 * @param parameter 模拟参数对象，包含了模拟的各种配置信息。
 *
 * @details
 * 函数的主要逻辑如下：
 * - 根据参数计算 Boltzmann 常数乘以温度的值 par_KBT。
 * - 使用嵌套循环遍历晶格点的邻居点（范围由 n1nbr 和 n2nbr 决定）。
 * - 对于每个晶格点，计算 Benergy 和 Aenergy 两种能量项。
 * - 使用这些能量项计算 Del_B 和 Del_A，分别代表 B 和 A 原子的能量因子。
 * - 将 Del_B 和 Del_A 存储在 BList(i, j) 和 AList(i, j) 中，以便后续使用。
 *
 * @param Benergy 表示 B 原子的能量项。
 * @param Aenergy 表示 A 原子的能量项。
 * @param Del_B 表示 B 原子的能量因子。
 * @param Del_A 表示 A 原子的能量因子。
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
 * @brief 初始化晶格位点
 *
 * 此函数根据模拟参数配置初始化晶格位点。根据参数的 read_file 标志，
 * 可以选择从文件中读取晶格位点或使用 AbInitStatesArrays 函数创建晶格位点。
 *
 * @param parameter 模拟参数对象，包含了模拟的各种配置信息。
 *
 * @details
 * 函数的主要逻辑如下：
 * - 如果 parameter.read_file 标志为 true，则从文件中读取晶格位点并返回。
 * - 如果 parameter.read_file 标志为 false，则使用 AbInitStatesArrays 函数创建晶格位点。
 *
 * @param TofileStatesArray(parameter) 用于从文件中读取晶格位点。
 * @param AbInitStatesArrays(parameter) 用于创建晶格位点。
 */
void Base::InitStatesArray(SimulationParameters& parameter)
{
	
	if (parameter.read_file) {
		TofileStatesArray(parameter);
		return;
	}

	AbInitStatesArrays(parameter);
}


/**
 * @brief 从文件中读取晶格位点
 *
 * 此函数根据给定的文件路径从文件中读取晶格位点数组，并将其存储在相应的数据结构中。
 *
 * @param parameter 模拟参数对象，包含了模拟的各种配置信息，以及文件路径信息。
 *
 * @details
 * 函数的主要逻辑如下：
 * - 打开指定文件（parameter.filepath）以进行读取。
 * - 如果无法打开文件，输出错误消息并返回。
 * - 从文件逐行读取数据，每行包含五个整数。
 * - 解析每行数据，将对应的值分配给晶格点（Site）和 V 原子的位置信息（vacx、vacy、vacz 和 vacnm）。
 * - 关闭文件。
 *
 * @param Site(d, a, b, c) 用于存储从文件中读取的状态数组。
 * @param vacx, vacy, vacz, vacnm 用于存储从文件中读取的 V 原子的位置信息。
 */
void Base::TofileStatesArray(SimulationParameters& parameter)
{
		std::ifstream fin(parameter.filepath, std::ios::in);
		if (!fin.is_open())
		{
			//std::cout << "无法找到这个文件！" << std::endl;
			return ;
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


/**
 * @brief 使用 Ab 模型初始化晶格位点数组
 *
 * 此函数根据 Ab 模型的规则初始化晶格位点数组，包括设置晶格点的原子类型和位置。
 *
 * @param parameter 模拟参数对象，包含了模拟的各种配置信息。
 *
 * @details
 * 函数的主要逻辑如下：
 * - 根据模拟参数计算不同原子种类的数量。
 *  STATE 2: vacancy, 0: A atom, 1: B atom, 3: C atom, 4: D atom
 * - 使用 GenerateLatticePoints 函数随机生成晶格点坐标和不同原子种类，并将它们添加到晶格中。
 *
 * @param CalculateAtomNumbers(parameter) 用于计算不同原子种类的数量。
 * @param GenerateLatticePoints(parameter) 用于生成晶格点和不同原子种类。
 */
void Base::AbInitStatesArrays(SimulationParameters& parameter)
{
	// 计算不同原子种类的数量
	CalculateAtomNumbers(parameter);

	// 生成晶格点和不同原子种类
	GenerateLatticePoints(parameter);

}


/**
 * @brief 分割字符串为子串
 *
 * 此函数将输入的字符串按照指定的分隔符进行分割，返回一个包含子串的字符串向量。
 *
 * @param s 要分割的输入字符串。
 * @param delim 用于分割子串的分隔符。
 * @return 包含分割后子串的字符串向量。
 *
 * @details
 * 函数的主要逻辑如下：
 * - 从字符串 s 中查找第一个非分隔符字符的位置（lastPos）。
 * - 从 lastPos 位置开始查找下一个分隔符字符的位置（pos）。
 * - 如果找到分隔符，将 lastPos 和 pos 之间的子串添加到字符串向量 tokens 中。
 * - 更新 lastPos 和 pos 位置，继续查找下一个子串。
 * - 当无法找到更多分隔符时，结束循环并返回 tokens。
 *
 * @param tokens 用于存储分割后子串的字符串向量。
 */
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


/**
 * @brief 初始化周期性边界条件
 *
 * 此函数根据模拟参数配置初始化周期性边界条件的数据结构，用于处理晶格点的周期性边界情况。
 *
 * @param parameter 模拟参数对象，包含了模拟的各种配置信息。
 *
 * @details
 * 函数的主要逻辑如下：
 * - 使用循环分别处理 x、y 和 z 方向上的周期性边界条件。
 * - 对于每个方向，将 -1 到 parameter.nx/ny/nz 范围内的整数进行处理。
 *   - 如果当前整数为 -1，则将其映射为对应的周期性边界位置。
 *   - 否则，将当前整数取模操作后的结果添加到相应的 modla/modlb/modlc 向量中。
 * - 最后，将 modlx、modly 和 modlz 分别指向 modla、modlb 和 modlc 中的合适位置。
 *
 * @param modla, modlb, modlc 分别用于存储 x、y 和 z 方向上的周期性边界条件。
 * @param modlx, modly, modlz 分别指向对应方向上的周期性边界条件向量的起始位置。
 */
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

	// 将 modlx、modly、modlz 指向对应方向上的周期性边界条件向量的起始位置
	modlx = &modla[1]; modly = &modlb[1]; modlz = &modlc[1];
}

/**
 * @brief BCC结构  循环遍历 B 原子的邻居
 *
 * 此函数用于循环遍历 B 原子的邻居信息，包括计算 NN1 和 NN2 邻居数量，
 * 并将结果存储在适当的数据结构中，以便在模拟中使用。
 *
 * @param parameter 模拟参数对象，包含了模拟的各种配置信息。
 *
 * @details
 * 函数的主要逻辑如下：
 * - 使用 OpenMP 并行化来处理每个 nm 值（原子种类）。
 * - 根据 nm 值设置 nmn。
 * - 使用三重嵌套循环遍历晶格点的三个维度。
 * - 对于每个晶格点，计算 NN1 和 NN2 邻居数量。
 * - 调用 CountNeighbors 函数将计算结果存储在适当的数据结构中。
 *
 * @param nmn 用于表示当前原子种类的标志。
 * @param count1 存储 NN1 邻居数量的向量。
 * @param count2 存储 NN2 邻居数量的向量。
 * @param CountNeighbors(nm, i, j, k, count1, count2) 用于将计算结果存储在数据结构中。
 */
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
					
					// 将计算结果存储在适当的数据结构中
					CountNeighbors(nm, i, j, k, count1, count2);

				}
			}
		}
	}
}

//@brief FCC结构  循环遍历 B 原子的邻居
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


/**
 * @brief 此函数计算bcc晶格中给定点 (i, j, k,nm) 的最近邻点。
 *
 * 此函数用于计算指定晶格点的一邻居（NN1）的数量，根据不同原子种类和晶格位置。
 *
 * @param parameter 模拟参数对象，包含了模拟的各种配置信息。
 * @param nm 当前原子种类的标志，用于确定邻居的计算规则。
 * @param nmn 当前原子种类的标志（相邻种类），用于确定邻居的计算规则。
 * @param i 晶格点的 x 坐标。
 * @param j 晶格点的 y 坐标。
 * @param k 晶格点的 z 坐标。
 * @return 包含不同原子种类数量的向量，下标表示不同种类，值表示该种类的邻居数量。
 *
 * @details
 * 函数的主要逻辑如下：
 * - 使用循环迭代计算 8 个可能的 NN1 邻居位置。
 * - 根据当前原子种类（nm）和晶格位置（i, j, k），计算每个邻居的坐标。
 * - 查询邻居位置的原子种类（site_value）。
 * - 如果邻居位置上有合法的原子种类（0 到 4），则在相应的计数器中增加数量。
 * - 返回包含不同原子种类数量的向量，用于表示 NN1 邻居数量。
 *
 * @note
 * 此函数假设了存在名为 modlx、modly 和 modlz 的全局向量，用于处理周期性边界条件。
 * v1nbr_bcc 数组包含了不同原子种类的 NN1 邻居相对位置信息。
 */
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


/**
 * @brief 记录晶格点的 NN1 和 NN2 邻居数量
 *
 * 此函数用于记录指定晶格点的 NN1 和 NN2 邻居数量，将结果存储在适当的数据结构中。
 *
 * @param nm 当前原子种类的标志，用于确定邻居的计算规则。
 * @param i 晶格点的 x 坐标。
 * @param j 晶格点的 y 坐标。
 * @param k 晶格点的 z 坐标。
 * @param count1 包含 NN1 邻居数量的向量，下标表示不同种类，值表示该种类的 NN1 邻居数量。
 * @param count2 包含 NN2 邻居数量的向量，下标表示不同种类，值表示该种类的 NN2 邻居数量。
 *
 * @details
 * 函数的主要逻辑如下：
 * - 根据指定的原子种类（nm）和晶格位置（i, j, k）存储 NN1 和 NN2 邻居数量。
 * - 将 NN1 邻居数量存储在 NN1 数据结构中，用于后续的模拟计算。
 * - 将 NN2 邻居数量存储在 NN2 数据结构中，用于后续的模拟计算。
 *
 * @note
 * 此函数假设存在名为 NN1 和 NN2 的全局数据结构，用于存储晶格点的 NN1 和 NN2 邻居数量。
 * count1 和 count2 向量包含了不同原子种类的邻居数量信息。
 */
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














