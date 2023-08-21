#include "kmc_jump.h"


//JumpBase
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

//direct KMC
template<typename T>
void JumpBase::SumAndChoice(T* s)
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

//reject KMC
template<typename T>
void JumpBase::SumAndChoice2(T* s)
{
	double max_w = w(0);
	for (int i = 1; i < s->n1nbr; i++) {
		if (w(i) > max_w) {
			max_w = w(i);
		}
	}

	double r = u(e) * s->n1nbr;

	int j = r + 1;

	for (int i = 0; i < s->n1nbr; i++) {
		if (j - r < (w(i) / max_w)) {
			path = j - 1;
			break;
		}
	}

}

//单空位第二邻居 
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

//多空位BCC 第一
template<typename T>
void JumpBase::UpdateMutipleN1BeforeNeighbor(T* s,int x, int y, int z, int nf, int vacnm, int type)
{
	for (int ii = 0; ii < 8; ii++) {
		int xx = s->modlx[x + v1nbr_bcc[nf][0][ii]];
		int yy = s->modly[y + v1nbr_bcc[nf][1][ii]];
		int zz = s->modlz[z + v1nbr_bcc[nf][2][ii]];

		s->NN1(type, s->vacnm, xx, yy, zz) = s->NN1(type, s->vacnm, xx, yy, zz) - 1;
		s->NN1(2, s->vacnm, xx, yy, zz) = s->NN1(2, s->vacnm, xx, yy, zz) + 1;
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


template<typename T>
void JumpBase::ExchangeSite(T* s,int nf)
{
	s->Site(s->vacnm, s->vacx, s->vacy, s->vacz) = s->Site(nf, x(path), y(path), z(path));
	s->vacx = x(path); s->vacy = y(path); s->vacz = z(path); s->vacnm = nf;
	s->Site(s->vacnm, s->vacx, s->vacy, s->vacz) = 2;
}

template<typename T>
void JumpBase::runBSimulationLoop(int& nf, T* s)
{
	
		CalculatedBEnergy(nf, s);
		SumAndChoice(s);
		UpdateNeighbor();
		ExchangeSite(s,nf);

}

template<typename T>
void JumpBase::runFSimulationLoop(blitz::Array<int, 1>& fnf, T* s)
{

		CalculatedFEnergy(fnf, s);
		SumAndChoice(s);
		UpdateNeighbor();
		ExchangeSite(s, fnf(path));
	
}


// BCC 
void BccJump::CalculatedEnergy(int inv, int nm, int inv0)
{
	int dn1 = s->NN1(s->vacnm, s->vacx, s->vacy, s->vacz) - s->NN1(nf, x(inv), y(inv), z(inv));
	int dn2 = s->NN2(s->vacnm, s->vacx, s->vacy, s->vacz) - s->NN2(nf, x(inv), y(inv), z(inv));

	if (s->Site(nf, x(inv), y(inv), z(inv)) == 1) {

		w(inv) = (s->BList(dn1, dn2));
		//std::cout << w(inv) << std::endl;
	}
	else {

		w(inv) = (s->AList(dn1, dn2));
		//std::cout << w(inv) << std::endl;
	}
	
}

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


std::vector<int> MultiBccJump::ReturnList(int inv)
{
	int a = s->NN1(0, nf, x(inv), y(inv), z(inv)); int b = s->NN1(1, nf, x(inv), y(inv), z(inv)); int c = s->NN1(2, nf, x(inv), y(inv), z(inv));
	int d = s->NN1(0, s->vacnm, s->vacx, s->vacy, s->vacz); int e = s->NN1(1, s->vacnm, s->vacx, s->vacy, s->vacz); int f = s->NN1(2, s->vacnm, s->vacx, s->vacy, s->vacz);
	
	int aa = s->NN2(0, nf, x(inv), y(inv), z(inv)); int bb = s->NN2(1, nf, x(inv), y(inv), z(inv)); int cc = s->NN2(2, nf, x(inv), y(inv), z(inv));
	int dd = s->NN2(0, s->vacnm, s->vacx, s->vacy, s->vacz); int ee = s->NN2(1, s->vacnm, s->vacx, s->vacy, s->vacz); int ff = s->NN2(2, s->vacnm, s->vacx, s->vacy, s->vacz);

	std::vector<int> list = { a,b,c,d,e,f,aa,bb,cc,dd,ee,ff };
	return list;
}

//MultiBcc
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

void MultiBccJump::UpdateNeighbor()
{
	if (s->Site(nf, x(path), y(path), z(path)) == 0) {
		// 0 - Fe    NN1 里是 1 

		UpdateMutipleN1BeforeNeighbor(s, x(path), y(path), z(path), nf, s->vacnm, 1);
	
		UpdateMutipleN2BeforeNeighbor(s, x(path), y(path), z(path), nf, 1);
	
		UpdateMutipleN1AfterNeighbor(s, s->vacx, s->vacy, s->vacz, nf, s->vacnm, 1);

		UpdateMutipleN2AfterNeighbor(s, s->vacx, s->vacy, s->vacz, s->vacnm, 1);

	}
	else if (s->Site(nf, x(path), y(path), z(path)) == 1) {
		//1 - Cu    NN储存在0位置

		UpdateMutipleN1BeforeNeighbor(s, x(path), y(path), z(path), nf, s->vacnm, 0);

		UpdateMutipleN2BeforeNeighbor(s, x(path), y(path), z(path), nf, 0);

		UpdateMutipleN1AfterNeighbor(s, s->vacx, s->vacy, s->vacz, nf, s->vacnm, 0);
	
		UpdateMutipleN2AfterNeighbor(s, s->vacx, s->vacy, s->vacz, s->vacnm, 0);

	}
}

void MultiBccJump::runSimulationLoop()
{
	MultiBccJump::runBSimulationLoop(nf, s);
}

//MultiFcc
std::vector<int> MultiFccJump::ReturnList(int inv,int nm)
{
	int a = s->NN1(0, nm, x(inv), y(inv), z(inv)); int b = s->NN1(1, nm, x(inv), y(inv), z(inv)); int c = s->NN1(2, nm, x(inv), y(inv), z(inv));
	int d = s->NN1(0, s->vacnm, s->vacx, s->vacy, s->vacz); int e = s->NN1(1, s->vacnm, s->vacx, s->vacy, s->vacz); int f = s->NN1(2, s->vacnm, s->vacx, s->vacy, s->vacz);

	int aa = s->NN2(0, nm, x(inv), y(inv), z(inv)); int bb = s->NN2(1, nm, x(inv), y(inv), z(inv)); int cc = s->NN2(2, nm, x(inv), y(inv), z(inv));
	int dd = s->NN2(0, s->vacnm, s->vacx, s->vacy, s->vacz); int ee = s->NN2(1, s->vacnm, s->vacx, s->vacy, s->vacz); int ff = s->NN2(2, s->vacnm, s->vacx, s->vacy, s->vacz);

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
		UpdateMutipleFN1BeforeNeighbor(s, x(path), y(path), z(path), fnf(path),1);
	
		UpdateMutipleN2BeforeNeighbor(s, x(path), y(path), z(path), fnf(path), 1);
		
		UpdateMutipleFN1AfterNeighbor(s, s->vacx, s->vacy, s->vacz,s->vacnm, 1);
	
		UpdateMutipleN2AfterNeighbor(s, s->vacx, s->vacy, s->vacz, s->vacnm, 1);
	
	}
	else if (s->Site(fnf(path), x(path), y(path), z(path)) == 1)
	{
		UpdateMutipleFN1BeforeNeighbor(s, x(path), y(path), z(path), fnf(path), 0);

		UpdateMutipleN2BeforeNeighbor(s, x(path), y(path), z(path), fnf(path), 0);

		UpdateMutipleFN1AfterNeighbor(s, s->vacx, s->vacy, s->vacz, s->vacnm, 0);

		UpdateMutipleN2AfterNeighbor(s, s->vacx, s->vacy, s->vacz, s->vacnm, 0);
		
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
			int a = s->NN1(0, nf, x(inv), y(inv), z(inv)); int b = s->NN1(1, nf, x(inv), y(inv), z(inv)); int c = s->NN1(2, nf, x(inv), y(inv), z(inv));
			int h = s->NN1(3, nf, x(inv), y(inv), z(inv)); int i = s->NN1(4, nf, x(inv), y(inv), z(inv));
			int d = s->NN1(0, s->vacnm, s->vacx, s->vacy, s->vacz); int e = s->NN1(1, s->vacnm, s->vacx, s->vacy, s->vacz); int f = s->NN1(2, s->vacnm, s->vacx, s->vacy, s->vacz);
			int j = s->NN1(3, s->vacnm, s->vacx, s->vacy, s->vacz); int k = s->NN1(4, s->vacnm, s->vacx, s->vacy, s->vacz);

			int l = d - 1 - a; int m = e - b; int n = f + a + 2 - c - d; int o = b - e; int p = c - f - 1;
			int q = h - j; int r = i - k; int t = j - h; int u = k - i;

			int aa = s->NN2(0, nf, x(inv), y(inv), z(inv)); int bb = s->NN2(1, nf, x(inv), y(inv), z(inv)); int cc = s->NN2(2, nf, x(inv), y(inv), z(inv));
			int hh = s->NN2(3, nf, x(inv), y(inv), z(inv)); int ii = s->NN2(4, nf, x(inv), y(inv), z(inv));
			int dd = s->NN2(0, s->vacnm, s->vacx, s->vacy, s->vacz); int ee = s->NN2(1, s->vacnm, s->vacx, s->vacy, s->vacz); int ff = s->NN2(2, s->vacnm, s->vacx, s->vacy, s->vacz);
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
			int a = s->NN1(0, nf, x(inv), y(inv), z(inv)); int b = s->NN1(1, nf, x(inv), y(inv), z(inv)); int c = s->NN1(2, nf, x(inv), y(inv), z(inv));
			int h = s->NN1(3, nf, x(inv), y(inv), z(inv)); int i = s->NN1(4, nf, x(inv), y(inv), z(inv));
			int d = s->NN1(0, s->vacnm, s->vacx, s->vacy, s->vacz); int e = s->NN1(1, s->vacnm, s->vacx, s->vacy, s->vacz); int f = s->NN1(2, s->vacnm, s->vacx, s->vacy, s->vacz);
			int j = s->NN1(3, s->vacnm, s->vacx, s->vacy, s->vacz); int k = s->NN1(4, s->vacnm, s->vacx, s->vacy, s->vacz);

			int l = e - 1 - b; int m = d - a; int n = f + b + 2 - e - c; int o = a - d; int p = c - f - 1;
			int q = h - j; int r = i - k; int t = j - h; int u = k - i;

			int aa = s->NN2(0, nf, x(inv), y(inv), z(inv)); int bb = s->NN2(1, nf, x(inv), y(inv), z(inv)); int cc = s->NN2(2, nf, x(inv), y(inv), z(inv));
			int hh = s->NN2(3, nf, x(inv), y(inv), z(inv)); int ii = s->NN2(4, nf, x(inv), y(inv), z(inv));
			int dd = s->NN2(0, s->vacnm, s->vacx, s->vacy, s->vacz); int ee = s->NN2(1, s->vacnm, s->vacx, s->vacy, s->vacz); int ff = s->NN2(2, s->vacnm, s->vacx, s->vacy, s->vacz);
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
			int a = s->NN1(0, nf, x(inv), y(inv), z(inv)); int b = s->NN1(1, nf, x(inv), y(inv), z(inv)); int c = s->NN1(2, nf, x(inv), y(inv), z(inv));
			int h = s->NN1(3, nf, x(inv), y(inv), z(inv)); int i = s->NN1(4, nf, x(inv), y(inv), z(inv));
			int d = s->NN1(0, s->vacnm, s->vacx, s->vacy, s->vacz); int e = s->NN1(1, s->vacnm, s->vacx, s->vacy, s->vacz); int f = s->NN1(2, s->vacnm, s->vacx, s->vacy, s->vacz);
			int j = s->NN1(3, s->vacnm, s->vacx, s->vacy, s->vacz); int k = s->NN1(4, s->vacnm, s->vacx, s->vacy, s->vacz);

			int l = j - 1 - h; int m = e - b; int n = f + h + 2 - c - j; int o = b - e; int p = c - f - 1;
			int q = a - d; int r = i - k; int t = d - a; int u = k - i;

			int aa = s->NN2(0, nf, x(inv), y(inv), z(inv)); int bb = s->NN2(1, nf, x(inv), y(inv), z(inv)); int cc = s->NN2(2, nf, x(inv), y(inv), z(inv));
			int hh = s->NN2(3, nf, x(inv), y(inv), z(inv)); int ii = s->NN2(4, nf, x(inv), y(inv), z(inv));
			int dd = s->NN2(0, s->vacnm, s->vacx, s->vacy, s->vacz); int ee = s->NN2(1, s->vacnm, s->vacx, s->vacy, s->vacz); int ff = s->NN2(2, s->vacnm, s->vacx, s->vacy, s->vacz);
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
			int a = s->NN1(0, nf, x(inv), y(inv), z(inv)); int b = s->NN1(1, nf, x(inv), y(inv), z(inv)); int c = s->NN1(2, nf, x(inv), y(inv), z(inv));
			int h = s->NN1(3, nf, x(inv), y(inv), z(inv)); int i = s->NN1(4, nf, x(inv), y(inv), z(inv));
			int d = s->NN1(0, s->vacnm, s->vacx, s->vacy, s->vacz); int e = s->NN1(1, s->vacnm, s->vacx, s->vacy, s->vacz); int f = s->NN1(2, s->vacnm, s->vacx, s->vacy, s->vacz);
			int j = s->NN1(3, s->vacnm, s->vacx, s->vacy, s->vacz); int k = s->NN1(4, s->vacnm, s->vacx, s->vacy, s->vacz);

			int l = k - 1 - i; int m = e - b; int n = f + i + 2 - c - k; int o = b - e; int p = c - f - 1;
			int q = a - d; int r = h - j; int t = d - a; int u = j - h;

			int aa = s->NN2(0, nf, x(inv), y(inv), z(inv)); int bb = s->NN2(1, nf, x(inv), y(inv), z(inv)); int cc = s->NN2(2, nf, x(inv), y(inv), z(inv));
			int hh = s->NN2(3, nf, x(inv), y(inv), z(inv)); int ii = s->NN2(4, nf, x(inv), y(inv), z(inv));
			int dd = s->NN2(0, s->vacnm, s->vacx, s->vacy, s->vacz); int ee = s->NN2(1, s->vacnm, s->vacx, s->vacy, s->vacz); int ff = s->NN2(2, s->vacnm, s->vacx, s->vacy, s->vacz);
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

		UpdateMutipleN1BeforeNeighbor(s, x(path), y(path), z(path), nf, s->vacnm, 1);
	
		UpdateMutipleN2BeforeNeighbor(s, x(path), y(path), z(path), nf, 1);
	
		UpdateMutipleN1AfterNeighbor(s, s->vacx, s->vacy, s->vacz, nf, s->vacnm, 1);

		UpdateMutipleN2AfterNeighbor(s, s->vacx, s->vacy, s->vacz, s->vacnm, 1);
		
	}
	else if (s->Site(nf, x(path), y(path), z(path)) == 1) {

		UpdateMutipleN1BeforeNeighbor(s, x(path), y(path), z(path), nf, s->vacnm, 0);
	
		UpdateMutipleN2BeforeNeighbor(s, x(path), y(path), z(path), nf, 0);
	
		UpdateMutipleN1AfterNeighbor(s, s->vacx, s->vacy, s->vacz, nf, s->vacnm, 0);
	
		UpdateMutipleN2AfterNeighbor(s, s->vacx, s->vacy, s->vacz, s->vacnm, 0);
	
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
		int a = s->NN1(0, nf, x(inv), y(inv), z(inv)); int b = s->NN1(1, nf, x(inv), y(inv), z(inv)); int c = s->NN1(2, nf, x(inv), y(inv), z(inv));
		int h = s->NN1(3, nf, x(inv), y(inv), z(inv)); int i = s->NN1(4, nf, x(inv), y(inv), z(inv));
		int d = s->NN1(0, s->vacnm, s->vacx, s->vacy, s->vacz); int e = s->NN1(1, s->vacnm, s->vacx, s->vacy, s->vacz); int f = s->NN1(2, s->vacnm, s->vacx, s->vacy, s->vacz);
		int j = s->NN1(3, s->vacnm, s->vacx, s->vacy, s->vacz); int k = s->NN1(4, s->vacnm, s->vacx, s->vacy, s->vacz);

		int l = d - 1 - a; int m = e - b; int n = f + a + 2 - c - d; int o = b - e; int p = c - f - 1;
		int q = h - j; int r = i - k; int t = j - h; int u = k - i;

		int aa = s->NN2(0, nf, x(inv), y(inv), z(inv)); int bb = s->NN2(1, nf, x(inv), y(inv), z(inv)); int cc = s->NN2(2, nf, x(inv), y(inv), z(inv));
		int hh = s->NN2(3, nf, x(inv), y(inv), z(inv)); int ii = s->NN2(4, nf, x(inv), y(inv), z(inv));
		int dd = s->NN2(0, s->vacnm, s->vacx, s->vacy, s->vacz); int ee = s->NN2(1, s->vacnm, s->vacx, s->vacy, s->vacz); int ff = s->NN2(2, s->vacnm, s->vacx, s->vacy, s->vacz);
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
		int a = s->NN1(0, nf, x(inv), y(inv), z(inv)); int b = s->NN1(1, nf, x(inv), y(inv), z(inv)); int c = s->NN1(2, nf, x(inv), y(inv), z(inv));
		int h = s->NN1(3, nf, x(inv), y(inv), z(inv)); int i = s->NN1(4, nf, x(inv), y(inv), z(inv));
		int d = s->NN1(0, s->vacnm, s->vacx, s->vacy, s->vacz); int e = s->NN1(1, s->vacnm, s->vacx, s->vacy, s->vacz); int f = s->NN1(2, s->vacnm, s->vacx, s->vacy, s->vacz);
		int j = s->NN1(3, s->vacnm, s->vacx, s->vacy, s->vacz); int k = s->NN1(4, s->vacnm, s->vacx, s->vacy, s->vacz);

		int l = e - 1 - b; int m = d - a; int n = f + b + 2 - e - c; int o = a - d; int p = c - f - 1;
		int q = h - j; int r = i - k; int t = j - h; int u = k - i;

		int aa = s->NN2(0, nf, x(inv), y(inv), z(inv)); int bb = s->NN2(1, nf, x(inv), y(inv), z(inv)); int cc = s->NN2(2, nf, x(inv), y(inv), z(inv));
		int hh = s->NN2(3, nf, x(inv), y(inv), z(inv)); int ii = s->NN2(4, nf, x(inv), y(inv), z(inv));
		int dd = s->NN2(0, s->vacnm, s->vacx, s->vacy, s->vacz); int ee = s->NN2(1, s->vacnm, s->vacx, s->vacy, s->vacz); int ff = s->NN2(2, s->vacnm, s->vacx, s->vacy, s->vacz);
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
		int a = s->NN1(0, nf, x(inv), y(inv), z(inv)); int b = s->NN1(1, nf, x(inv), y(inv), z(inv)); int c = s->NN1(2, nf, x(inv), y(inv), z(inv));
		int h = s->NN1(3, nf, x(inv), y(inv), z(inv)); int i = s->NN1(4, nf, x(inv), y(inv), z(inv));
		int d = s->NN1(0, s->vacnm, s->vacx, s->vacy, s->vacz); int e = s->NN1(1, s->vacnm, s->vacx, s->vacy, s->vacz); int f = s->NN1(2, s->vacnm, s->vacx, s->vacy, s->vacz);
		int j = s->NN1(3, s->vacnm, s->vacx, s->vacy, s->vacz); int k = s->NN1(4, s->vacnm, s->vacx, s->vacy, s->vacz);

		int l = j - 1 - h; int m = e - b; int n = f + h + 2 - c - j; int o = b - e; int p = c - f - 1;
		int q = a - d; int r = i - k; int t = d - a; int u = k - i;

		int aa = s->NN2(0, nf, x(inv), y(inv), z(inv)); int bb = s->NN2(1, nf, x(inv), y(inv), z(inv)); int cc = s->NN2(2, nf, x(inv), y(inv), z(inv));
		int hh = s->NN2(3, nf, x(inv), y(inv), z(inv)); int ii = s->NN2(4, nf, x(inv), y(inv), z(inv));
		int dd = s->NN2(0, s->vacnm, s->vacx, s->vacy, s->vacz); int ee = s->NN2(1, s->vacnm, s->vacx, s->vacy, s->vacz); int ff = s->NN2(2, s->vacnm, s->vacx, s->vacy, s->vacz);
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
		int a = s->NN1(0, nf, x(inv), y(inv), z(inv)); int b = s->NN1(1, nf, x(inv), y(inv), z(inv)); int c = s->NN1(2, nf, x(inv), y(inv), z(inv));
		int h = s->NN1(3, nf, x(inv), y(inv), z(inv)); int i = s->NN1(4, nf, x(inv), y(inv), z(inv));
		int d = s->NN1(0, s->vacnm, s->vacx, s->vacy, s->vacz); int e = s->NN1(1, s->vacnm, s->vacx, s->vacy, s->vacz); int f = s->NN1(2, s->vacnm, s->vacx, s->vacy, s->vacz);
		int j = s->NN1(3, s->vacnm, s->vacx, s->vacy, s->vacz); int k = s->NN1(4, s->vacnm, s->vacx, s->vacy, s->vacz);

		int l = k - 1 - i; int m = e - b; int n = f + i + 2 - c - k; int o = b - e; int p = c - f - 1;
		int q = a - d; int r = h - j; int t = d - a; int u = j - h;

		int aa = s->NN2(0, nf, x(inv), y(inv), z(inv)); int bb = s->NN2(1, nf, x(inv), y(inv), z(inv)); int cc = s->NN2(2, nf, x(inv), y(inv), z(inv));
		int hh = s->NN2(3, nf, x(inv), y(inv), z(inv)); int ii = s->NN2(4, nf, x(inv), y(inv), z(inv));
		int dd = s->NN2(0, s->vacnm, s->vacx, s->vacy, s->vacz); int ee = s->NN2(1, s->vacnm, s->vacx, s->vacy, s->vacz); int ff = s->NN2(2, s->vacnm, s->vacx, s->vacy, s->vacz);
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
		UpdateMutipleFN1BeforeNeighbor(s, x(path), y(path), z(path), fnf(path), 1);
		
		UpdateMutipleN2BeforeNeighbor(s, x(path), y(path), z(path), fnf(path), 1);
		
		UpdateMutipleFN1AfterNeighbor(s, s->vacx, s->vacy, s->vacz, s->vacnm, 1);
	
		UpdateMutipleN2AfterNeighbor(s, s->vacx, s->vacy, s->vacz, s->vacnm, 1);

	}
	else if (s->Site(fnf(path), x(path), y(path), z(path)) == 1)
	{
		UpdateMutipleFN1BeforeNeighbor(s, x(path), y(path), z(path), fnf(path), 0);
		
		UpdateMutipleN2BeforeNeighbor(s, x(path), y(path), z(path), fnf(path), 0);
	
		UpdateMutipleFN1AfterNeighbor(s, s->vacx, s->vacy, s->vacz, s->vacnm, 0);
		
		UpdateMutipleN2AfterNeighbor(s, s->vacx, s->vacy, s->vacz, s->vacnm, 0);

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


void ParaBccJump::Sector(int i, int j, int k, int idx, int idy, int idz) {

	double t1 = 0.001;
	double t2 = 0;

	// 随机选择s->Vaclists(i, j, k)中的一个进行mc
	int randomNumber1 = rand() % (s->Vaclists(i, j, k).size());

	auto& vac = s->Vaclists(i, j, k)[randomNumber1];

	//while (t1>=t2) {
		//for (auto& vac : s->Vaclists(i, j, k)) {
			// 处理每个 VAC
		int vacx = vac.x;
		int vacy = vac.y;
		int vacz = vac.z;
		int vacnm = vac.nm;

		int nf = 0; int path = 0;//局部变量优先级更高
		blitz::Array<int, 1> x, y, z; blitz::Array<double, 1> w;
		x.resize(8); y.resize(8); z.resize(8); w.resize(8);
		x = 0; y = 0; z = 0; w = 0;


		for (int step = 0; step < 1000; step++) {

			nf = (vac.nm == 0) ? 1 : 0;

    		CalculatedEnergys(vacx, vacy, vacz, vacnm, nf, x, y, z, w);
			Choice(path, w);
			UpdateNeighbors(vacx, vacy, vacz, vacnm, nf, x, y, z, w, path);
			ExchangeSites(vac.x, vac.y, vac.z, vac.nm, nf, x, y, z, w, path);

			//如果空位跳出区域
			if (Border(vac.x, vac.y, vac.z, vac.nm, i, j, k, idx, idy, idz)) {
				
				int nexti = vac.x / idx; int nextj = vac.y / idy; int nextk = vac.z / idz;
				//add vac
				std::cout << "before" << std::endl;
				for (int a = 0; a < s->Vaclists(nexti, nextj, nextk).size(); a++) {
					
					std::cout << s->Vaclists(nexti, nextj, nextk)[a].x << " " <<
						s->Vaclists(nexti, nextj, nextk)[a].y << " " <<
						s->Vaclists(nexti, nextj, nextk)[a].z << " " <<
						s->Vaclists(nexti, nextj, nextk)[a].nm << std::endl;
				}
				std::cout << "after" << std::endl;
				s->Vaclists(nexti, nextj, nextk).push_back(vac);
				
				for (int a = 0; a < s->Vaclists(nexti, nextj, nextk).size(); a++) {

					std::cout << s->Vaclists(nexti, nextj, nextk)[a].x << " " <<
						s->Vaclists(nexti, nextj, nextk)[a].y << " " <<
						s->Vaclists(nexti, nextj, nextk)[a].z << " " <<
						s->Vaclists(nexti, nextj, nextk)[a].nm << std::endl;
				}

				//del vac
				std::cout << "before" << std::endl;
				for (int a = 0; a < s->Vaclists(i, j, k).size(); a++) {

					std::cout << s->Vaclists(i, j, k)[a].x << " " <<
						s->Vaclists(i, j, k)[a].y << " " <<
						s->Vaclists(i, j, k)[a].z << " " <<
						s->Vaclists(i, j, k)[a].nm << std::endl;
				}
				//del vac
				auto it = std::find_if(s->Vaclists(i, j, k).begin(), s->Vaclists(i, j, k).end(),
					[&](const Base::VAC& other) { return vac.x == other.x && vac.y == other.y && vac.z == other.z && vac.nm == other.nm; });

				s->Vaclists(i, j, k).erase(it);

				std::cout << "after" << std::endl;
				for (int a = 0; a < s->Vaclists(i, j, k).size(); a++) {

					std::cout << s->Vaclists(i, j, k)[a].x << " " <<
						s->Vaclists(i, j, k)[a].y << " " <<
						s->Vaclists(i, j, k)[a].z << " " <<
						s->Vaclists(i, j, k)[a].nm << std::endl;
				}

				// 选择一个新的vac
				if (s->Vaclists(i, j, k).empty()) {
					break;
				}
				int randomNumber2 = rand() % (s->Vaclists(i, j, k).size());
				vac = s->Vaclists(i, j, k)[randomNumber2];

				vacx = vac.x;
				vacy = vac.y;
				vacz = vac.z;
				vacnm = vac.nm;

			};

			//t2 = t2 + 0.001;
		}
		//}
	//}

}


void ParaBccJump::Divide()
{	

	int idx = s->p->nx / s->p->parallelx;
	int idy = s->p->ny / s->p->parallely;
	int idz = s->p->nz / s->p->parallelz;

	#pragma omp parallel for
	for (int i = 0; i < s->p->parallelx; i += 2) {
		for (int j = 0; j < s->p->parallely; j += 2) {
			for (int k = 0; k < s->p->parallelz; k += 2) {
				
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
}
//
void ParaBccJump::CalculatedEnergys(int vacx, int vacy, int vacz, int vacnm, int& nf,
									blitz::Array<int, 1> x, blitz::Array<int, 1> y,
									blitz::Array<int, 1> z, blitz::Array<double, 1> w)
{
	for (int inv = 0; inv < 8; inv++) {

		x(inv) = s->modlx[vacx + v1nbr_bcc[vacnm][0][inv]];
		y(inv) = s->modly[vacy + v1nbr_bcc[vacnm][1][inv]];
		z(inv) = s->modlz[vacz + v1nbr_bcc[vacnm][2][inv]];

		if (s->Site(nf, x(inv), y(inv), z(inv)) == 1)
		{
			
			int a = s->NN1(0, nf, x(inv), y(inv), z(inv)); int b = s->NN1(1, nf, x(inv), y(inv), z(inv)); int c = s->NN1(2, nf, x(inv), y(inv), z(inv));
			int h = s->NN1(3, nf, x(inv), y(inv), z(inv)); int i = s->NN1(4, nf, x(inv), y(inv), z(inv));
			int d = s->NN1(0, vacnm, vacx, vacy, vacz); int e = s->NN1(1, vacnm, vacx, vacy, vacz); int f = s->NN1(2, vacnm, vacx, vacy, vacz);
			int j = s->NN1(3, vacnm, vacx, vacy, vacz); int k = s->NN1(4, vacnm,vacx, vacy, vacz);

			int l = d - 1 - a; int m = e - b; int n = f + a + 2 - c - d; int o = b - e; int p = c - f - 1;
			int q = h - j; int r = i - k; int t = j - h; int u = k - i;

			int aa = s->NN2(0, nf, x(inv), y(inv), z(inv)); int bb = s->NN2(1, nf, x(inv), y(inv), z(inv)); int cc = s->NN2(2, nf, x(inv), y(inv), z(inv));
			int hh = s->NN2(3, nf, x(inv), y(inv), z(inv)); int ii = s->NN2(4, nf, x(inv), y(inv), z(inv));
			int dd = s->NN2(0,vacnm, vacx, vacy, vacz); int ee = s->NN2(1, vacnm, vacx, vacy, vacz); int ff = s->NN2(2, vacnm, vacx, vacy, vacz);
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
			int a = s->NN1(0, nf, x(inv), y(inv), z(inv)); int b = s->NN1(1, nf, x(inv), y(inv), z(inv)); int c = s->NN1(2, nf, x(inv), y(inv), z(inv));
			int h = s->NN1(3, nf, x(inv), y(inv), z(inv)); int i = s->NN1(4, nf, x(inv), y(inv), z(inv));
			int d = s->NN1(0, vacnm, vacx, vacy, vacz); int e = s->NN1(1, vacnm, vacx, vacy, vacz); int f = s->NN1(2, vacnm, vacx, vacy, vacz);
			int j = s->NN1(3, vacnm, vacx, vacy, vacz); int k = s->NN1(4, vacnm, vacx, vacy, vacz);

			int l = e - 1 - b; int m = d - a; int n = f + b + 2 - e - c; int o = a - d; int p = c - f - 1;
			int q = h - j; int r = i - k; int t = j - h; int u = k - i;

			int aa = s->NN2(0, nf, x(inv), y(inv), z(inv)); int bb = s->NN2(1, nf, x(inv), y(inv), z(inv)); int cc = s->NN2(2, nf, x(inv), y(inv), z(inv));
			int hh = s->NN2(3, nf, x(inv), y(inv), z(inv)); int ii = s->NN2(4, nf, x(inv), y(inv), z(inv));
			int dd = s->NN2(0, vacnm, vacx, vacy, vacz); int ee = s->NN2(1, vacnm, vacx, vacy, vacz); int ff = s->NN2(2, vacnm, vacx, vacy, vacz);
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
			int a = s->NN1(0, nf, x(inv), y(inv), z(inv)); int b = s->NN1(1, nf, x(inv), y(inv), z(inv)); int c = s->NN1(2, nf, x(inv), y(inv), z(inv));
			int h = s->NN1(3, nf, x(inv), y(inv), z(inv)); int i = s->NN1(4, nf, x(inv), y(inv), z(inv));
			int d = s->NN1(0, vacnm, vacx, vacy, vacz); int e = s->NN1(1, vacnm, vacx, vacy, vacz); int f = s->NN1(2, vacnm, vacx, vacy, vacz);
			int j = s->NN1(3, vacnm, vacx, vacy, vacz); int k = s->NN1(4, vacnm, vacx, vacy, vacz);

			int l = j - 1 - h; int m = e - b; int n = f + h + 2 - c - j; int o = b - e; int p = c - f - 1;
			int q = a - d; int r = i - k; int t = d - a; int u = k - i;

			int aa = s->NN2(0, nf, x(inv), y(inv), z(inv)); int bb = s->NN2(1, nf, x(inv), y(inv), z(inv)); int cc = s->NN2(2, nf, x(inv), y(inv), z(inv));
			int hh = s->NN2(3, nf, x(inv), y(inv), z(inv)); int ii = s->NN2(4, nf, x(inv), y(inv), z(inv));
			int dd = s->NN2(0, vacnm, vacx, vacy, vacz); int ee = s->NN2(1, vacnm, vacx, vacy, vacz); int ff = s->NN2(2, vacnm, vacx, vacy, vacz);
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
			int a = s->NN1(0, nf, x(inv), y(inv), z(inv)); int b = s->NN1(1, nf, x(inv), y(inv), z(inv)); int c = s->NN1(2, nf, x(inv), y(inv), z(inv));
			int h = s->NN1(3, nf, x(inv), y(inv), z(inv)); int i = s->NN1(4, nf, x(inv), y(inv), z(inv));
			int d = s->NN1(0, vacnm, vacx, vacy, vacz); int e = s->NN1(1, vacnm, vacx, vacy, vacz); int f = s->NN1(2, vacnm, vacx, vacy, vacz);
			int j = s->NN1(3, vacnm, vacx, vacy, vacz); int k = s->NN1(4, vacnm, vacx, vacy, vacz);

			int l = k - 1 - i; int m = e - b; int n = f + i + 2 - c - k; int o = b - e; int p = c - f - 1;
			int q = a - d; int r = h - j; int t = d - a; int u = j - h;

			int aa = s->NN2(0, nf, x(inv), y(inv), z(inv)); int bb = s->NN2(1, nf, x(inv), y(inv), z(inv)); int cc = s->NN2(2, nf, x(inv), y(inv), z(inv));
			int hh = s->NN2(3, nf, x(inv), y(inv), z(inv)); int ii = s->NN2(4, nf, x(inv), y(inv), z(inv));
			int dd = s->NN2(0, vacnm, vacx, vacy, vacz); int ee = s->NN2(1, vacnm, vacx, vacy, vacz); int ff = s->NN2(2, vacnm, vacx, vacy, vacz);
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
//
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
//
//
void ParaBccJump::UpdateNeighbors(int vacx, int vacy, int vacz, int vacnm, int& nf,
	blitz::Array<int, 1> x, blitz::Array<int, 1> y,
	blitz::Array<int, 1> z, blitz::Array<double, 1> w,int& path)
{
	if (s->Site(nf, x(path), y(path), z(path)) == 0) {

		for (int ii = 0; ii < 8; ii++) {
			int xx = s->modlx[x(path) + v1nbr_bcc[nf][0][ii]];
			int yy = s->modly[y(path) + v1nbr_bcc[nf][1][ii]];
			int zz = s->modlz[z(path) + v1nbr_bcc[nf][2][ii]];

			s->NN1(1, vacnm, xx, yy, zz) = s->NN1(1,vacnm, xx, yy, zz) - 1;
			s->NN1(2, vacnm, xx, yy, zz) = s->NN1(2,vacnm, xx, yy, zz) + 1;
		}


		for (int jj = 0; jj < 6; jj++) {
			int xx = s->modlx[x(path) + v2nbr_bcc[0][jj]];
			int	yy = s->modly[y(path) + v2nbr_bcc[1][jj]];
			int	zz = s->modlz[z(path) + v2nbr_bcc[2][jj]];

			s->NN2(1, nf, xx, yy, zz) = s->NN2(1, nf, xx, yy, zz) - 1;
			s->NN2(2, nf, xx, yy, zz) = s->NN2(2, nf, xx, yy, zz) + 1;
		}



		for (int ii = 0; ii < 8; ii++) {
			int xx = s->modlx[vacx + v1nbr_bcc[vacnm][0][ii]];
			int yy = s->modly[vacy + v1nbr_bcc[vacnm][1][ii]];
			int zz = s->modlz[vacz + v1nbr_bcc[vacnm][2][ii]];

			s->NN1(1, nf, xx, yy, zz) = s->NN1(1, nf, xx, yy, zz) + 1;
			s->NN1(2, nf, xx, yy, zz) = s->NN1(2, nf, xx, yy, zz) - 1;
		}


		for (int jj = 0; jj < 6; jj++) {
			int xx = s->modlx[vacx + v2nbr_bcc[0][jj]];
			int	yy = s->modly[vacy + v2nbr_bcc[1][jj]];
			int	zz = s->modlz[vacz + v2nbr_bcc[2][jj]];

			s->NN2(1, vacnm, xx, yy, zz) = s->NN2(1, vacnm, xx, yy, zz) + 1;
			s->NN2(2, vacnm, xx, yy, zz) = s->NN2(2, vacnm, xx, yy, zz) - 1;

		}
	}
	else if (s->Site(nf, x(path), y(path), z(path)) == 1) {

		for (int ii = 0; ii < 8; ii++) {

			int xx = s->modlx[x(path) + v1nbr_bcc[nf][0][ii]];

			int yy = s->modly[y(path) + v1nbr_bcc[nf][1][ii]];

			int zz = s->modlz[z(path) + v1nbr_bcc[nf][2][ii]];


			s->NN1(0, vacnm, xx, yy, zz) = s->NN1(0, vacnm, xx, yy, zz) - 1;

			s->NN1(2, vacnm, xx, yy, zz) = s->NN1(2, vacnm, xx, yy, zz) + 1;

		}


		for (int jj = 0; jj < 6; jj++) {
			int xx = s->modlx[x(path) + v2nbr_bcc[0][jj]];
			int	yy = s->modly[y(path) + v2nbr_bcc[1][jj]];
			int	zz = s->modlz[z(path) + v2nbr_bcc[2][jj]];

			s->NN2(0, nf, xx, yy, zz) = s->NN2(0, nf, xx, yy, zz) - 1;
			s->NN2(2, nf, xx, yy, zz) = s->NN2(2, nf, xx, yy, zz) + 1;

		}



		for (int ii = 0; ii < 8; ii++) {
			int xx = s->modlx[vacx + v1nbr_bcc[vacnm][0][ii]];
			int yy = s->modly[vacy + v1nbr_bcc[vacnm][1][ii]];
			int zz = s->modlz[vacz + v1nbr_bcc[vacnm][2][ii]];

			s->NN1(0, nf, xx, yy, zz) = s->NN1(0, nf, xx, yy, zz) + 1;
			s->NN1(2, nf, xx, yy, zz) = s->NN1(2, nf, xx, yy, zz) - 1;

		}


		for (int jj = 0; jj < 6; jj++) {
			int xx = s->modlx[vacx + v2nbr_bcc[0][jj]];
			int	yy = s->modly[vacy + v2nbr_bcc[1][jj]];
			int	zz = s->modlz[vacz + v2nbr_bcc[2][jj]];

			s->NN2(0, vacnm, xx, yy, zz) = s->NN2(0, vacnm, xx, yy, zz) + 1;
			s->NN2(2, vacnm, xx, yy, zz) = s->NN2(2, vacnm, xx, yy, zz) - 1;
		}
	}
	else if (s->Site(nf, x(path), y(path), z(path)) == 3) {

		for (int ii = 0; ii < 8; ii++) {

			int xx = s->modlx[x(path) + v1nbr_bcc[nf][0][ii]];

			int yy = s->modly[y(path) + v1nbr_bcc[nf][1][ii]];

			int zz = s->modlz[z(path) + v1nbr_bcc[nf][2][ii]];


			s->NN1(3, vacnm, xx, yy, zz) = s->NN1(3, vacnm, xx, yy, zz) - 1;

			s->NN1(2, vacnm, xx, yy, zz) = s->NN1(2, vacnm, xx, yy, zz) + 1;

		}


		for (int jj = 0; jj < 6; jj++) {
			int xx = s->modlx[x(path) + v2nbr_bcc[0][jj]];
			int	yy = s->modly[y(path) + v2nbr_bcc[1][jj]];
			int	zz = s->modlz[z(path) + v2nbr_bcc[2][jj]];

			s->NN2(3, nf, xx, yy, zz) = s->NN2(3, nf, xx, yy, zz) - 1;
			s->NN2(2, nf, xx, yy, zz) = s->NN2(2, nf, xx, yy, zz) + 1;

		}



		for (int ii = 0; ii < 8; ii++) {
			int xx = s->modlx[vacx + v1nbr_bcc[vacnm][0][ii]];
			int yy = s->modly[vacy + v1nbr_bcc[vacnm][1][ii]];
			int zz = s->modlz[vacz + v1nbr_bcc[vacnm][2][ii]];

			s->NN1(3, nf, xx, yy, zz) = s->NN1(3, nf, xx, yy, zz) + 1;
			s->NN1(2, nf, xx, yy, zz) = s->NN1(2, nf, xx, yy, zz) - 1;

		}


		for (int jj = 0; jj < 6; jj++) {
			int xx = s->modlx[s->vacx + v2nbr_bcc[0][jj]];
			int	yy = s->modly[s->vacy + v2nbr_bcc[1][jj]];
			int	zz = s->modlz[s->vacz + v2nbr_bcc[2][jj]];

			s->NN2(3, vacnm, xx, yy, zz) = s->NN2(3, vacnm, xx, yy, zz) + 1;
			s->NN2(2, vacnm, xx, yy, zz) = s->NN2(2, vacnm, xx, yy, zz) - 1;
		}
	}
	else if (s->Site(nf, x(path), y(path), z(path)) == 4) {

		for (int ii = 0; ii < 8; ii++) {

			int xx = s->modlx[x(path) + v1nbr_bcc[nf][0][ii]];

			int yy = s->modly[y(path) + v1nbr_bcc[nf][1][ii]];

			int zz = s->modlz[z(path) + v1nbr_bcc[nf][2][ii]];


			s->NN1(4, vacnm, xx, yy, zz) = s->NN1(4, vacnm, xx, yy, zz) - 1;

			s->NN1(2, vacnm, xx, yy, zz) = s->NN1(2, vacnm, xx, yy, zz) + 1;

		}


		for (int jj = 0; jj < 6; jj++) {
			int xx = s->modlx[x(path) + v2nbr_bcc[0][jj]];
			int	yy = s->modly[y(path) + v2nbr_bcc[1][jj]];
			int	zz = s->modlz[z(path) + v2nbr_bcc[2][jj]];

			s->NN2(4, nf, xx, yy, zz) = s->NN2(4, nf, xx, yy, zz) - 1;
			s->NN2(2, nf, xx, yy, zz) = s->NN2(2, nf, xx, yy, zz) + 1;

		}



		for (int ii = 0; ii < 8; ii++) {
			int xx = s->modlx[vacx + v1nbr_bcc[vacnm][0][ii]];
			int yy = s->modly[vacy + v1nbr_bcc[vacnm][1][ii]];
			int zz = s->modlz[vacz + v1nbr_bcc[vacnm][2][ii]];

			s->NN1(4, nf, xx, yy, zz) = s->NN1(4, nf, xx, yy, zz) + 1;
			s->NN1(2, nf, xx, yy, zz) = s->NN1(2, nf, xx, yy, zz) - 1;

		}


		for (int jj = 0; jj < 6; jj++) {
			int xx = s->modlx[vacx + v2nbr_bcc[0][jj]];
			int	yy = s->modly[vacy + v2nbr_bcc[1][jj]];
			int	zz = s->modlz[vacz + v2nbr_bcc[2][jj]];

			s->NN2(4, vacnm, xx, yy, zz) = s->NN2(4, vacnm, xx, yy, zz) + 1;
			s->NN2(2, vacnm, xx, yy, zz) = s->NN2(2, vacnm, xx, yy, zz) - 1;
		}
	}
}
//
//
////判断要交换的点是否越界，不越界，交换 
////越界 把vac接入到越界的vaclists，然后重新选择一个空位，如果没有空位，此进程结束模拟。
//
//
//
void ParaBccJump::ExchangeSites(int& vacx, int& vacy, int& vacz, int& vacnm, int& nf,
	blitz::Array<int, 1> x, blitz::Array<int, 1> y,
	blitz::Array<int, 1> z, blitz::Array<double, 1> w, int& path)
{
	s->Site(vacnm, vacx, vacy, vacz) = s->Site(nf, x(path), y(path), z(path));
	vacx = x(path); vacy = y(path); vacz = z(path); vacnm = nf;
	s->Site(vacnm, vacx, vacy,vacz) = 2;

}

//
////检查是否越界
bool ParaBccJump::Border(int& vacx, int& vacy, int& vacz, int& vacnm,
	int i, int j, int k,int idx,int idy,int idz) {

	if(vacx >= i * idx && vacx < (i + 1) * idx && vacy >= j * idy && vacy < (j + 1) * idy && vacz >= k * idz && vacz < (k + 1) * idz){
		return false;
	}
	else {
		return true;
	}
}

// 重新赋值空位


