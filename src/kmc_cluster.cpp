#include "kmc_cluster.h"


vector<double> BCCmethod::ToCoodinate(int a, int b, int c, int d) {
	
	vector<double> tempvec;
	if (d == 0) {
		tempvec.push_back(a); tempvec.push_back(b); tempvec.push_back(c);
		return tempvec;
	}
	else {
		tempvec.push_back(a + 0.5); tempvec.push_back(b + 0.5); tempvec.push_back(c + 0.5);
		return tempvec;
	}
}


void BCCmethod::Cluster(std::unique_ptr<Base>& s) {

	unordered_map<int, unordered_map<int, unordered_map<int, unordered_map<int, int>>>> Bsite;
	vector<int> Num;

	//find Blist
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < this->ny; j++) {
			for (int k = 0; k < this->nz; k++) {
				for (int n = 0; n < 2; n++) {
					if (s->Site(n,i,j,k) == 1) {
						Bsite[i][j][k][n] = 1;
					}
				}
			}
		}
	}

	//多维map查找
	unordered_map<int, unordered_map<int, unordered_map<int, unordered_map<int, int>>>>::iterator p1;
	unordered_map<int, unordered_map<int, unordered_map<int, int>>>::iterator p2;
	unordered_map<int, unordered_map<int, int>>::iterator p3;
	unordered_map<int, int>::iterator p4;

	for (p1 = Bsite.begin(); p1 != Bsite.end(); p1++)
	{
		for (p2 = p1->second.begin(); p2 != p1->second.end(); p2++)
		{
			for (p3 = p2->second.begin(); p3 != p2->second.end(); p3++)
			{
				for (p4 = p3->second.begin(); p4 != p3->second.end(); p4++)
				{
					//if (Bsite[p1->first][p2->first][p3->first][p4->first] == 1) {
					if (p4->second == 1) {
						int nmn; int cluster = 1; p4->second = 0;
						if (p4->first == 0) {
							nmn = 1;
						}else{
							nmn = 0;
						}

						for (int inv = 0; inv < 8; inv++) {

							int x = (s->modlx[p1->first + v1nbr_bcc[p4->first][0][inv]]);
							int y = (s->modly[p2->first + v1nbr_bcc[p4->first][1][inv]]);
							int z = (s->modlz[p3->first + v1nbr_bcc[p4->first][2][inv]]);

							if (Bsite.find(x) != Bsite.end()) {
								if (Bsite[x].find(y) != Bsite[x].end()) {
									if (Bsite[x][y].find(z) != Bsite[x][y].end()) {
										if (Bsite[x][y][z].find(nmn) != Bsite[x][y][z].end()) {
											if (Bsite[x][y][z][nmn] == 1) {

												cluster = cluster + 1;
												cluster = Deep(cluster, x, y, z, nmn, s, Bsite);

											}
										}
									}
								}
							}
						}
						Num.emplace_back(cluster);					
					}
				}
			}
		}
	}

	this->Size(Num);
	this->Pairs(s);
	//CrystalMethod<T>::Size(Num);

}


int  BCCmethod::Deep(int cluster,int  a, int b,int  c, int nmn, std::unique_ptr<Base>& s,unordered_map<int, unordered_map<int, unordered_map<int, unordered_map<int, int>>>> &Bsite){

	Bsite[a][b][c][nmn] = 0;// 这里查询过了，下次跳过

	int nms;
	switch (nmn) {
	case 0:
		nms = 1;
	case 1:
		nms = 0;
	}

	for (int inv = 0; inv < 8; inv++) {

		int x = (s->modlx[a + v1nbr_bcc[nmn][0][inv]]);
		int y = (s->modly[b + v1nbr_bcc[nmn][1][inv]]);
		int z = (s->modlz[c + v1nbr_bcc[nmn][2][inv]]);

		if (Bsite.find(x) != Bsite.end()) {
			if (Bsite[x].find(y) != Bsite[x].end()) {
				if (Bsite[x][y].find(z) != Bsite[x][y].end()) {
					if (Bsite[x][y][z].find(nms) != Bsite[x][y][z].end()) {
						if (Bsite[x][y][z][nms] == 1) {

							cluster = cluster + 1;
							cluster = Deep(cluster, x, y, z, nms, s, Bsite);

						}
					}
				}
			}
		}
	}

	return cluster;
}


void CrystalMethod::Size(vector<int>& Num) {

	int alone = count(Num.begin(), Num.end(), 1);
	cout << alone << endl;

	map<int, int>tempmap;
	for (int i = 0; i < Num.size(); i++) {
		//存在Key
		if (tempmap.count(Num[i]) > 0) {
			tempmap[Num[i]] = tempmap[Num[i]] + 1;
		}
		else {
			tempmap[Num[i]] = 1;
		}
	}

	vector<double> Allradius;
	int number = 0;
	map<int, int>::iterator it;
	for (it = tempmap.begin(); it != tempmap.end(); it++)
	{
		if(it->first > radius_start){
			//半径
			double r = (cbrt(it->first) * 2.78 * sqrt(3)) / 40;	//per_radius
			for(int i= 0; i < it->second; i++){
				Allradius.push_back(r);
			}
			//数密度
			number = number + it->second;
		}
	}
	//倒叙排序
	reverse(Allradius.begin(), Allradius.end());
	if (!Allradius.empty()) {
		//最大半径
		MaxRadius.push_back(Allradius.front());
		//cout << Allradius.front() <<endl;
		//平均半径
		double sum = 0;
		for (int i = 0; i < Allradius.size(); i++) {
			sum += Allradius[i];
		}
		
		AverageRadius.push_back(sum/ Allradius.size());
		//数密度
		NumbericalDensity.push_back(number / (nx * ny * nz));
	}
	else {
		//最大半径
		MaxRadius.push_back(0);
		//cout << Allradius[0]<<endl;
		//平均半径
		AverageRadius.push_back(0);
		//数密度
		NumbericalDensity.push_back(0);
	}
	//孤立原子
	Alone.push_back(alone);
	//原子分布
	DictRadius.push_back(tempmap);
}


void CrystalMethod::OutSize(string path)
{
	// 打开输出文件流
	ofstream outFile1(path+"/output1.dat");
	outFile1 << "Alone" << " " << "AverageRadius" << " " << "MaxRadius" << " " << "NumbericalDensity" << " " <<
		"A-B" << " " << "A-V" << " " << "A-A" << " " << "B-B" << " " << "B-V" << " " << "V-V" << endl;
	// 写入变量
	for (int i = 0; i < Alone.size(); i++) {
		outFile1 << Alone[i] << " "<<AverageRadius[i] << " " << MaxRadius[i] << " "<< NumbericalDensity[i] <<" "<<TPair[i]["01"] << " " << TPair[i]["02"] << " " << TPair[i]["00"] << " " << TPair[i]["11"] << " " << TPair[i]["12"] << " " << TPair[i]["22"] <<endl;
	}

	outFile1.close();

	ofstream outFile2(path+"/output2.dat");
	for (int i = 0; i < Alone.size(); i++) {
		for (auto it = DictRadius[i].begin(); it != DictRadius[i].end(); it++) {
			outFile2 << it->first << " ";
		}
		outFile2 << "value";

		for (auto it = DictRadius[i].begin(); it != DictRadius[i].end(); it++) {
			outFile2 << it->second << " ";
		}
		outFile2 << endl;
	}

	// 关闭文件流
	outFile2.close();

}


void BCCmethod::Pairs(std::unique_ptr<Base>& s) {
	
	//深拷贝
	blitz::Array<int, 4> arr_copy(s->Site.shape());
	arr_copy = s->Site;

	map<string, int> Pair;
	Pair["01"] = 0, Pair["02"], Pair["00"], Pair["11"], Pair["12"], Pair["22"] = 0, 0, 0, 0, 0, 0;
	int nmn = 0;
	for (int nm = 0; nm < 2; nm++) {
		if (nm == 0) {
			nmn = 1;
		}
		else {
			nmn = 0;
		}

		for (int k = 0; k < this->nz; k++) {
			for (int j = 0; j < this->ny; j++) {
				for (int i = 0; i < this->nx; i++) {
					for (int inv = 0; inv < 8; inv++) {

						int ii = s->modlx[i + v1nbr_bcc[nm][0][inv]];
						int jj = s->modly[j + v1nbr_bcc[nm][1][inv]];
						int kk = s->modlz[k + v1nbr_bcc[nm][2][inv]];

						switch (arr_copy(nm,i,j,k)){
						case 0:
							switch (arr_copy(nmn,ii,jj,kk)) {
							case 0:
								Pair["00"] = Pair["00"] + 1;
							case 1:
								Pair["01"] = Pair["01"] + 1;
							case 2:
								Pair["02"] = Pair["02"] + 1;
							case 3:
								continue;
							}

						case 1:
							switch (arr_copy(nmn, ii, jj, kk)) {
							case 0:
								Pair["01"] = Pair["01"] + 1;
							case 1:
								Pair["11"] = Pair["11"] + 1;
							case 2:
								Pair["12"] = Pair["12"] + 1;
							case 3:
								continue;
							}
						case 2:
							switch (arr_copy(nmn, ii, jj, kk)) {
							case 0:
								Pair["02"] = Pair["02"] + 1;
							case 1:
								Pair["12"] = Pair["12"] + 1;
							case 2:
								Pair["22"] = Pair["22"] + 1;
							case 3:
								continue;
							}
						}

					}
					arr_copy(nm, i, j, k) = 3;  // pass
					}
				}
			}
		}
	
	this->TPair.emplace_back(Pair);
}


void BCCmethod::Print(std::unique_ptr<Base>& obj,int a,std::string filepackage) {
		
	// Create directory if it does not exist
		
		ofstream outFile;	//定义ofstream对象outFile
	
	    outFile.open(filepackage +"/POSCAR"+std::to_string(a), ios::out);	//打开文件
	
	    outFile << "Fe" << endl;
	    outFile << "1.0000000000000000" << endl;
	    outFile << this->nx << "\t" << 0.0000000000000000 << "\t" << 0.0000000000000000 << endl;
	    outFile << 0.0000000000000000 << "\t" << this->nx << "\t" << 0.0000000000000000 << endl;
	    outFile << 0.0000000000000000 << "\t" << 0.0000000000000000 << "\t" << this->nx << endl;
	    outFile << "Fe" <<"\t"<< "Cu" << endl;
	    outFile << obj->Anum << "\t" << obj->Bnum << endl;
	    outFile << "Cartesian"<< endl;
	
	    for (int i = 0; i < this->nx; i++) {
	        for (int j = 0; j < this->ny; j++) {
	            for (int k = 0; k < this->nz; k++) {
	                for (int n = 0; n < 2; n++) {
	                    if (obj->Site(n,i,j,k) == 0) {
	                        vector<double> tempvec = ToCoodinate(i, j, k, n);
	                        outFile << tempvec[0] << "\t" << tempvec[1] << "\t" << tempvec[2] << endl; //写入操作
	                    }                     
	                }
	            }
	        }
	    }
	
	    for (int i = 0; i < this->nx; i++) {
	        for (int j = 0; j < this->ny; j++) {
	            for (int k = 0; k < this->nz; k++) {
	                for (int n = 0; n < 2; n++) {
	                    if (obj->Site(n, i, j, k) == 1) {
	                        vector<double> tempvec = ToCoodinate(i, j, k, n);
	                        outFile << tempvec[0] << "\t" << tempvec[1] << "\t" << tempvec[2] << endl; //写入操作
	                    }
	                }
	            }
	        }
	    }
	
	    outFile.close();	//关闭文件
}



//void OutMethods::MSD(unordered_map<int, unordered_map<int, unordered_map<int, unordered_map<int, int>>>> Site)
//{
//	
//	//msd
//	int Fe = 0; int Cu = 0; int V = 0;
//	for (int i = 0; i < nx; i++) {
//		for (int j = 0; j < ny; j++) {
//			for (int k = 0; k < nz; k++) {
//				for (int n = 0; n < Length; n++) {
//
//					if (Site[i][j][k][n] == 0) {
//						vector<double> Ord1 = OriginIndex[PastIndex[i][j][k][n]];
//						vector<double> Ord2 = ToCoodinate(i, j, k, n);
//						Fe = CalculateMsd(Ord1, Ord2) + Fe;
//
//					}
//					else if (Site[i][j][k][n] == 1) {
//						vector<double> Ord1 = OriginIndex[PastIndex[i][j][k][n]];
//						vector<double> Ord2 = ToCoodinate(i, j, k, n);
//						Cu = CalculateMsd(Ord1, Ord2) + Cu;
//					}
//					else {
//						vector<double> Ord1 = OriginIndex[PastIndex[i][j][k][n]];
//						vector<double> Ord2 = ToCoodinate(i, j, k, n);
//						V = CalculateMsd(Ord1, Ord2) + V;
//					}
//
//				}
//			}
//		}
//	}
//	Felist.push_back(Fe); Culist.push_back(Cu); Vlist.push_back(V);
//	/*cout << Fe << endl;
//	cout << Cu << endl;
//	cout << V << endl;*/
//
//}


//double OutMethods::CalculateMsd(vector<double> Num1, vector<double> Num2) {
//
//	double val = sqrt(pow(2.87 * (Num1[0] - Num2[0]), 2) + pow(2.87 * (Num1[1] - Num2[1]), 2) + pow(2.87 * (Num1[2] - Num2[2]), 2));
//	return val;
//}



vector<double> FCCmethod::ToCoodinate(int a, int b, int c, int n)
{

	vector<double> tempvec;

	if (n == 0) {
		tempvec.push_back(a); tempvec.push_back(b); tempvec.push_back(c);
		return tempvec;;
	}
	else if (n == 1) {
		tempvec.push_back(a + 0.5); tempvec.push_back(b + 0.5); tempvec.push_back(c);
		return tempvec;
	}
	else if (n == 2) {
		tempvec.push_back(a + 0.5); tempvec.push_back(b); tempvec.push_back(c + 0.5);
		return tempvec;
	}
	else {
		tempvec.push_back(a); tempvec.push_back(b + 0.5); tempvec.push_back(c + 0.5);
		return tempvec;
	}


}


void FCCmethod::Cluster(std::unique_ptr<Base>& s)
{
	unordered_map<int, unordered_map<int, unordered_map<int, unordered_map<int, int>>>> Bsite;
	vector<int> Num;

	//find Blist
	for (int i = 0; i < this->nx; i++) {
		for (int j = 0; j < this->ny; j++) {
			for (int k = 0; k < this->nz; k++) {
				for (int n = 0; n < 4; n++) {
					if (s->Site(n, i, j, k) == 1) {
						Bsite[i][j][k][n] = 1;
					}
				}
			}
		}
	}

	//多维map查找
	unordered_map<int, unordered_map<int, unordered_map<int, unordered_map<int, int>>>>::iterator p1;
	unordered_map<int, unordered_map<int, unordered_map<int, int>>>::iterator p2;
	unordered_map<int, unordered_map<int, int>>::iterator p3;
	unordered_map<int, int>::iterator p4;

	for (p1 = Bsite.begin(); p1 != Bsite.end(); p1++)
	{
		for (p2 = p1->second.begin(); p2 != p1->second.end(); p2++)
		{
			for (p3 = p2->second.begin(); p3 != p2->second.end(); p3++)
			{
				for (p4 = p3->second.begin(); p4 != p3->second.end(); p4++)
				{
					//if (Bsite[p1->first][p2->first][p3->first][p4->first] == 1) {
					if (p4->second == 1) {
						int cluster = 1; p4->second = 0;

						for (int nmn = 0; nmn < 4; nmn++) {
							if (nmn != p4->first) {
								for (int inv = 0; inv < 4; inv++) {

									int x = s->modlx[p1->first + v1nbr_fcc[0][p4->first][nmn][inv]];
									int y = s->modly[p2->first + v1nbr_fcc[1][p4->first][nmn][inv]];
									int z = s->modlz[p3->first + v1nbr_fcc[2][p4->first][nmn][inv]];

									if (Bsite.find(x) != Bsite.end()) {
										if (Bsite[x].find(y) != Bsite[x].end()) {
											if (Bsite[x][y].find(z) != Bsite[x][y].end()) {
												if (Bsite[x][y][z].find(nmn) != Bsite[x][y][z].end()) {
													if (Bsite[x][y][z][nmn] == 1) {

														cluster = cluster + 1;
														cluster = Deep(cluster, x, y, z, nmn, s, Bsite);

													}
												}
											}
										}
									}
								}
							}
						}

						Num.emplace_back(cluster);
					}
				}
			}
		}
	}

	this->Size(Num);
}

int  FCCmethod::Deep(int cluster, int  x, int y, int  z, int nm, std::unique_ptr<Base>& s, unordered_map<int, unordered_map<int, unordered_map<int, unordered_map<int, int>>>>& Bsite) {

	Bsite[x][y][z][nm] = 0;// 这里查询过了，下次跳过


	for (int nmn = 0; nmn < 4; nmn++) {
			if (nmn != nm) {
				for (int inv = 0; inv < 4; inv++) {

					int xx = s->modlx[x + v1nbr_fcc[0][nm][nmn][inv]];
					int yy = s->modly[y + v1nbr_fcc[1][nm][nmn][inv]];
					int zz = s->modlz[z + v1nbr_fcc[2][nm][nmn][inv]];

					if (Bsite.find(xx) != Bsite.end()) {
						if (Bsite[xx].find(yy) != Bsite[xx].end()) {
							if (Bsite[xx][yy].find(zz) != Bsite[xx][yy].end()) {
								if (Bsite[xx][yy][zz].find(nmn) != Bsite[xx][yy][zz].end()) {
									if (Bsite[xx][yy][zz][nmn] == 1) {

										cluster = cluster + 1;
										cluster = Deep(cluster, xx, yy, zz, nmn, s, Bsite);
									}
								}
							}
						}
					}
				}
			}
		
	}

	return cluster;
}

void FCCmethod::Print(std::unique_ptr<Base>& obj,int a,std::string filename)
{
	ofstream outFile;	//定义ofstream对象outFile

	outFile.open(filename + "/POSCAR" + std::to_string(a), ios::out);	//打开文件

	outFile << "Fe" << endl;
	outFile << "1.0000000000000000" << endl;
	outFile << this->nx << "\t" << 0.0000000000000000 << "\t" << 0.0000000000000000 << endl;
	outFile << 0.0000000000000000 << "\t" << this->nx << "\t" << 0.0000000000000000 << endl;
	outFile << 0.0000000000000000 << "\t" << 0.0000000000000000 << "\t" << this->nx << endl;
	outFile << "Fe" << "\t" << "Cu" << endl;
	outFile << obj->Anum << "\t" << obj->Bnum << endl;
	outFile << "Cartesian" << endl;

	for (int i = 0; i < this->nx; i++) {
		for (int j = 0; j < this->ny; j++) {
			for (int k = 0; k < this->nz; k++) {
				for (int n = 0; n < 4; n++) {
					if (obj->Site(n, i, j, k) == 0) {
						vector<double> tempvec = ToCoodinate(i, j, k, n);
						outFile << tempvec[0] << "\t" << tempvec[1] << "\t" << tempvec[2] << endl; //写入操作
					}

				}
			}
		}
	}

	for (int i = 0; i < this->nx; i++) {
		for (int j = 0; j < this->ny; j++) {
			for (int k = 0; k < this->nz; k++) {
				for (int n = 0; n < 2; n++) {
					if (obj->Site(n, i, j, k) == 1) {
						vector<double> tempvec = ToCoodinate(i, j, k, n);
						outFile << tempvec[0] << "\t" << tempvec[1] << "\t" << tempvec[2] << endl; //写入操作
					}
				}
			}
		}
	}

	outFile.close();	//关闭文件

}

void FCCmethod::Pairs(std::unique_ptr<Base>& s)
{
	//深拷贝
	blitz::Array<int, 4> arr_copy(s->Site.shape());
	arr_copy = s->Site;

	map<string, int> Pair;
	Pair["01"] = 0, Pair["02"], Pair["00"], Pair["11"], Pair["12"], Pair["22"] = 0, 0, 0, 0, 0, 0;

	for (int nm = 0; nm < 4; nm++) {
		for (int k = 0; k < this->nz; k++) {
			for (int j = 0; j < this->ny; j++) {
				for (int i = 0; i < this->nx; i++) {
					for (int nmn = 0; nmn < 4; nmn++) {
						if (nmn != nm) {
							for (int inv = 0; inv < 4; inv++) {

								int ii = s->modlx[i + v1nbr_fcc[0][nm][nmn][inv]];
								int jj = s->modly[j + v1nbr_fcc[1][nm][nmn][inv]];
								int kk = s->modlz[k + v1nbr_fcc[2][nm][nmn][inv]];

								switch (arr_copy(nm, i, j, k)) {
								case 0:
									switch (arr_copy(nmn, ii, jj, kk)) {
									case 0:
										Pair["00"] = Pair["00"] + 1;
									case 1:
										Pair["01"] = Pair["01"] + 1;
									case 2:
										Pair["02"] = Pair["02"] + 1;
									case 3:
										continue;
									}

								case 1:
									switch (arr_copy(nmn, ii, jj, kk)) {
									case 0:
										Pair["01"] = Pair["01"] + 1;
									case 1:
										Pair["11"] = Pair["11"] + 1;
									case 2:
										Pair["12"] = Pair["12"] + 1;
									case 3:
										continue;
									}
								case 2:
									switch (arr_copy(nmn, ii, jj, kk)) {
									case 0:
										Pair["02"] = Pair["02"] + 1;
									case 1:
										Pair["12"] = Pair["12"] + 1;
									case 2:
										Pair["22"] = Pair["22"] + 1;
									case 3:
										continue;
									}
								}

							}
							arr_copy(nm, i, j, k) = 3;  // pass
						}
					}
				}
			}

			this->TPair.emplace_back(Pair);
		}

	}


}





void SizeBCCmethod::Cluster(std::unique_ptr<Base>& s)
{
	
	//MultiSizeBcc* s = dynamic_cast<MultiSizeBcc*>(ss.get());

	unordered_map<int, unordered_map<int, unordered_map<int, unordered_map<int, int>>>> Bsite;
	vector<int> Num;

	//find Blist
	for (int i = 0; i < this->nx; i++) {
		for (int j = 0; j < this->ny; j++) {
			for (int k = 0; k < this->nz; k++) {
				for (int n = 0; n < 2; n++) {			
					if (s->Site(n, i, j, k) != 0 ) {
						Bsite[i][j][k][n] = 1;
					}
				}
			}
		}
	}

	//多维map查找
	unordered_map<int, unordered_map<int, unordered_map<int, unordered_map<int, int>>>>::iterator p1;
	unordered_map<int, unordered_map<int, unordered_map<int, int>>>::iterator p2;
	unordered_map<int, unordered_map<int, int>>::iterator p3;
	unordered_map<int, int>::iterator p4;

	for (p1 = Bsite.begin(); p1 != Bsite.end(); p1++)
	{
		for (p2 = p1->second.begin(); p2 != p1->second.end(); p2++)
		{
			for (p3 = p2->second.begin(); p3 != p2->second.end(); p3++)
			{
				for (p4 = p3->second.begin(); p4 != p3->second.end(); p4++)
				{
					//if (Bsite[p1->first][p2->first][p3->first][p4->first] == 1) {
					if (p4->second == 1) {
						int nmn; int cluster = 1; p4->second = 0;
						if (p4->first == 0) {
							nmn = 1;
						}
						else {
							nmn = 0;
						}

						for (int inv = 0; inv < 8; inv++) {

							int x = (s->modlx[p1->first + v1nbr_bcc[p4->first][0][inv]]);
							int y = (s->modly[p2->first + v1nbr_bcc[p4->first][1][inv]]);
							int z = (s->modlz[p3->first + v1nbr_bcc[p4->first][2][inv]]);

							if (Bsite.find(x) != Bsite.end()) {
								if (Bsite[x].find(y) != Bsite[x].end()) {
									if (Bsite[x][y].find(z) != Bsite[x][y].end()) {
										if (Bsite[x][y][z].find(nmn) != Bsite[x][y][z].end()) {
											if (Bsite[x][y][z][nmn] == 1) {

												cluster = cluster + 1;
												cluster = Deep(cluster, x, y, z, nmn, s, Bsite);

											}
										}
									}
								}
							}
						}
						Num.emplace_back(cluster);
					}
				}
			}
		}
	}

	this->Size(Num);
}

void SizeBCCmethod::Pairs(std::unique_ptr<Base>& s)
{
	//深拷贝
	blitz::Array<int, 4> arr_copy(s->Site.shape());
	arr_copy = s->Site;

	map<string, int> Pair;
	Pair["01"] = 0, Pair["02"], Pair["00"], Pair["11"], Pair["12"], Pair["22"] = 0, 0, 0, 0, 0, 0;
	Pair["03"] = 0, Pair["04"], Pair["13"], Pair["14"], Pair["23"], Pair["24"] = 0, 0, 0, 0, 0, 0;
	Pair["33"] = 0, Pair["34"], Pair["44"]= 0, 0, 0;
	
	int nmn = 0;
	for (int nm = 0; nm < 2; nm++) {
		if (nm == 0) {
			nmn = 1;
		}
		else {
			nmn = 0;
		}

		for (int k = 0; k < this->nz; k++) {
			for (int j = 0; j < this->ny; j++) {
				for (int i = 0; i < this->nx; i++) {
					for (int inv = 0; inv < 8; inv++) {

						int ii = s->modlx[i + v1nbr_bcc[nm][0][inv]];
						int jj = s->modly[j + v1nbr_bcc[nm][1][inv]];
						int kk = s->modlz[k + v1nbr_bcc[nm][2][inv]];

						switch (arr_copy(nm, i, j, k)) {
						case 0:
							switch (arr_copy(nmn, ii, jj, kk)) {
							case 0:
								Pair["00"] = Pair["00"] + 1;
							case 1:
								Pair["01"] = Pair["01"] + 1;
							case 2:
								Pair["02"] = Pair["02"] + 1;
							case 3:
								Pair["03"] = Pair["03"] + 1;
							case 4:
								Pair["04"] = Pair["04"] + 1;
							case 5:
								continue;
							}
						case 1:
							switch (arr_copy(nmn, ii, jj, kk)) {
							case 0:
								Pair["01"] = Pair["01"] + 1;
							case 1:
								Pair["11"] = Pair["11"] + 1;
							case 2:
								Pair["12"] = Pair["12"] + 1;
							case 3:
								Pair["13"] = Pair["13"] + 1;
							case 4:
								Pair["14"] = Pair["14"] + 1;
							case 5:
								continue;
							}
						case 2:
							switch (arr_copy(nmn, ii, jj, kk)) {
							case 0:
								Pair["02"] = Pair["02"] + 1;
							case 1:
								Pair["12"] = Pair["12"] + 1;
							case 2:
								Pair["22"] = Pair["22"] + 1;
							case 3:
								Pair["23"] = Pair["23"] + 1;
							case 4:
								Pair["24"] = Pair["24"] + 1;
							case 5:
								continue;
							}
						case 3:
							switch (arr_copy(nmn, ii, jj, kk)) {
							case 0:
								Pair["03"] = Pair["03"] + 1;
							case 1:
								Pair["13"] = Pair["13"] + 1;
							case 2:
								Pair["23"] = Pair["23"] + 1;
							case 3:
								Pair["33"] = Pair["33"] + 1;
							case 4:
								Pair["34"] = Pair["34"] + 1;
							case 5:
								continue;
							}
						case 4:
							switch (arr_copy(nmn, ii, jj, kk)) {
							case 0:
								Pair["04"] = Pair["04"] + 1;
							case 1:
								Pair["14"] = Pair["14"] + 1;
							case 2:
								Pair["24"] = Pair["24"] + 1;
							case 3:
								Pair["34"] = Pair["34"] + 1;
							case 4:
								Pair["44"] = Pair["44"] + 1;
							case 5:
								continue;
							}
						}

					}
					arr_copy(nm, i, j, k) = 5;  // pass
				}
			}
		}
	}

	this->TPair.emplace_back(Pair);
	
}

int SizeBCCmethod::Deep(int cluster, int a, int b, int c, int nmn, std::unique_ptr<Base>& s, unordered_map<int, unordered_map<int, unordered_map<int, unordered_map<int, int>>>>& Bsite)
{
	Bsite[a][b][c][nmn] = 0;// 这里查询过了，下次跳过

	int nms;
	switch (nmn) {
	case 0:
		nms = 1;
	case 1:
		nms = 0;
	}

	for (int inv = 0; inv < 8; inv++) {

		int x = (s->modlx[a + v1nbr_bcc[nmn][0][inv]]);
		int y = (s->modly[b + v1nbr_bcc[nmn][1][inv]]);
		int z = (s->modlz[c + v1nbr_bcc[nmn][2][inv]]);

		if (Bsite.find(x) != Bsite.end()) {
			if (Bsite[x].find(y) != Bsite[x].end()) {
				if (Bsite[x][y].find(z) != Bsite[x][y].end()) {
					if (Bsite[x][y][z].find(nms) != Bsite[x][y][z].end()) {
						if (Bsite[x][y][z][nms] == 1) {

							cluster = cluster + 1;
							cluster = Deep(cluster, x, y, z, nms, s, Bsite);

						}
					}
				}
			}
		}
	}

	return cluster;
}

void SizeBCCmethod::Print(std::unique_ptr<Base>& obj, int a, std::string filepackage)
{

	// Create directory if it does not exist

	ofstream outFile;	//定义ofstream对象outFile

	outFile.open(filepackage + "/"+"POSCAR" + std::to_string(a), ios::out);	//打开文件

	outFile << "Fe" << endl;
	outFile << "1.0000000000000000" << endl;
	outFile << this->nx << "\t" << 0.0000000000000000 << "\t" << 0.0000000000000000 << endl;
	outFile << 0.0000000000000000 << "\t" << this->nx << "\t" << 0.0000000000000000 << endl;
	outFile << 0.0000000000000000 << "\t" << 0.0000000000000000 << "\t" << this->nx << endl;
	outFile << "Fe" << "\t" << "Cu" << "\t" << "Al" << "\t" << "Ag" << endl;
	outFile << obj->Anum << "\t" << obj->Bnum << "\t" << obj->Cnum << "\t" << obj->Dnum << endl;
	outFile << "Cartesian" << endl;

	for (int i = 0; i < this->nx; i++) {
		for (int j = 0; j < this->ny; j++) {
			for (int k = 0; k < this->nz; k++) {
				for (int n = 0; n < 2; n++) {
					if (obj->Site(n, i, j, k) == 0) {
						vector<double> tempvec = ToCoodinate(i, j, k, n);
						outFile << tempvec[0] << "\t" << tempvec[1] << "\t" << tempvec[2] << endl; //写入操作
					}
				}
			}
		}
	}

	for (int i = 0; i < this->nx; i++) {
		for (int j = 0; j < this->ny; j++) {
			for (int k = 0; k < this->nz; k++) {
				for (int n = 0; n < 2; n++) {
					if (obj->Site(n, i, j, k) == 1) {
						vector<double> tempvec = ToCoodinate(i, j, k, n);
						outFile << tempvec[0] << "\t" << tempvec[1] << "\t" << tempvec[2] << endl; //写入操作
					}
				}
			}
		}
	}
	for (int i = 0; i < this->nx; i++) {
		for (int j = 0; j < this->ny; j++) {
			for (int k = 0; k < this->nz; k++) {
				for (int n = 0; n < 2; n++) {
					if (obj->Site(n, i, j, k) == 3) {
						vector<double> tempvec = ToCoodinate(i, j, k, n);
						outFile << tempvec[0] << "\t" << tempvec[1] << "\t" << tempvec[2] << endl; //写入操作
					}
				}
			}
		}
	}
	for (int i = 0; i < this->nx; i++) {
		for (int j = 0; j < this->ny; j++) {
			for (int k = 0; k < this->nz; k++) {
				for (int n = 0; n < 2; n++) {
					if (obj->Site(n, i, j, k) == 4) {
						vector<double> tempvec = ToCoodinate(i, j, k, n);
						outFile << tempvec[0] << "\t" << tempvec[1] << "\t" << tempvec[2] << endl; //写入操作
					}
				}
			}
		}
	}
	outFile.close();	//关闭文件
}

vector<double> SizeBCCmethod::ToCoodinate(int a, int b, int c, int d)
{
	vector<double> tempvec;
	if (d == 0) {
		tempvec.push_back(a); tempvec.push_back(b); tempvec.push_back(c);
		return tempvec;
	}
	else {
		tempvec.push_back(a + 0.5); tempvec.push_back(b + 0.5); tempvec.push_back(c + 0.5);
		return tempvec;
	}
}
