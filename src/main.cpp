#include <chrono> 

#include "kmc_input.h"
#include "kmc_init.h"
#include "kmc_jump.h"
#include "kmc_cluster.h"
#include <ctime>


std::string getInputFileName(int argc, const char* argv[]);
std::unique_ptr<Base> createLattice(SimulationParameters& parameter);
std::unique_ptr<JumpBase> createJump(std::unique_ptr<Base>& lattice);
std::unique_ptr<CrystalMethod> createMethod(SimulationParameters& parameter);
//void runSimulationLoop(SimulationParameters& parameter,std::unique_ptr<Base>& lattice, std::unique_ptr<JumpBase>& jump, std::unique_ptr<CrystalMethod>& method);

int main(int argc, const char* argv[]) {

	time_t t1;
	time_t t2;
	struct tm tm1, tm2;
	double seconds;
	//const std::string input_file_name = getInputFileName(argc, argv);
	
	time(&t1);//获取现在的时间

	SimulationParameters parameter("input.yaml");
	
	auto lattice = createLattice(parameter);

	auto method = createMethod(parameter);

	//method->Print(lattice, 1, ".");

	auto jump = createJump(lattice);

	int num = 0;
	for (int i = 0; i < parameter.par_step; i++) {

		jump->runSimulationLoop();

	
		if (i % parameter.outposcar == 0) {
			method->Print(lattice,num,parameter.filepackage);
			num++;
		}
	}
	//method->OutSize(parameter.filepackage);

	//auto method = createMethod(parameter);
	
	//runSimulationLoop(parameter,lattice, jump, method);
	time(&t2);//获取现在的时间
	seconds = difftime(t2, t1);//返回double类型
	std::cout << "difftime(): " << seconds << std::endl;

	return 0;
}



std::string getInputFileName(int argc, const char* argv[]) {
	if (argc < 2) {
		//std::cerr << "用法：" << argv[0] << " <输入文件名>";
		exit(1);
	}
	return std::string(argv[1]);
}


std::unique_ptr<Base> createLattice(SimulationParameters& parameter) {
	return TypeBuilder::create(parameter);
}

std::unique_ptr<JumpBase> createJump(std::unique_ptr<Base>& lattice) {
	return JumpBuilder::create(lattice);
}

std::unique_ptr<CrystalMethod> createMethod(SimulationParameters& parameter) {
	return MethodBuilder::create(parameter);
}

/*void runSimulationLoop(SimulationParameters& parameter,std::unique_ptr<Base>& lattice, std::unique_ptr<JumpBase>& jump, std::unique_ptr<CrystalMethod>& method)*/
//	for (int i = 0; i < 1000000; i++) {
//		//jump->CalculatedEnergy();
//		//jump->SumAndChoice();
//		jump->UpdateNeighbor();
//		jump->ExchangeSite();
//		if (i % parameter.step_log == 0) {
//			method->Cluster(lattice);
//		}
//		/*if (i % parameter.outposcar == 0) {
//			method->Print(lattice, i / parameter.outposcar, parameter.filepackage);
//			method->Pairs(lattice);
//		}*/
//	}
//	//method->OutSize(parameter.filepackage);
//}