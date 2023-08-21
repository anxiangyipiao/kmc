#include <chrono> 

#include "kmc_input.h"
#include "kmc_init.h"
#include "kmc_jump.h"
#include "kmc_cluster.h"

std::string getInputFileName(int argc, const char* argv[]);
std::unique_ptr<Base> createLattice(SimulationParameters& parameter);
std::unique_ptr<JumpBase> createJump(std::unique_ptr<Base>& lattice);
std::unique_ptr<CrystalMethod> createMethod(SimulationParameters& parameter);
//void runSimulationLoop(SimulationParameters& parameter,std::unique_ptr<Base>& lattice, std::unique_ptr<JumpBase>& jump, std::unique_ptr<CrystalMethod>& method);

int main(int argc, const char* argv[]) {

	
	const std::string input_file_name = getInputFileName(argc, argv);
	
	SimulationParameters parameter(input_file_name);

	auto lattice = createLattice(parameter);

	auto method = createMethod(parameter);

	// method->Print(lattice, init, parameter.filepackage);

	auto jump = createJump(lattice);

	for (int i = 0; i < parameter.par_step; i++) {

		jump->runSimulationLoop();

		// if (i % parameter.step_log == 0) {
		// 	method->Cluster(lattice);
		// }
		if (i % parameter.outposcar == 0) {
			method->Print(lattice,i,parameter.filepackage);
		}
	}
	

	//auto method = createMethod(parameter);
	
	//runSimulationLoop(parameter,lattice, jump, method);
	

	return 0;
}



std::string getInputFileName(int argc, const char* argv[]) {
	if (argc < 2) {
		std::cerr << "�÷���" << argv[0] << " <�����ļ���>\n";
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
