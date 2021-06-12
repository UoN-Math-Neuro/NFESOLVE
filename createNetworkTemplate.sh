#!/bin/bash

SimName=$1;
mkdir ${SimName};
Makefile=${SimName}"/Makefile";

printf "PROGRAM_NAME := "${SimName}"\nDEPEND := "${SimName}"_Driver.cpp "${SimName}".cpp\n" >> ${Makefile};
printf 'NFESOLVE_DIR := $(HOME)/NFESOLVE-master/\nNFESOLVE_INCDIR := $(NFESOLVE_DIR)/include\nNFESOLVE_LIBDIR := $(NFESOLVE_DIR)/lib\n' >> ${Makefile};
printf 'CC := g++\nCFLAGS := -std=c++11 -O2 -static-libgcc -static-libstdc++\nOMP := -fopenmp\nINC := -I $(NFESOLVE_INCDIR)\n' >> ${Makefile};
printf "LIB_INCLUDES := -larmadillo\n" >> ${Makefile};
printf 'all: $(PROGRAM_NAME)\n$(PROGRAM_NAME) : $(DEPEND)\n\t' >> ${Makefile};
printf '$(CC) $(CFLAGS) $(INC) -o $(PROGRAM_NAME) $^ -L$(NFESOLVE_LIBDIR) -lNFESOLVE $(LIB_INCLUDES)\n' >> ${Makefile};
printf 'clean :\n\trm -f $(PROGRAM_NAME)' >> ${Makefile};

SimCPP=${SimName}/${SimName}.cpp;
printf '#include "'${SimName}'.hpp"\n' >> ${SimCPP};
printf ${SimName}"::"$SimName"(const arma::vec& parameters, arma::mat& connectivity)\n{\n" >> ${SimCPP};
printf "mParameters = parameters;\nmpConnectivity = &connectivity;\n}\n" >> ${SimCPP};
printf "void "${SimName}"::ComputeF(const double t, const arma::vec& u, arma::vec& F) const\n{\nF=-u;\n}" >> ${SimCPP};

SimHPP=${SimName}/${SimName}.hpp;
printf '#ifndef '${SimName}'HEADERDEF\n#define '${SimName}'HEADERDEF\n#include "NFESOLVE.hpp"\n#include <armadillo>\nclass '${SimName}" : public ODEInterface\n{\n" >> ${SimHPP};
printf "public:\n"${SimName}"(const arma::vec& parameters, arma::mat& connectivity);\n" >> ${SimHPP};
printf "void ComputeF(const double t, const arma::vec& u, arma::vec& F) const;\nprivate:\narma::vec mParameters;\n" >> ${SimHPP};
printf 'arma::mat* mpConnectivity;\n};\n#endif' >> ${SimHPP};

Driver=${SimName}/${SimName}_Driver.cpp;
printf '#include "NFESOLVE.hpp"\n#include "'${SimName}'.hpp"\n#include <armadillo>\n#include <iostream>\n#include <string>\n' >> ${Driver};
printf "int main(int argc, char* argv[])\n{\n" >> ${Driver};
printf 'double initialTime = 0.0;\ndouble finalTime = 100;\ndouble stepSize = 0.01;\nstd::string outputFileName = "'${SimName}'.dat";\n' >> ${Driver};
printf "int saveGap = 1;\nint printGap = 1;\n" >> ${Driver};
printf "double tau = 1.0;\narma::vec P(1); P.fill(tau);\n" >> ${Driver};
printf 'arma::mat W;\nW.load("C.txt");\narma::uword N = W.n_rows;\n' >> ${Driver};
printf "arma::vec IC(2); IC(0) = 1.0; IC(1) = -1.0;\n"$SimName" prob(P, W);\n" >> ${Driver};
printf "RungeKutta4Solver solver(prob, IC, initialTime, finalTime, stepSize, outputFileName, saveGap, printGap);\nsolver.Solve();\nreturn 0;\n}" >> ${Driver};

printf "1, 0\n0, 1" >> ${SimName}"/C.txt"