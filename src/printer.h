#ifndef _PRINTER_HPP_
#define _PRINTER_HPP_

#include<fstream>
#include<string>
#include"fdaPDE.h"
#include <iomanip>

class printer{


// std::string directory = "/Users/giuliopn/PACSworkspace3/PACSworkspace/GAM_tests/debugging_output/CPP/";
// std::string directory = "/Users/giuliopn/PACSworkspace4/PACS_fdaPDE/GAM_tests/debugging_output/CPP/";

// std::string directory = "/home/alb/Scrivania/PACS/Git_Folder/debugging_output/CPP/";
public:
	static int counter;
static void SaveMatrixXr(std::string& name_txt, const MatrixXr & mat){
		const UInt numRows = mat.rows();
		const UInt numCols = mat.cols();
		constexpr UInt CONSTPREC = 16;


		std::string directory = "/Users/giuliopn/PACSworkspace5/debugging_output/CPP/";
		name_txt = directory + name_txt;

		std::ofstream FileDataMatrix(name_txt);
		if(FileDataMatrix.is_open()){
			for(int i=0; i<numRows; i++){
				for(int j=0; j<numCols; j++){
					FileDataMatrix<< std::setprecision(20) << mat(i,j);
					if( j< numCols -1 )FileDataMatrix<< ',';
				}
			FileDataMatrix<<'\n';
			}
		}
		FileDataMatrix.close();
}

template<typename T>
static void SaveDimension(std::string& name_txt, const std::vector<T> & vec){
		const UInt numelem = vec.size();
		constexpr UInt CONSTPREC = 16;


		std::string directory = "/Users/giuliopn/PACSworkspace5/debugging_output/CPP/";
		name_txt = directory + name_txt;

		std::ofstream FileDataMatrix(name_txt);
		if(FileDataMatrix.is_open()){
				FileDataMatrix<< vec.size();
		}
		FileDataMatrix.close();
}


template<typename T>
static void STDvector(std::string& name_txt, const std::vector<T> & vec){
		const UInt numelem = vec.size();
		constexpr UInt CONSTPREC = 16;

		std::string directory = "/Users/giuliopn/PACSworkspace5/debugging_output/CPP/";
		name_txt = directory + name_txt;

		std::ofstream File(name_txt);
		if(File.is_open()){
			for(int k=0; k<numelem; k++){
					File<< std::scientific <<vec[k]<<'\n';
			}
		}
		File.close();
}


static void saveVectorXr(std::string& name_txt, const VectorXr & vect){
   	const UInt size = vect.size();
		constexpr UInt CONSTPREC = 16;

	    std::string directory =  "/Users/giuliopn/PACSworkspace5/debugging_output/CPP/";
		name_txt = directory + name_txt;

	std::ofstream File(name_txt);
	if(File.is_open()){
		for(int k=0; k<size; k++){
				File<< std::scientific <<vect[k]<<'\n';
		}
	}
	File.close();
}

static void milestone(std::string name_txt){
	counter++;
	std::string directory =  "/Users/giuliopn/PACSworkspace5/debugging_output/CPP/";
	name_txt = directory + name_txt;
	std::ofstream File(name_txt);
	File<<counter;
	File.close();
}

static void variableInt(std::string name_txt, int x){
	std::string directory =  "/Users/giuliopn/PACSworkspace5/debugging_output/CPP/";
	name_txt = directory + name_txt;
	std::ofstream File(name_txt);
	File<<x;
	File.close();
}

};

int printer::counter = 0;

#endif
