#pragma once
#include <iostream>
#include <string>
#include <fstream>
#include <cassert>
#include <vector>
#include <utility>
#include "system.h"


class csvHandler {
public:
	csvHandler(std::string filename, std::string comment) {
		assert(filename != "");
		assert(comment.at(0) == '#'  && "Comment format is wrong- every line starts with hash");
		m_fileName = filename;
		m_path = generatePath(m_fileName);
		m_comment = comment;
	}
	void WriteMatrixToFile(int** matrix, int NrSamplingLengths) {
		std::ofstream OutFile;
		OutFile.open(m_path, std::fstream::trunc);
		for (int i = 0; i < NrSamplingLengths; i++) {
			OutFile << matrix[i][0];
			for (int j = 1; j < NrSamplingLengths; j++) {
				OutFile << "," << matrix[i][j];
			}
			OutFile << "\n";
		}
		OutFile.close();
	}

	void WriteToCSV(std::vector<std::pair<std::string, std::vector<double>>> &dataStruct) {
		std::ofstream OutFile;
		OutFile.open(m_path, std::fstream::trunc);

		//adds the comment to the top of the csv file
		OutFile << m_comment << "\n";

		int largestColumnSize = 0;
		std::vector<int> columnsizes(dataStruct.size());
		// adds the header to the csv file and finds the length of the longest column
		for (int i = 0; i < dataStruct.size(); i++) {
			OutFile << dataStruct.at(i).first;
			if (i != dataStruct.size() - 1) {
				OutFile << ",";
			}
			columnsizes[i] = dataStruct.at(i).second.size();
			if (columnsizes[i] > largestColumnSize) {
				largestColumnSize = columnsizes[i];
			}
		}
		OutFile << "\n";

		// add the data structure to the csv file 
		double item;
		for (int i = 0; i < largestColumnSize; i++) {
			for (int j = 0; j < dataStruct.size(); j++) {
				if (i < columnsizes[j]) {
					item = dataStruct.at(j).second.at(i);
					OutFile << item;
				}
				if (j != dataStruct.size() - 1) {
					OutFile << ",";
				}
			}
			OutFile << "\n";
		}
		OutFile.close();
	}
void WriteToCSV(std::vector<std::pair<std::string, std::vector<int>>>& dataStruct) {
	std::ofstream OutFile;
	OutFile.open(m_path, std::fstream::trunc);

	//adds the comment to the top of the csv file
	OutFile << m_comment << "\n";

	int largestColumnSize = 0;
	std::vector<int> columnsizes(dataStruct.size());
	// adds the header to the csv file and finds the length of the longest column
	for (int i = 0; i < dataStruct.size(); i++) {
		OutFile << dataStruct.at(i).first;
		if (i != dataStruct.size() - 1) {
			OutFile << ",";
		}
		columnsizes[i] = dataStruct.at(i).second.size();
		if (columnsizes[i] > largestColumnSize) {
			largestColumnSize = columnsizes[i];
		}
	}
	OutFile << "\n";

	// add the data structure to the csv file 
	double item;
	for (int i = 0; i < largestColumnSize; i++) {
		for (int j = 0; j < dataStruct.size(); j++) {
			if (i < columnsizes[j]) {
				item = dataStruct.at(j).second.at(i);
				OutFile << item;
			}
			if (j != dataStruct.size() - 1) {
				OutFile << ",";
			}
		}
		OutFile << "\n";
	}
	OutFile.close();
}


private:
	std::string generatePath(std::string fileName) {
		return "..//..//..//Output//" + fileName + ".csv";
	}

	std::string m_fileName = "";
	std::string m_path = "";
	std::string m_comment = "";
};