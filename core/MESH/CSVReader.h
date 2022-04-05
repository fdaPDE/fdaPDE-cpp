#ifndef __CSV_READER_H__
#define __CSV_READER_H__

#include <cstddef>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>
#include <iostream>
#include <Eigen/Core>
#include <set>

// a super simple class to parse and read simple .csv files
// (works with output of write.csv() from R)
class CSVReader{

 private:
  // read only stream
  std::ifstream file;
  // split a string in a vector of strings according to a separator
  std::vector<std::string> split(std::string input_, std::string separator);
  // remove char j from string s
  std::string remove(std::string s, const char j);

  std::vector<std::string> columnNames_;

 public:
  CSVReader(const std::string& file_) {
    // opens the file
    file.open(file_);
  }

  // free resources
  ~CSVReader() { file.close(); }

  // parse the file returning a vector for each column
  template <typename T>
    std::unordered_map<std::string, std::vector<T>> parseFile();

  // an utility for converting parsed files into Eigen matrix
  template <typename T>
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> toEigen(std::unordered_map<std::string, std::vector<T>>& input_);
};

std::vector<std::string> CSVReader::split(std::string input_, std::string separator) {
  std::vector<std::string> parsedString;

  // keep string position
  size_t j = 0;
  while(j != std::string::npos){
    j = input_.find(separator);
    // store string in token list
    parsedString.push_back(input_.substr(0, j));
    input_.erase(0, j+separator.length());
  }
  
  return parsedString;
}

std::string CSVReader::remove(std::string s, const char j) {
  std::string result = "";
  for(char c : s){ 
    if (c != j) result += c;
  }
  return result;
}

// parse a file by converting its entries according to type T (commonly int,
// double). Observe that the produced map will contain only entries of the
// specified type
template <typename T>
std::unordered_map<std::string, std::vector<T>> CSVReader::parseFile(){
  std::unordered_map<std::string, std::vector<T>> parsedFile;
  // set column by reading first line
  std::string line;
  getline(file, line); // read line from stream
  std::vector<std::string> columnNames = split(line, ",");
  
  // set header
  for(size_t j = 1; j < columnNames.size(); ++j){
    columnNames_.push_back(remove(columnNames[j], '\"'));
    parsedFile[remove(columnNames[j], '\"')];
  }
  
  // read file until end line by line
  while(getline(file, line)){
    // split CSV line in tokens
    std::vector<std::string> parsedLine = split(line, ",");

    // fill data structure
    for(size_t j = 1; j < columnNames.size(); ++j){
      // convert string to type T
      std::istringstream ss(parsedLine[j]);
      T dataPoint;
      ss >> dataPoint;
      // insert value
      parsedFile[remove(columnNames[j], '\"')].push_back(dataPoint);
    }
  }

  return parsedFile;
}

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> CSVReader::toEigen(std::unordered_map<std::string, std::vector<T>>& input_) {
  // get number of columns
  unsigned int numColumns = input_.size();
  // get number of rows (equal for every column)
  unsigned int numRows = input_[columnNames_[0]].size();
  
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> result;
  result.resize(numRows, numColumns);
  
  // fill eigen matrix
  size_t i = 0;
  for(std::string key : columnNames_){
    for(size_t j = 0; j < input_[key].size(); ++j){
      result(j,i) = input_[key][j];
    }
    ++i;
  }
  
  return result;
}

#endif // __CSV_READER_H__
