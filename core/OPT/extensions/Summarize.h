#ifndef __SUMMARIZE_H__
#define __SUMMARIZE_H__

#include <chrono>

class Summarize {

private:
  // used for timing
  std::chrono::steady_clock::time_point begin{};
  std::chrono::steady_clock::time_point end{};

  // prints an SVector (a.k.a. eigen vector) by row instead of by column
  std::string rowPrint(const std::vector<double>& input){
    std::string output = "[" + std::to_string(input[0]);

    for(size_t j = 1; j < input.size(); ++j){
      output += ", " + std::to_string(input[j]);
    }

    output += "]";
    return output;
  }
  
public:
  Summarize() = default;
  
  template <typename Optimizer, typename Objective>
  bool endOptimization(Optimizer& opt, Objective& obj){
    // stop timer
    end = std::chrono::steady_clock::now();
    // compute elapsed time
    auto time = std::chrono::duration_cast<std::chrono::microseconds>(end - begin);
    
    // transform SVector in standard vector for easy of printing
    const std::vector<double> v(&opt.getXnew()[0], &opt.getXnew()[0] + opt.getXnew().cols() - 1);
    
    std::cout << "Summary of optimization routine" << std::endl;
    std::cout << "   number of iterations:     " << opt.getNumIteration() << std::endl;
    std::cout << "   optimium point found:     " << rowPrint(v) << std::endl;
    std::cout << "   objective optimum value:  " << obj(opt.getXnew()) << std::endl;
    std::cout << "   l2 error:                 " << opt.getError() << std::endl;
    std::cout << "   time (milliseconds):      " << opt.getTime().count()/1000 << std::endl;
    
    return false;
  }
};

#endif // __SUMMARIZE_H__
