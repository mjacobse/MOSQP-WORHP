#include "algorithm/MOSQP.hpp"
#include "test_problems/TestProblems.hpp"
#include <iostream>


int main()
{    
    mosqp::MONLP &test = test_problems::GE3();
    mosqp::MOSQP mosqp(test);
    mosqp::ParetoFront front = mosqp.Solve();
    std::cout << "===========================================================" << std::endl;
    std::cout << "Objective evaluations: " << static_cast<double>(test.GetNumEvalF()) / test.GetNumObjectives() << std::endl;
    std::cout << "Gradient evaluations: " << static_cast<double>(test.GetNumEvalDF()) / test.GetNumObjectives() << std::endl;
    std::cout << "D2F evaluations: " << static_cast<double>(test.GetNumEvalD2F()) / test.GetNumObjectives() << std::endl;
    std::cout << "D2G evaluations: " << test.GetNumEvalD2G() << std::endl;
    std::cout << "Constraint evaluations: " << test.GetNumEvalG() << std::endl;
    std::cout << "Jacobian evaluations: " << test.GetNumEvalDG() << std::endl << std::endl;
    
    system("pause");
    return 0;
}
