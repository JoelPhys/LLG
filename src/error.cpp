#include "../inc/error.h"
#include <iostream>
#include <cstdlib>
#include <string>
#include <sstream>

namespace error {

    void Info(std::string f,int l)
    {
        std::cerr << "\n***Error in File: " << f << ", at line: " << l <<std::flush;
    }
    void Message(std::string em)
    {
        std::cerr << ". " << em << "***\n" << std::endl;
        exit(0);
    }
}