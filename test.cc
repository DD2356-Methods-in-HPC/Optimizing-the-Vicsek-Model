#include <iostream>
#include <assert.h>

#include "simulation.h"


void test_add() {
    std::cout << "Testing add: 1+2" << std::endl;
    int res = add(1, 2);
    assert(res == 3);
}


int main() {

    std::cout << "Running tests..." << std::endl;
    test_add();
    std::cout << "All tests passed!" << std::endl;

    
  
}