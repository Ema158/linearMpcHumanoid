#include <iostream>
#include "controller/controller.hpp"
#include "controller/robotParameters.hpp"
#include "controller/robotInfo.hpp"
#include <iostream>

int main() {
    robotInfo nao;
    std::vector<linkInertia> T;// = nao.initialConfiguration();
    T = nao.getLinks();
    for (int i=0;i<nao.getNumBodies();i++){
            std::cout << T[i].com << std::endl;
    }
    
    return 0;
}