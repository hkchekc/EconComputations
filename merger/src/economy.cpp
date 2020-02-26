#include <vector>
#include "authority.hpp"
#include "economy.hpp"

using namespace std;

 Economy::Economy(int no_firm, double k1, double k2){
     
}

double bellman(){
    int rowLoc;
    double consumerVal, socialVal, producerVal;

    consumerVal = priceNew - priceOld;
    producerVal = valNew - valOld;
    socialVal = consumerVal + ;

    return consumerVal, socialVal;
}