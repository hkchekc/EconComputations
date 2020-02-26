#include <vector>
#include "authority.hpp"

using namespace std;

class Economy{
public:
    vector<double> ratioanlity; // size 4, representing the four quantites
    Authority authority;
    int no_firm;
    Economy(int nk, double k1, double k2);
    double consumerVal, socialVal;
    double bellman();
}