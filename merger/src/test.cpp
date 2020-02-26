# include <iostream>
#include <stdio.h>
#include <cmath>
#include <math.h>
#include <vector>
#include <algorithm>
# include <array>
using namespace std;

const int len = 50;
double a[len] = {1010.0, 5., 200.};
string testStr = "the skju";
array<double, len> anotherArr = {1,2,3,4,5};
string test;

double loop(int iter, double (*func)(double)){
    double test_double = 1.5;
    double second_double = 0.0;
    for (;iter>0;iter--){
        second_double =  (*func)(second_double);
        cout<< to_string(second_double) << "\n";
    }
    return second_double;
}

double add_up(double dou){
    dou += 0.2345;
    return dou;
}

vector<double>* add_dim(vector<double>* arr){
    vector<double>* vec;
    vec = new vector<double>;
    vec = 0;
    return vec;
}


int main(){
    struct Point {
    int x;
    int y;
    double third_double;
    };

    Point* p;      // declare pointer to a Point struct

    p = new Point; // dynamically allocate a Point
    p->x = 12;  // set the field values.
    p->y = 34;
    p->third_double = loop(20, add_up) ;
    cout<< to_string((*p).third_double) << "\n";
    return 0;
}


