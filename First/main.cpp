#include <cmath>
#include <iostream>
#include <string>

template <class T>
void findULP(const std::string& precision)
{
    T Epsilon = 1.;
    while (1 + Epsilon / 2 != 1 ) {
        Epsilon /= 2;
    }
    std::cout << precision <<" precision ULP: " << Epsilon << std::endl;

    //std::cout << 1 + Epsilon << std::endl << 1 + Epsilon + Epsilon/2 << std::endl << std::endl;
}

template <class T>
void findMaxNum(const std::string& precision)
{
    T MaxNum = 1.;
    unsigned int count = 0;
    while (!std::isinf(MaxNum*2))
    {
        MaxNum*=2;
        count++;
    }
    std::cout << precision <<" precision MaxNum: " << std::scientific << MaxNum << std::endl;
    std::cout << precision <<" precision std MaxNum: " << std::scientific << std::numeric_limits<T>::max() << std::endl;
    std::cout << precision <<" precision Count: "<< count << std::endl;
}

template <class T>
void findMinNum(const std::string& precision)
{
    T MinNum = 1.;
    unsigned int count = 0;
    while (MinNum/2 != 0)
    {
        MinNum/=2;
        count++;
    }
    std::cout << precision <<" precision MinNum: " << std::scientific << MinNum << std::endl;
    std::cout << precision <<" precision std MinNum: " << std::scientific << std::numeric_limits<T>::min() << std::endl;
    std::cout << precision <<" precision Count: "<< count << std::endl;
}

template <class T>
void findSum(const std::string& precision)
{
    std::cout << precision << " precision Sums: " << std::endl;
    T fromGreaterToLesserSum = 0.;
    for (size_t n = 10000; n >= 1; --n) {
        fromGreaterToLesserSum += pow(-1, n) / n;
    }
    std::cout << fromGreaterToLesserSum << std::endl;

    T fromLesserToGreaterSum = 0.;
    for (size_t n = 1; n <= 10000; ++n) {
        fromLesserToGreaterSum += pow(-1, n) / n;
    }
    std::cout << fromLesserToGreaterSum << std::endl;

    T fromGreaterToLesserPositiveSum = 0.;
    T fromGreaterToLesserNegativeSum = 0.;
    for (size_t n = 10000; n >= 1;)
    {
        fromGreaterToLesserNegativeSum += 1. / n;
        --n;
        fromGreaterToLesserPositiveSum += -1. / n;
        --n;
    }
    std::cout << fromGreaterToLesserNegativeSum + fromGreaterToLesserPositiveSum << std::endl;

    T fromLesserToGreaterPositiveSum = 0.;
    T fromLesserToGreaterNegativeSum = 0.;
    for (size_t n = 1; n <= 10000;)
    {
        fromLesserToGreaterNegativeSum += -1. / n;
        ++n;
        fromLesserToGreaterPositiveSum += 1. / n;
        ++n;
    }
    std::cout << fromLesserToGreaterNegativeSum + fromLesserToGreaterPositiveSum << std::endl;
}
template <class T>
void findAll(const std::string& precision)
{
    findULP<T>(precision);
    findMinNum<T>(precision);
    findMaxNum<T>(precision);
    std::cout << std::endl;
}

int main()
{
    std::cout.precision(64);

    std::cout << "a)" << std::endl;

    findAll<float>("Single");
    findAll<double>("Double");

    std::cout.precision(16);

    std::cout << "b)" << std::endl;
    findSum<float>("Single");
    findSum<double>("Double");

    return 0;
}
