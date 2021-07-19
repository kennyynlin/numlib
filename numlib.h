#ifndef NUMLIB_NUMLIB_H
#define NUMLIB_NUMLIB_H
#include <vector>

class dataset
{
    private:
    std::vector<double> xData;
    std::vector<double> yData;
    bool isEqual() {return xData.size() == yData.size();}
    //Least Square
    double calSumA(int power);
    double calSumB(int power);
    double calSigma(std::vector<double> coef);
    std::vector<double> getPolyCoef(int m);
    //Gregory Newton
    double comb(double s, int k);
    double calPloy(std::vector<std::vector<double>> diffTable, double s);
    std::vector<std::vector<double>> calDiffTable();

    public:
    dataset();
    dataset(std::vector<double> x, std::vector<double> y);
    dataset(dataset &d);
    std::vector<double> getXData(){return xData;}
    std::vector<double> getYData(){return yData;}
    void resetXYData(std::vector<double> xData, std::vector<double> yData);
    //Least Square
    double leastSquareY(double x);
    std::vector<double> leastSquareCoef();
    //Gregory Newton
    double gregoryNewton(double x);
};

#endif //NUMLIB_NUMLIB_H
