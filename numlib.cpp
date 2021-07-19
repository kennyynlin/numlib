#include "numlib.h"
#include <vector>
#include <cmath>
#include <cfloat>
#include <omp.h>

dataset::dataset()
{
    xData = {0};
    yData = {0};
}

dataset::dataset(std::vector<double> xData, std::vector<double> yData)
{
    if (isEqual())
    {
        this->xData = xData;
        this->yData = yData;
    }
}
dataset::dataset(dataset &d)
{
    xData = d.xData;
    yData = d.yData;
}

void dataset::resetXYData(std::vector<double> xData, std::vector<double> yData)
{
    if (isEqual())
    {
        this->xData = xData;
        this->yData = yData;
    }
}
double dataset::leastSquareY(double x)
{
    double minSigma = DBL_MAX, tempSigma, y = 0;
    std::vector<double> tempCoef, bestCoef;

    for(size_t i = 1 ; i < xData.size() ; i++)
    {
        tempCoef = getPolyCoef((int)i);
        tempSigma = calSigma(tempCoef);
        if(tempSigma < minSigma)
        {
            minSigma = tempSigma;
            bestCoef = tempCoef;
        }
    }

    for(size_t i = 0 ; i < bestCoef.size() ; i++)
    {
        y = y + bestCoef[i]* pow(x, i);
    }

    return y;
}
std::vector<double> dataset::leastSquareCoef()
{
    double minSigma = DBL_MAX, tempSigma;
    std::vector<double> tempCoef, bestCoef;

    for(size_t i = 1 ; i < xData.size() ; i++)
    {
        tempCoef = getPolyCoef((int)i);
        tempSigma = calSigma(tempCoef);
        if(tempSigma < minSigma)
        {
            minSigma = tempSigma;
            bestCoef = tempCoef;
        }
    }
    return bestCoef;
}
double dataset::gregoryNewton(double x)
{
    double h = xData[1] - xData[0];
    double s = (x - xData[0])/h;
    double p;
    std::vector<std::vector<double>> diffTable;

    diffTable = calDiffTable();
    p = calPloy(diffTable, s);

    return p;
}
double dataset::calSumA(int power)
{
    double sum = 0;
    for(size_t i = 0 ; i < xData.size() ; i++)
    {
        sum = sum + pow(xData[i], power);
    }
    return sum;
}
double dataset::calSumB(int power)
{
    double sum = 0;
    for(int i = 0 ; i < xData.size() ; i++)
    {
        sum = sum + yData[i]*pow(xData[i], power);
    }
    return sum;
}
double dataset::calSigma(std::vector<double> coef)
{
    double sigma = 0, sum;
    for(int i = 0 ; i < xData.size()  ; i++)
    {
        sum = 0;
        for(int j = 0 ; j < coef.size() ; j++)
        {
            sum = sum + coef[j]*pow(xData[i], j);
        }

        sigma = sigma + (yData[i] - sum)*(yData[i] - sum);
    }
    sigma = sqrt(sigma/((int)xData.size() - ((int)coef.size() - 1)));
    return sigma;
}
std::vector<double> dataset::getPolyCoef(int m)
{
    std::vector<double> coef;
    std::vector<std::vector<double>> augMatrix(m+1, std::vector<double> (m+2, 0) );

    for(int i = m, ii = 0 ; i >= 0 && ii < m + 1; i--, ii++)
    {
        for(int j = m, jj = 0 ; j >= 0 && jj < m + 1 ; j--, jj++)
        {
            augMatrix[ii][jj] = calSumA(i + j);
        }
        augMatrix[ii][m+1] = calSumB(i);
    }

    //partial pivoting
    for (int i = 0; i < m + 1 ; i++)
    {
        for (int k = i + 1 ; k < m + 1; k++)
        {
            if (abs(augMatrix[i][i]) < augMatrix[k][i])
            {
                for (int j = 0; j <= m + 1 ; j++)
                {
                    double temp = augMatrix[i][j];
                    augMatrix[i][j] = augMatrix[k][j];
                    augMatrix[k][j] = temp;
                }
            }
        }
    }

    double d;
    for (int k = 0 ; k < m + 1 ; k++)
    {
        for (int i = 0; i < m + 1; i++)
        {
            if (k!=i)
            {
                d = augMatrix[i][k] / augMatrix[k][k];
                for(int j = 0; j < m + 2; j++)
                {
                    augMatrix[i][j] = augMatrix[i][j] - d*augMatrix[k][j];
                }
            }
        }
    }

    for (int i = m ; i > 0 ; i--)
    {
        coef.push_back(augMatrix[i][m+1] / augMatrix[i][i]);
    }

    return coef;
}
double dataset::comb(double s, int k)
{
    double denominator = 1, numerator = 1;
    if (k == 0)
        return 1;
    else
    {
        for (int i = 0; i < k; i++)
        {
            numerator = numerator*(s - i);
            denominator = denominator*(k - i);
        }

        return numerator / denominator;
    }
}
double dataset::calPloy(std::vector<std::vector<double>> diffTable, double s)
{
    double p = 0;

    for (size_t i = 0; i < xData.size(); i++)
    {
        p = p + comb(s, (int)i) * diffTable[0][i];
    }

    return p;
}
std::vector<std::vector<double>> dataset::calDiffTable()
{
    size_t n = xData.size();
    std::vector<std::vector<double>> diffTable(n, std::vector<double>(n, 0));

    for (size_t i = 0; i < n ; i++)
    {
        diffTable[i][0] = yData[i];
    }

    for(size_t i = 1 ; i < n ; i++)
    {
        #pragma omp parallel for num_threads(2) default(shared)
        for(size_t j = 0 ; j< n - i ; j++)
        {
            diffTable[j][i] = diffTable[j+1][i-1] - diffTable[j][i-1];
        }
    }

    return diffTable;
}