# numlib
This is a C++ numerical method library contains Least Square Polynomial Approximation, Gregory-Newton Interpolation, and Numerical Integration.

## Usage
```C++
#include "numlib.h"
```
### Least Square Polynomial Approximation
```C++
std::vector<double> xdata{1, 2, 3, 4, 5};
std::vector<double> ydata{1, 8, 9, 64, 125};
dataset ds(xdata, ydata);

# returns polynomial coefficient
std::vector<double> coef;
coef = ds.leastSquareCoef();
```
coef should be {0, 0, 0, 1} which represents 0*x^0+0*x^1+0*x^2+1*x^3.
```C++
# returns y value
double yValue;
yValue = ds.leastSquareY(6);
```
yValue should be 6^3 = 216.
### Gregory-Newton Interpolation
Please refer to the explanation of [Gregory-Newton Forward Difference](https://nptel.ac.in/content/storage2/courses/122104019/numerical-analysis/Rathish-kumar/rathish-oct31/fratnode8.html).
```C++
#returns interpolation value
yValue = ds.gregoryNewton(6);
```
yValue should also be 216.

## Contributing
Pull requests are welcome! Feel free to contact me via email.
