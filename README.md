# SCATCUB_SRBF

The main function of this repository is SCATCUB_SRBF.m which consists in the algorithm described in the paper

_Roberto Cavoretto, Alessandra De Rossi, Giacomo Elefante and Alvise Sommariva - "Adaptive RBF cubature by scattered data on spherical polygons"_

Moreover there are some demos to reproduce the various numerical experiments included in the paper.

**It is necessary the use of the package CHEBFUN to run the code**

In particular:
 ```
1. demo_adaptive  - demo for the numerical experiments on the adaptive algorithm for integration. It might include two input parameters f_type, domain_type
                    domain_type : to choose the spherical polygon and it takes values 1,2,3 - 1: Africa region; 2: Australian region; 3: America region
                    f_type      : to choose the integrand function it goes from 1 to 33
                                  the test functions of the paper are:     1. Africa region    : f1: 28, f2: 27, f3: 29
                                                                           2. Australia region : f1: 30, f2: 27, f3: 31
                                                                           3. America regions  : f1: 32, f2: 27, f3: 33
```
```
  2. demo_cond   - demo for the numerical experiments on the conditioning of the Vandermonde-like matrix. It might include 4 input parameters rbf_type, domain, do_plot, isdeg
                   rbf_type    : to choose the rbf basis and it takes value 1 or 2 - 1:  RP; 2: TPS 
                   domain      : to choose the domain of the experiments, it takes in input a polyshape polygons
                   do_plot     : to reproduce a plot of the experiments - 1: yes; 0: no
                   is_deg      : to express if the coordinates of the polyshape polygon are expressed in degree - 1: degree; 0: radiant
```
```
3. doPlotErrSphere - demo to reproduce the plot of interpolation error between the function and the rbf interpolant with 1600 centers. It might include 3 input parameters domain_type, f_type, rbf_type
                     domain_type : to choose the spherical polygon and it takes values 1,2,3 - 1: Africa region; 2: Australian region; 3: America region
                     f_type      : to choose the integrand function it goes from 1 to 33
                                   the test functions of the paper are:     1. Africa region    : f1: 28, f2: 27, f3: 29
                                                                            2. Australia region : f1: 30, f2: 27, f3: 31
                                                                            3. America regions  : f1: 32, f2: 27, f3: 33
                     rbf_type    : to choose the rbf basis and it takes value 1 or 2 - 1:  RP; 2: TPS
```
```
4. doPlotErrVSLOOCV  - demo to reproduce the integration error of the rbf varying the exponents of the rbf. It might include 3 input paramenters domain_type,rbf_type, f_type
                       domain_type : to choose the spherical polygon and it takes values 1,2,3 - 1: Africa region; 2: Australian region; 3: America region
                       f_type      : to choose the integrand function it goes from 1 to 33
                                     the test functions of the paper are:     1. Africa region    : f1: 28, f2: 27, f3: 29
                                                                              2. Australia region : f1: 30, f2: 27, f3: 31
                                                                              3. America regions  : f1: 32, f2: 27, f3: 33
                       rbf_type    : to choose the rbf basis and it takes value 1 or 2 - 1:  RP; 2: TPS
```
```
5. doPlotSphere  - demo to reproduce the graphic of the function evaluating it in the spherical polygon. It might include 2 input parameters domain_type,f_type
                   domain_type : to choose the spherical polygon and it takes values 1,2,3 - 1: Africa region; 2: Australian region; 3: America region
                   f_type      : to choose the integrand function it goes from 1 to 33
                                 the test functions of the paper are:     1. Africa region    : f1: 28, f2: 27, f3: 29
                                                                          2. Australia region : f1: 30, f2: 27, f3: 31
                                                                          3. America regions  : f1: 32, f2: 27, f3: 33
   ```
