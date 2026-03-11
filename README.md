# Computer-assisted-proof-for-Stokes-constants-in-generalized-standard-maps
Complementary code for Chapter 5 titled "Stokes constants: a computer assisted proof" of my doctoral thesis with title "Splitting of separatrices in generalized standard maps"

This project uses the CAPD::DynSys library.

Official repository: https://github.com/CAPDgroup/CAPD

CAPD main webpage: http://capd.ii.uj.edu.pl

# User's guide
Clone the respository:
```
git clone https://github.com/DidacGilRams/Computer-assisted-proof-for-Stokes-constants-in-generalized-standard-maps
```
Enter the repository, clone, compile and install CAPD:
```
cd Computer-assisted-proof-for-Stokes-constants-in-generalized-standard-maps
cd CODE
make setup
```
Run the programs:
```
./bin/PI1_N4
```
Compile the programs after modification:
```
make
```
# Content of the programs
```Stokes()``` executed in ```main()``` produces the respective computer assisted proofs.
```Stokes()``` contains the parameters to play with described below:
- PI1_N4.cpp: deals with the polynomial case with one monomial using an initial aproximation with 4 terms. Contains the parameters ```dk0,fk0,P,L,rho``` and ```gamma```. By default, ```dk0=2,fk0=1,P=1000,L=1000,rho=8.215``` and ```gamma=0.949327```.
- PI1_N5.cpp: deals with the polynomial case with one monomial using an initial aproximation with 5 terms. Contains the parameters ```dk0,,fk0,P,L,rho``` and ```gamma```. By default, ```dk0=2,fk0=1,P=10000,L=500,rho=8.171``` and ```gamma=0.949327```.
- PI2_3_N7.cpp: deals with the polynomial case with two monomials of degrees 2 and 3 using an initial aproximation with 7 terms. Contains the parameters ```a_2,fk0,P,L,rho``` and ```gamma```. By default, ```a_2=6.02,fk0=1,P=10000,L=200,rho=9.155``` and ```gamma=0.949327```.
- TI1_N4.cpp: deals with the trigonometric case with one monomial using an initial aproximation with 4 terms. Contains the parameters ```dk0,P,L,rho``` and ```gamma```. By default, ```dk0=1,P=20000,L=1000,rho=5.475``` and ```gamma=0.949327```.
