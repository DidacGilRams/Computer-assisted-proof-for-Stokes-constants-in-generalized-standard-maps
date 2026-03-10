#include "capd/capdlib.h"
#include <iostream>
#include <chrono>
#include <ctime>
using namespace std;
using namespace capd;
#include "capd/matrixAlgorithms/lib.h"
#include "capd/matrixAlgorithms/floatMatrixAlgorithms.hpp"
using namespace capd::alglib;
using namespace capd::matrixAlgorithms;

//----------------------------------- Trigonometric case with one monomial (Section 5.5.1) ----------------------------------------------------------

//This code gives an interval for the Stokes constant p_2^{-1}(d_{k_0}) for every degree d_{k_0} (has to be introduced in the Stokes() function).
//This version uses a first approximation for phi_0 with 4 terms, that is, psi_0=psi_0^0+psi_0^1+psi_0^2+psi_0^3+psi_0^4

//<<<<<< Implementation of some functions used in the program >>>>>>

interval angle(interval x,interval y)// Function that computes the angle on the plane
{
    if((y>=0) && (2*x>=y))           return  atan(y/x);
	if((y>=0) && (abs(x)<=2*abs(y))) return -atan(x/y)+interval(interval(0.5))*interval::pi();
	if((x<=0) && (2*abs(x)>=abs(y))) return  atan(y/x)-interval::pi();
	if((y<=0) && (abs(x)<=2*abs(y))) return -atan(x/y)-interval(interval(0.5))*interval::pi();
	if((y<0) && (2*x>=-y))           return  atan(y/x);
	if((x>=0) && (2*abs(x)>=abs(y))) return  atan(y/x);
	cout << "x = " << x << endl;
	cout << "y = " << y << endl;
	cout << "failed to compute the angle." << endl;
	abort();
	return interval(0);
}

interval mod2Pi(interval a)// Function that gives the 2pi modulus
{
	if(a>interval::pi()) return mod2Pi(a-interval(2)*interval::pi());
	if(a<-interval::pi()) return mod2Pi(a+interval(2)*interval::pi());
	return a;
}

interval angle(IVector z)// Function that gives the angle of a complex number
{
        return angle(z[0],z[1]);
}

interval mod(IVector z)// Modulus of a complex number
{
	return sqrt(power(z[0],2)+power(z[1],2));
}

IVector cprod(IVector z1,IVector z2)// Complex product
{
	IVector y(2);
	y[0]=z1[0]*z2[0]-z1[1]*z2[1];
	y[1]=z1[0]*z2[1]+z1[1]*z2[0];
	return y;
}

IVector log(IVector z)// Complex logarithm
{
        interval r=mod(z);
		IVector w(2);
        w[0] = log(r);
        w[1] = angle(z);
        return w;
}

IVector exp(IVector z)// Complex exponential
{
		interval r=exp(z[0]);
        IVector w(2);
        w[0]=r*cos(z[1]);
        w[1]=r*sin(z[1]);
        return w;
}

IVector power(IVector z,interval n)// Complex power interval n of z
{
	interval r,theta;
	r=mod(z);
	r=power(r,n);
	theta=angle(z);
	z[0]=r*cos(n*theta);
	z[1]=r*sin(n*theta);
	return z;
}

int floor(interval n)//Floor of the left bound of an interval
{
	int p=0;
	for (int i=1;i<leftBound(n)+1;i++)
	{
		p=i;
	}
	return p;
}

interval fac(int n)// Factorial
{
	interval j=interval(1);
	if (n<=0)
	{
		return 1;
	}
	else
	{
		for(int i=1;i<=n;i++)
		{
			j=j*interval(i);
		}
		return j;
	}
}

interval dfac(int n)// Double factorial
{
	interval j=interval(n);
	if (n<=0)
	{
		return 1;
	}
	else
	{
		for(int i=1;2*i<n;i++)
		{
			j=j*interval(n-2*i);
		}
		return j;
	}
}

interval binom(interval a,int k)//Generalized binomial
{
	interval j=a,x=a/interval(k);
	for(int i=1;i<k;i++)
	{
		x=x*((a-i)/(k-i));
	}
	return x;
}

void split(IVector x,IVector &d,IVector &psi)// Function that split vector x (4 dimensional) in vectors d and psi (two dimensional)
{
	d[0]=x[0];
	d[1]=x[1];
	psi[0]=x[2];
	psi[1]=x[3];
}

void dsplit(IVector dx,IVector &d,IVector &psi,IVector &dd,IVector &dpsi)// Function that split vector dx (8 dimensional) in vectors d, psi, dd, dpsi (two dimensional)
{
	d[0]=dx[0];
	d[1]=dx[1];
	psi[0]=dx[2];
	psi[1]=dx[3];
	dd[0]=dx[4];
	dd[1]=dx[5];
	dpsi[0]=dx[6];
	dpsi[1]=dx[7];
}

IVector combine(IVector d,IVector psi)// Inverse of split
{
	IVector x(4);
	x[0]=d[0];
	x[1]=d[1];
	x[2]=psi[0];
	x[3]=psi[1];
	return x;
}

IVector dcombine(IVector d,IVector psi,IVector dd,IVector dpsi)// Inverse of dsplit
{
	IVector x(8);
	x[0]=d[0];
	x[1]=d[1];
	x[2]=psi[0];
	x[3]=psi[1];
	x[4]=dd[0];
	x[5]=dd[1];
	x[6]=dpsi[0];
	x[7]=dpsi[1];
	return x;
}

//<<<<<< Obtaining rho_0 for the existence of solutions of the inner equation (Section 5.2) >>>>>>

//Norms of the inverse operators

interval nu0(int dk0)//Definition of nu_0 in the trigonometric case
{
	return interval(2.)/interval(dk0); 
}

int pNn(int dk0,int N,int n)//Definition of p_N^n
{
	interval x;
	x=interval(N-n)/interval(dk0);
	return floor(x)+1;
}

interval B(interval nu)//Defines the constants B(nu) used to bound the inverse S (Proposition 5.2.8)
{
	int n=floor(nu);
	interval x,Pi=interval::pi();
	x=interval(dfac(n-2)/dfac(n-1));
	if (power(-1,n)==1){
		x=x;
	}
	else{
		x=x*Pi/interval(2);
	}
	return x;
}

interval S(interval rho,interval gamma,interval nu)// Bound of the first order inverse operator S(Proposition 5.2.8)
{
	interval s=interval(0.),k;
	for(int i=1;i<=4;i++)
	{
		k=(interval(4.)+power(2*gamma+rho,2))/(power(interval(2-i),2)+power(2*gamma+rho,2));
		k=power(k,nu/interval(2.));
		s=s+k;
	}
	return (s/rho)+B(nu);
}

interval G(interval rho,interval gamma)// Bound of the second order inverse operator G (Proposition 5.2.8)
//Here for the trigonometric case kappa=0 and since we took psi00+psi01+psi02+psi03+psi04 as approximation nu_N=10
{
	interval x;
	x=S(rho,gamma,interval(3));
	x=x*(S(rho,gamma,interval(9))+S(rho,gamma,interval(12)));
	return x;
}

//The following fuctions are used to compute the norm of the first iteration.

interval E0(int dk0,IVector b,interval rho)//The sum of E_N^n for all 0<=n<=N (Lemma 5.2.2)
{ 
	interval psi00,psi01,psi02,psi03,psi04;
	psi00=1-power(rho,-2);
	psi00=log(psi00)*power(rho,12);
	for(int i=1;i<6;i++)
	{
		psi00=psi00+power(rho,12-(2*i))/i;
	}
	psi00=nu0(dk0)*abs(psi00);
	psi01=power(rho-1,-2)+power(rho+1,-2)-2/power(rho,2);
	psi01=psi01*power(rho,12);
	for(int i=1;i<5;i++)
	{
		psi01=psi01-2*power(rho,12-2*(i+1))*(2*(i+1)-1);
	}
	psi01=abs(b[0]*psi01);
	psi02=power(rho-1,-4)+power(rho+1,-4)-2/power(rho,4);
	psi02=psi02*power(rho,12);
	for(int i=1;i<4;i++)
	{
		psi02=psi02-2*power(rho,12-2*(i+2))*(fac(2*(i+2)-1)/(fac(3)*fac(2*i)));
	}
	psi02=abs(b[1]*psi02);
	psi03=power(rho-1,-6)+power(rho+1,-6)-2/power(rho,6);
	psi03=psi03*power(rho,12);
	for(int i=1;i<3;i++)
	{
		psi03=psi03-2*power(rho,12-2*(i+3))*(fac(2*(i+3)-1)/(fac(5)*fac(2*i)));
	}
	psi03=abs(b[2]*psi03);
	psi04=power(rho-1,-8)+power(rho+1,-8)-2/power(rho,8);
	psi04=psi04*power(rho,12);
	for(int i=1;i<2;i++)
	{
		psi04=psi04-2*power(rho,12-2*(i+4))*(fac(2*(i+4)-1)/(fac(7)*fac(2*i)));
	}
	psi04=abs(b[3]*psi04);
	return psi00+psi01+psi02+psi03+psi04;
}

interval CN(int dk0,IVector b,int a,interval rho)//Definition of C_N(a,rho) (see Lemma 5.2.4)
{
	interval c=interval(0.);
	for(int n=1;n<=4;n++)
	{
		c=c+abs(b[n-1])*power(rho,a*nu0(dk0)-2*n);//Absolute values of psi_0^n evaluated at rho
	}
	return c;
}

interval f0(int dk0,IVector b,interval rho)//The sum of F_N^k (Lemma 5.2.5)
{
	interval c,psi01,psi02,psi03,psi04;
	psi01=abs(b[0])*power(rho,-2);//Absolute values of psi0i evaluated at rho
	psi02=abs(b[1])*power(rho,-4);
	psi03=abs(b[2])*power(rho,-6);
	psi04=abs(b[3])*power(rho,-8);
	c=CN(dk0,b,0,rho);
	c=exp(dk0*c)*power(dk0*c,5)/fac(5);//Part of the remainder
	c=c+pow(dk0,2)*(psi01*psi04+psi02*(psi03+psi04)+power(psi03+psi04,2)/2);//Quadratic part
	c=c+pow(dk0,3)*(power(psi01,2)*(psi03+psi04)/2+psi01*power(psi02+psi03+psi04,2)/2+power(psi02+psi03+psi04,3)/fac(3));//Cubic part
	c=c+pow(dk0,4)*(power(psi01,3)*(psi02+psi03+psi04)/fac(3)+power(psi01,2)*power(psi02+psi03+psi04,2)/interval(4)+psi01*power(psi02+psi03+psi04,3)/fac(3)+power(psi02+psi03+psi04,4)/fac(4));//4th power part
	c=nu0(dk0)*c*power(rho,10);//(1/rho^2)*c multiplied by rho^10 coming from the norm 
	return c;
}

interval F0(int dk0,IVector b,interval rho,interval gamma)//Computation of G*E_0^b (Corollaries 5.2.7 and 5.2.9)
{
	return G(rho,gamma)*(E0(dk0,b,rho)+f0(dk0,b,rho));
}

//The following are used to compute the Lipschidz constant (Lemma 5.2.11)

interval A(int dk0,IVector b,interval rho)//Boundig the norm of the linear term A (Lemma 5.2.10)
{
	return nu0(dk0)*pow(dk0,2)*exp(dk0*CN(dk0,b,0,rho))*CN(dk0,b,0,rho);
}

interval R(int dk0,IVector b,interval rho,interval M0)//Bounding higher order terms R^b (Lemma 5.2.11)
{
	return nu0(dk0)*dk0*exp(dk0*CN(dk0,b,0,rho))*(exp((dk0*M0)/power(rho,10))-interval(1.));//multiplied by rho^2 from the inverse G
}

interval M0(int dk0,IVector b,interval rho,interval gamma)//Function M_0(rho,gamma)
{
	return 2*F0(dk0,b,rho,gamma);
}

interval Lip(int dk0,IVector b,interval rho,interval gamma)//The Lipschidz constant for the existence (Lemma 5.2.11)
{
	return G(rho,gamma)*(A(dk0,b,rho)+R(dk0,b,rho,M0(dk0,b,rho,gamma)));
}


//<<<<<< Obtaining rho_1 for the existence of zeta_1 and explicit formula for psi^u-psi^s (Section 5.3) >>>>>>

//Norms of the inverse operators

interval c(interval gamma)//Definition of c(gamma) involved in the domains (see definition of the domain  D_{rho,gamma}^delta)
{
	interval x;
	x=power(gamma,2)+1;
	x=power(x,interval(0.5));
	x=gamma/x;
	return x;
}

interval Sd(interval rho,interval gamma,interval nu)// Bound of the first order inverse operator S^d  (Proposition 5.3.2)
{
	interval r,x1,x2,Pi=interval::pi();
	r=rho+(3*gamma)/2;
	if (leftBound(nu)>0){
		x1=power(1-c(gamma)/2,2)/power(rho+gamma*(5-c(gamma))/2,2);
		x1=power(1+x1,nu/2);
		x1=x1/(2*Pi*power(c(gamma),2)*r);
		x2=(3-2*c(gamma))/(4*(power(1-c(gamma)/2,2)+power(rho+gamma*(5-c(gamma))/2,2)));
		x2=power(1-x2,-nu/2);
		x2=(1+power(c(gamma),-1))*B(nu)*x2;
		x1=x1+x2;
	}
	else{
		if (leftBound(nu)<=-1){
			x1=1-(abs(nu)-1)/(2*Pi*r);
			x1=power(c(gamma),-1)+power(x1,-1);
			x1=1+x1/(2*Pi*r);
			x2=power(1-c(gamma)/2,2)+power(rho+gamma*(5-c(gamma))/2,2);
			x2=1+(5-2*c(gamma))/(4*x2);
			x2=power(x2,abs(nu)/2);
			x1=x1*x2;
		}
		else{
			x1=(1+power(c(gamma),2))/(2*Pi*power(c(gamma),2)*r);
			x1=x1+power(abs(nu)*c(gamma),-1);
		}
	}
	x1=power(1+power(r,-1),abs(nu))*x1+power(r,-1);
	return x1;
}

//In order to define the bound for the last inverse, we need to define the the following two bounds related to eta_1^d and eta_2^d

interval eta1(int dk0,IVector b,interval rho,interval gamma)//Computation of eta_1^b (Lemma 5.3.3)
{
	interval r,x;
	r=(rho+3*gamma/2);
	x=nu0(dk0);
	for(int n=1;n<=4;n++)
	{
		x=x+2*n*abs(b[n-1])*power(rho,-2*n);//Absolute values of dpsi_0^n evaluated at rho
	}
	x=x+(2*M0(dk0,b,rho,gamma))/(c(gamma)*power(rho,9));
	return x+(M0(dk0,b,rho+gamma/2,gamma))/(power(r,10));
}

interval eta2(int dk0,IVector b,interval rho,interval gamma)//Computation of eta_2^b (Corollary 5.3.6)
{
	interval r,x,eta;
	r=(rho+3*gamma/2);
	eta=eta1(dk0,b,rho,gamma);
	x=1+power(r,-1);
	x=Sd(rho,gamma,-interval(3))*eta*power((2*nu0(dk0))-eta,-2)*x;
	return x;
}

interval Gd(int dk0,IVector b,interval rho,interval gamma)// Bound of the second order inverse operator G^d (Corollary 5.3.7)
//Here for the trigonometric case kappa=0 and since we took psi00+psi01+psi02+psi03+psi04 as approximation nu=10
{
	return eta1(dk0,b,rho,gamma)*eta2(dk0,b,rho,gamma)*(Sd(rho,gamma,interval(10))+Sd(rho,gamma,interval(13)));
}

interval Ad(int dk0,IVector b,interval rho, interval gamma)//Computation of A^{d,b} (Lemma 5.3.1)
{
	interval r,c;
	r=(rho+3*gamma/2);
	c=CN(dk0,b,0,r)+M0(dk0,b,rho+gamma/2,gamma)*power(r,-10);
	return nu0(dk0)*2*pow(dk0,2)*M0(dk0,b,rho+gamma/2,gamma)*exp(dk0*c);
}

interval Lipd(int dk0,IVector b,interval rho,interval gamma)//Bound for the operator G^d[A^d*] (Lemma 5.3.8)
{
	return Gd(dk0,b,rho,gamma)*Ad(dk0,b,rho,gamma)*power(rho+3*gamma/2,-10);
}

interval zeta11(int dk0,IVector b,interval rho,interval gamma)//Computation of zeta_1^{1,b} (Lemma 5.3.8)
{
	return (eta1(dk0,b,rho,gamma)*Gd(dk0,b,rho,gamma)*Ad(dk0,b,rho,gamma))/(1-Lipd(dk0,b,rho,gamma));
}

//<<<<<< Rigorous computations of psi^u-psi^s and dpsi^s (Section 5.4) >>>>>>

//Initial aproximation of the manifold and its derivative

IVector psi0(int dk0,IVector b,IVector z)//Computation of psi0=psi00+psi01+psi02+psi03+psi04
{
	IVector x(2);
	x=-nu0(dk0)*log(z);
	for(int n=1;n<=4;n++)
	{
		x=x+b[n-1]*power(z,-2*n);//psi_0^n evaluated at z
	}
	return x;
}

IVector psi0(int dk0,IVector b,interval re,interval im)//Computation of psi0(re,im)
{
	IVector z(2);
	z[0]=re;
	z[1]=im;
	return psi0(dk0,b,z);
}

IVector dpsi0(int dk0,IVector b,IVector z)// Derivative of psi0=psi00+psi01+psi02+psi03+psi04
{
	IVector x(2);
	x=-nu0(dk0)*power(z,-1);
	for(int n=1;n<=4;n++)
	{
		x=x+(-2*n)*b[n-1]*power(z,-(2*n+1));//dpsi_0^n evaluated at z
	}
	return x;
}

IVector dpsi0(int dk0,IVector b,interval re,interval im)// Derivative of psi_0(re,im)
{
	IVector z(2);
	z[0]=interval(re);
	z[1]=im;
	return dpsi0(dk0,b,z);
}

//<<<<<< Iterative methods for psi^u, psi^s (Section 5.5.1) and derivative psi^s (Section 5.5.2) >>>>>>

//The following functions are defined to obtain the generating functions F and dF for the iterative methods in Section 5

int omega(int dk0)//Computation of Omega(d_{k_0})
{
	return 20*dk0;
}

IVector T0(int dk0,IVector b,IVector z)//T_0^Omega(z) for psi00+psi01+psi02+psi03+psi04 (Lemma 5.4.4)
{
	int Omega=omega(dk0);
	IVector psi00(2),psi01(2),psi02(2),psi03(2),psi04(2);
	psi00[0]=interval(0.);
	psi00[1]=psi00[0];
	psi01=psi00;
	psi02=psi00;
	psi03=psi00;
	psi04=psi00;
	for(int i=6;i<=pNn(dk0,Omega,0);i++)
	{
		psi00=psi00+power(z,-2*i)/(2*i);
	}
	psi00=nu0(dk0)*psi00;
	for(int i=5;i<=pNn(dk0,Omega,1);i++)
	{
		psi01=psi01+power(z,-2*(i+1))*(2*(i+1)-1);
	}
	psi01=b[0]*psi01;
	for(int i=4;i<=pNn(dk0,Omega,2);i++)
	{
		psi02=psi02+power(z,-2*(i+2))*(fac(2*(i+2)-1)/(fac(3)*fac(2*i)));
	}
	psi02=b[1]*psi02;
	for(int i=3;i<=pNn(dk0,Omega,3);i++)
	{
		psi03=psi03+power(z,-2*(i+3))*(fac(2*(i+3)-1)/(fac(5)*fac(2*i)));
	}
	psi03=b[2]*psi03;
	for(int i=2;i<=pNn(dk0,Omega,4);i++)
	{
		psi04=psi04+power(z,-2*(i+4))*(fac(2*(i+4)-1)/(fac(7)*fac(2*i)));
	}
	psi04=b[3]*psi04;
	return -2*(psi00+psi01+psi02+psi03+psi04);
}

IVector T0b(int dk0,IVector b,IVector z)//Sum of T_0^{n,b}(z) (Lemma 5.4.4 and Remark 5.4.5)
{
	int Omega=omega(dk0);
	interval psi00,psi01,psi02,psi03,psi04,c;
	IVector x(2);
	c=mod(z);
	psi00=nu0(dk0)*((2*pNn(dk0,Omega,0)-1)/fac(2*(pNn(dk0,Omega,0)+1)))*(power(c-1,-(2*pNn(dk0,Omega,0)))+power(c+1,-(2*pNn(dk0,Omega,0)))-2*power(c,-(2*pNn(dk0,Omega,0))));
	psi01=abs(b[0])*(2*pNn(dk0,Omega,1)+1)*(power(c-1,-(2+2*pNn(dk0,Omega,1)))+power(c+1,-(2+2*pNn(dk0,Omega,1)))-2*power(c,-(2+2*pNn(dk0,Omega,1))));
	psi02=abs(b[1])*binom(2*2+2*pNn(dk0,Omega,2)-1,2*pNn(dk0,Omega,2))*(power(c-1,-(2*2+2*pNn(dk0,Omega,2)))+power(c+1,-(2*2+2*pNn(dk0,Omega,2)))-2*power(c,-(2*2+2*pNn(dk0,Omega,2))));
	psi03=abs(b[2])*binom(2*3+2*pNn(dk0,Omega,3)-1,2*pNn(dk0,Omega,3))*(power(c-1,-(2*3+2*pNn(dk0,Omega,3)))+power(c+1,-(2*3+2*pNn(dk0,Omega,3)))-2*power(c,-(2*3+2*pNn(dk0,Omega,3))));
	psi04=abs(b[3])*binom(2*4+2*pNn(dk0,Omega,4)-1,2*pNn(dk0,Omega,4))*(power(c-1,-(2*4+2*pNn(dk0,Omega,4)))+power(c+1,-(2*4+2*pNn(dk0,Omega,2)))-2*power(c,-(2*4+2*pNn(dk0,Omega,4))));
	x[0]=interval(-1.,1.)*(psi00+psi01+psi02+psi03+psi04);
	x[1]=x[0];
	return x;
}

IVector T1(int dk0,IVector b,IVector psi,IVector z)//T_1^Omega for psi00+psi01+psi02+psi03+psi04 at point z (Lemma 5.4.4)
{
	int Omega=omega(dk0);
	IVector x(2),y(2),psi01(2),psi02(2),psi03(2),psi04(2);
	psi01=b[0]*power(z,-2);
	psi02=b[1]*power(z,-4);
	psi03=b[2]*power(z,-6);
	psi04=b[3]*power(z,-8);
	x=nu0(dk0)*power(z,-2);
	y=dk0*psi+pow(dk0,2)*(cprod(psi01,(psi04+psi))+cprod(psi02,(psi03+psi04+psi))+power(psi03+psi04+psi,2)/2);
	y=y+pow(dk0,3)*(cprod(power(psi01,2),(psi03+psi04+psi))/2+cprod(psi01,power(psi02+psi03+psi04+psi,2))/2+power(psi02+psi03+psi04+psi,3)/fac(3));
	y=y+pow(dk0,4)*(cprod(power(psi01,3),(psi02+psi03+psi04+psi))/interval(6)+cprod(power(psi01,2),power(psi02+psi03+psi04+psi,2))/interval(4)+cprod(psi01,power(psi02+psi03+psi04+psi,3))/interval(6)+power(psi02+psi03+psi04+psi,4)/fac(4));
	for(int k=5;k<=Omega;k++)
	{
		y=y+pow(dk0,k)*power(psi01+psi02+psi03+psi04+psi,k)/fac(k);
	}
	x=cprod(x,y);
	return x;
}

IVector T1b(int dk0,IVector b,IVector psi,IVector z)//T_1^{k,b} for psi00+psi01+psi02+psi03+psi04 at point z (Lemma 5.4.4)
{
	int Omega=omega(dk0);
	interval c=mod(z);
	IVector x(2);
	c=CN(dk0,b,0,c)+mod(psi);
	c=exp(dk0*c)*power(dk0*c,Omega+1);
	c=c*(nu0(dk0)/(fac(Omega+1)*power(mod(z),2)));
	x[0]=interval(-1.,1.)*c;
	x[1]=x[0];
	return x;
}

IVector F(int dk0,IVector b,IVector psi,IVector z)//Right hand side of the iterative method for psi_1^u and psi_1^s
{
	return T0(dk0,b,z)+T1(dk0,b,psi,z)+T0b(dk0,b,z)+T1b(dk0,b,psi,z);
}

IVector dF(int dk0,IVector b,IVector psi,IVector z)//Right hand side of the iterative method for the derivative of psi_1^s
{
	return nu0(dk0)*dk0*exp(dk0*(psi0(dk0,b,z)+psi));
}

// We treat Fu and Fs separately only for the sake of documentation.
// The function is the same, but the justification of the formulae 
// is written out differently.

IVector Fu(int dk0,IVector b,IVector x,IVector z)//Iterative function for psi_1^u
{
	// below the convention is
	// d=d(n)=psi(n)-psi(n-1)
	// psi=psi(n)
	IVector d(2),psi(2);
	split(x,d,psi);
	// below we comute d(n+1):
	d = F(dk0,b,psi,z)+d;
	// Note that 
	// d(n+1) = psi(n+1)-psi(n)
	// psi(n+1) = psi(n)+d(n+1):
	psi = psi+d;
	return combine(d,psi);
}

IVector Fs(int dk0,IVector b,IVector x,IVector z)//Iterative function for psi_1^s
{
	// below the convention is:
	// d=d(n)=psi(n)-psi(n+1)
	// psi=psi(n)
	IVector d(2),psi(2);
	split(x,d,psi);
	// below we compute d(n-1):
	d = F(dk0,b,psi,z)+d;
	// Note that
	// d(n-1) = psi(n-1) - psi(n)
	// psi(n-1) = psi(n) + d(n-1)
	psi = psi + d;
	return combine(d,psi);
}

IVector dFs(int dk0,IVector b,IVector dx,IVector z)//Iterative function for the derivative of psi_1^s
{
	// below the convention is:
	// d=d(n)=psi(n)-psi(n+1)
	// psi=psi(n)
	// dd=dd(n)=dpsi(n)-dpsi(n+1)
	// dpsi=dpsi(n)
	IVector d(2),psi(2),dd(2),dpsi(2);
	dsplit(dx,d,psi,dd,dpsi);
	// below we compute d(n-1) and dd(n-1):
	d = F(dk0,b,psi,z)+d;
	dd = cprod(dF(dk0,b,psi,z),dpsi)+dd;
	// Note that:
	// d(n-1) = psi(n-1) - psi(n) (same for dd(n-1)
	// psi(n-1) = psi(n) + d(n-1) (same for dpsi(n-1)
	psi = psi + d;
	dpsi = dpsi + dd;
	return dcombine(d,psi,dd,dpsi);
}

//Initial conditions for the iterative methods (Sections 51 and 52)

IVector psi10(int dk0,IVector b,interval n,interval rho,interval gamma)//Initial condition for psi_1^{u,s} (using Proposition 5.1.2)
{
	interval varrho=rho+2*gamma;
	IVector z({n,-varrho}),psi10(2);
	psi10[0]=interval(-1.,1.)*M0(dk0,b,rho,gamma)/(power(mod(z),10));
	psi10[1]=psi10[0];
	return psi10;
}

IVector d0_u(int dk0,IVector b,interval n,interval rho,interval gamma)//Initial condition for d^u
{
	return psi10(dk0,b,n,rho,gamma) - psi10(dk0,b,n-1,rho,gamma);
}

IVector d0_s(int dk0,IVector b,interval n,interval rho,interval gamma)//Initial condition for d^s 
{
	return psi10(dk0,b,n,rho,gamma) - psi10(dk0,b,n+interval(1.),rho,gamma);
}

IVector x0_u(int dk0,IVector b,interval n,interval rho,interval gamma)// Initial conditions unstable combined
{
	return combine(d0_u(dk0,b,n,rho,gamma),psi10(dk0,b,n,rho,gamma));
}

IVector x0_s(int dk0,IVector b,interval n,interval rho,interval gamma)// Initial conditions stable combined
{
	return combine(d0_s(dk0,b,n,rho,gamma),psi10(dk0,b,n,rho,gamma));
}

IVector dpsi10(int dk0,IVector b,interval n,interval rho,interval gamma)//Initial condition for the derivative of psi_1^s (see Remark 5.3.4)
{
	interval M1,varrho=rho+2*gamma;
	IVector z({n,-varrho}),dpsi10(2);
	M1=interval(2)*M0(dk0,b,rho,gamma)/c(gamma);
	M1=M1+(10*M0(dk0,b,rho+gamma/2,gamma))/(rho+3*gamma/2);//Here we use that we are computing the derivative at im z= -i(rho+2gamma)
	M1=M1/(power(mod(z),10));
	dpsi10[0]=interval(-1.,1.)*M1;
	dpsi10[1]=dpsi10[0];
	return dpsi0(dk0,b,z)+dpsi10;
}

IVector dd0_s(int dk0,IVector b,interval n,interval rho,interval gamma)//Initial condition for derivative d^s 
{
	return dpsi10(dk0,b,n,rho,gamma) - dpsi10(dk0,b,n+interval(1.),rho,gamma);
}

IVector dx0_s(int dk0,IVector b,interval n,interval rho,interval gamma)// Initial conditions derivative stable combined
{
	return dcombine(d0_s(dk0,b,n,rho,gamma),psi10(dk0,b,n,rho,gamma),dd0_s(dk0,b,n,rho,gamma),dpsi10(dk0,b,n,rho,gamma));//Using Cauchy estimates with r=1/2
}

//<<<<<< Definition of the Simpson sum and its error (Section 2.2) >>>>>>

IVector Wu(int dk0,IVector b,interval re,int n,interval rho,interval gamma)// Psi1 unstable at points re-1-irho and re-irho interated from -n
{
	interval varrho=-rho-2*gamma;//Height where to compute the difference (where we have defined the interval re(z) in [0,1])
	IVector z(2),x=x0_u(dk0,b,interval(n+re),rho,gamma),psi1(2),psi2(2),d(2);
	z[0]=interval(n+re);
	z[1]=varrho;
	for(int i=0;i<abs(n)-1;i++)
	{
		x=Fu(dk0,b,x,z);
		z[0]=z[0]+1;
	}
	split(x,d,psi1);
	x=Fu(dk0,b,x,z);
	split(x,d,psi2);
	x=combine(psi1,psi2);
	return x;
}

IVector dWs(int dk0,IVector b,interval re,int n,interval rho,interval gamma)// psi1, dpsi1 stable at points re-irho and re+1-irho interated from n
{
	interval varrho=-rho-2*gamma;//Height where to compute the difference (where we have defined the interval re(z) in [0,1])
	IVector z(2),x=dx0_s(dk0,b,interval(n+re),rho,gamma),psi1(2),psi2(2),d(2),dpsi1(2),dpsi2(2),dd(2);
	z[0]=interval(n+re);
	z[1]=varrho;
	for(int i=0;i<n-1;i++)
	{
		x=dFs(dk0,b,x,z);
		z[0]=z[0]-1;
	}
	dsplit(x,d,psi1,dd,dpsi1);
	x=dFs(dk0,b,x,z);
	dsplit(x,d,psi2,dd,dpsi2);
	x=dcombine(psi1,psi2,dpsi1,dpsi2);
	return x;
}

//<<<<<< Definition of the Simpson sum and its error (Section 5.1.2) >>>>>>

IVector SSum(int dk0,IVector b,int L,interval rho,interval gamma,int m)// Simpson sum of p2(z)*e^(2pi*z) between -ivarrho and 1-ivarrho
{
	interval s,varrho;
	IVector err(2),SS(2),x(4),y(2),dx(8),d(2),dd(2),psi1sL(2),psi1uL(2),psi1sR(2),psi1uR(2),dpsisL(2),dpsisR(2),z(2),deltaL(2),deltaR(2);
	varrho=rho+2*gamma;
	SS[0]=interval(0.);
	SS[1]=interval(0.);
	for(int i=0;i<m+1;i++)
	{
		s=interval(i)/interval(m);
		x=Wu(dk0,b,s,-L,rho,gamma);//Computing psi_1^u at -ivarrho+s-1 and -ivarrho+s
		split(x,psi1uL,psi1uR);
		dx=dWs(dk0,b,s-1,L,rho,gamma);//Computing psi_1^s and its derivative at -ivarrho+s-1 and -ivarrho+s
		dsplit(dx,psi1sR,psi1sL,dpsisR,dpsisL);
		deltaL=psi1uL-psi1sL;
		deltaR=psi1uR-psi1sR;
		z[1]=-varrho;
		z[0]=interval(s-1);
		err=power(mod(z),-11)*zeta11(dk0,b,rho,gamma)*interval(-1.,1.);
		y=cprod(deltaR,dpsisL+err);
		z[0]=interval(s);
		err=power(mod(z),-11)*zeta11(dk0,b,rho,gamma)*interval(-1.,1.);
		y=y-cprod(deltaL,dpsisR+err);
		z[0]=interval(0.);
		z[1]=2*interval::pi()*s;
		y=cprod(y,exp(z));
		if (i==0 || i==m){
			y=y;
		}
		else{
			if (pow(-1,i)==1){
				y=2*y;
			}
			else{
				y=interval(4)*y;
			}
		}
		SS=SS+y;
	}
	SS=SS/interval(3*m);
	return SS;
}

interval ESS(int dk0,IVector b,interval rho,interval gamma,int P)//Error for the Simpson sum (Lemma 5.1.6)
{
	interval r,a,x,y;
	a=c(gamma);
	r=rho+2*gamma-a/2;
	y=interval(0.);
	for (int i=0;i<5;i++)
	{
		y=y+(fac(4)/fac(i))*power(interval::pi(),i)*power(a,i-4);
	}
	x=eta1(dk0,b,rho,gamma)/r;
	x=x+(zeta11(dk0,b,rho,gamma)/power(r,11));
	x=(x*M0(dk0,b,rho+gamma/2,gamma))/power(r,10);
	x=(x*pow(2,4)*y)/(45*pow(P,4));
	return x;
}

//<<<<<< Computer assited proof >>>>>>

// Remark 1: The method is really slow if P increases since it has to compute two times de difference for each step (see function SSum). My recomendation is to check 
//           first the sum with P small (say 10) and see that we obtain a good interval for the parameters chosen (observe that increasing the L also makes the process slower).
//           Then, once we obtain a good interval, repite the process including the error ESS with a good P for it  
// Remark 2: Observe that increasing L we start getting worst intervals, this is due to the number of computations and the accumulated error.
// Remark 3: There is an optimal rho that gives better (smaller) intervals as a result and it is not always the smaller or the bigger since
//           in both cases (increasing or decreasing it) we can create a different kind of error. On one side, if we decrease rho,
//           the error comming from zeta_1^{1,b} in the SSum increase. On the other side, increasing rho the number of correct digids in the 
//           computation of delta decrease.
// Summary:  It is recomended to play with the parameters in order to get a better knowledge of their behaviour and implrove the obtained results
//           The results given in the thesis are NOT optimal and CAN be improved.

void Stokes()
{
	int dk0=1,P=20000,L=1000;
	interval rho=interval(5.475),gamma=interval(0.949327),err,varrho,nu0;//
	IVector SS,b(4);
	nu0=interval(2.)/interval(dk0);
	varrho=interval(rho)+2*interval(gamma);
	b[0]=-(nu0/8);
	b[1]=(nu0*7*13)/(pow(2,6)*pow(3,3));
	b[2]=(-nu0*11*29)/(pow(2,7)*pow(3,2)*5);
	b[3]=(nu0*2837*9467)/(pow(2,11)*pow(3,6)*pow(5,2)*7);
	cout << "-- Trigonometric --" << endl;
	cout << "d_{k_0} = " << dk0 << endl;
	cout << "P = " << P << endl;
	cout << "L = " << L << endl;
	cout << "rho = " << rho << endl;
	cout << "varrho = " << rho+2*gamma << endl;
	cout << "Lipschidz = " << Lip(dk0,b,interval(rho),interval(gamma)) << " <=0.5 " << endl;//It has to be less than 1/2 to have existence of psi_1 (see Section 5.2.3)
	cout << "gamma = " << gamma << endl;
	cout << "c(gamma) = " << c(interval(gamma)) << " <=0.688493 " << endl;//Condition on gamma (see Proposition 5.3.2).
	cout << "Lipschidz^d = " << Lipd(dk0,b,interval(rho),interval(gamma)) << " < 1 " << endl;//It has to be less than 1 to have existence of zeta_1 (see proof of Lemma 5.3.8).
	SS=SSum(dk0,b,L,interval(rho),interval(gamma),P);//Simpson sum without error. 
	cout << "Simpson sum  = (" << SS*exp(2*interval::pi()*varrho) <<")" << endl;
	err=ESS(dk0,b,rho,gamma,P);//Error Simpson sum. Here we see which P we need.
	cout << "ESS " << err*exp(2*interval::pi()*varrho) << endl;
	cout << "Stokes constant = " << (SS[1]+err*interval(-1.,1.))*exp(2*interval::pi()*varrho) << endl;
}


int main()
{
  	cout.precision(10);
	try
	{
		std::chrono::time_point<std::chrono::system_clock> start, end;
		start = std::chrono::system_clock::now();
		Stokes();
		end = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds = end - start;
		std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
	}
	catch(exception& e)
  	{
    		cout << "\n\nException caught: "<< e.what() << endl;
  	}
  return 0;
} 
