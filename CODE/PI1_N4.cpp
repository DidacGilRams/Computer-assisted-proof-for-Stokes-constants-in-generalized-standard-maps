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

//----------------------------------- Polynomial case with one monomial ----------------------------------------------------------

//This code gives an interval for the Stokes constant p_2^{-1}(d_{k_0}) for every degree d_{k_0} (has to be introduced in the test() function) when I={d_{k_0}}.

//<<<<<< Implementation of some functions used in the program >>>>>>

interval angle(interval x,interval y)// Function that computes the angle on the plane
{
    if((y>=0) && (2*x>=y))           return  atan(y/x);
	if((y>=0) && (abs(x)<=2*abs(y))) return -atan(x/y)+interval(0.5)*interval::pi();
	if((x<=0) && (2*abs(x)>=abs(y))) return  atan(y/x)-interval::pi();
	if((y<=0) && (abs(x)<=2*abs(y))) return -atan(x/y)-interval(0.5)*interval::pi();
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
	if(a>interval::pi()) return mod2Pi(a-interval(2.)*interval::pi());
	if(a<-interval::pi()) return mod2Pi(a+interval(2.)*interval::pi());
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
	interval j=interval(interval(1.));
	for(int i=1;i<=n;i++)
	{
		j=j*interval(i);
	}
	return j;
}

interval dfac(int n)// Double factorial
{
	interval j=interval(n);
	for(int i=1;2*i<n;i++)
	{
		j=j*interval(n-2*i);
	}
	return j;
}

interval binom(interval a,int k)//Generalized binomial
{
	if (k==0)
	{
		return interval(1);
	}
	else
	{
		interval j=a,x=a/interval(k);
		for(int i=1;i<k;i++)
		{
			x=x*((a-i)/(k-i));
		}
		return x;
	}
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

//<<<<<< Obtaining rho_0 for the existence of solutions of the inner equation (Section 3) >>>>>>

//Norms of the inverse operators

interval nu_0(int dk0)//Definition of nu_0 in the trigonometric case
{
	return interval(2.)/interval(dk0-1); 
}

int pNn(int dk0,int N,int n)
{
	interval x;
	x=interval(N-n)/interval(dk0-1);
	return floor(x)+1;
}


interval B(interval nu)//Defines the constants B(nu) from bounding the inverses S (Proposition 3.5)
{
	int n=floor(nu);
	interval x,Pi=interval::pi();
	x=interval(dfac(n-2)/dfac(n-1));
	if (power(-1,n)==1){
		x=x;
	}
	else{
		x=x*Pi/interval(2.);
	}
	return x;
}

interval S(interval rho,interval gamma,interval nu)// Bound of the first order inverse operator in a domain with rho aplied to a nu Banach space (Proposition 3.5)
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

interval G(int dk0,interval rho,interval gamma)// Bound of the second order inverse operator in a domain with rho aplied to a nu Banach space (Proposition 3.5)
//Here for the polynomial case kappa=interval(2.)/(dk0-1) and since we took psi00+psi01+psi02+psi03+psi04 as approximation nu_N=10+nu_0
{
	interval x;
	x=S(rho,gamma,interval(3)+2*nu_0(dk0));
	x=x*(S(rho,gamma,interval(9))+S(rho,gamma,interval(12)+2*nu_0(dk0)));
	return x;
}

//Here we define the coefficients b_i(d_{k_0}) for i={1,2,3,4} of the functions psi_0^i

interval b1(int dk0)//Coefficient of psi_0^1
{
	interval x,nu0=nu_0(dk0);
	x=-nu0*(nu0+interval(1.))*(nu0+interval(3.));
	return x/fac(4);
}

interval b2(int dk0)//Coefficient of psi_0^2
{
	interval x,nu0=nu_0(dk0);
	x=(nu0*(nu0+interval(1.))*(nu0+interval(2.))*(nu0+interval(3.))*(nu0+interval(4.))*(nu0+interval(5.)))/fac(6);
	x=x+((nu0+interval(2.))*(nu0+interval(3.))*(nu0+interval(4.))*(nu0+interval(5.))*b1(dk0))/fac(4);
	x=-interval(2.)*x;
	x=x+dk0*(nu0+interval(1.))*power(b1(dk0),2);
	x=x/(((nu0+interval(4.))*(nu0+interval(5.)))-((nu0+interval(1.))*(nu0+interval(2.))));
	return x;
}

interval b3(int dk0)//Coefficient of psi_0^3
{
	interval x,nu0=nu_0(dk0);
	x=(nu0*(nu0+interval(1.))*(nu0+interval(2.))*(nu0+interval(3.))*(nu0+interval(4.))*(nu0+interval(5.))*(nu0+interval(6.))*(nu0+interval(7.)))/fac(8);
	x=x+((nu0+interval(2.))*(nu0+interval(3.))*(nu0+interval(4.))*(nu0+interval(5.))*(nu0+interval(6.))*(nu0+interval(7.))*b1(dk0))/fac(6);
	x=x+((nu0+interval(4.))*(nu0+interval(5.))*(nu0+interval(6.))*(nu0+interval(7.))*b2(dk0))/fac(4);
	x=-interval(2.)*x;
	x=x+nu0*(nu0+interval(1.))*dk0*(dk0-interval(1.))*(b1(dk0)*b2(dk0)+((dk0-interval(2.))*power(b1(dk0),3))/fac(3));
	x=x/(((nu0+interval(6.))*(nu0+interval(7.)))-((nu0+1)*(nu0+2)));
	return x;
}

interval b4(int dk0)//Coefficient of psi_0^4
{
	interval x,nu0=nu_0(dk0);
	x=(nu0*(nu0+interval(1.))*(nu0+interval(2.))*(nu0+interval(3.))*(nu0+interval(4.))*(nu0+interval(5.))*(nu0+interval(6.))*(nu0+interval(7.))*(nu0+interval(8.))*(nu0+interval(9.)))/fac(10);
	x=x+((nu0+interval(2.))*(nu0+interval(3.))*(nu0+interval(4.))*(nu0+interval(5.))*(nu0+interval(6.))*(nu0+interval(7.))*(nu0+interval(8.))*(nu0+interval(9.))*b1(dk0))/fac(8);
	x=x+((nu0+interval(4.))*(nu0+interval(5.))*(nu0+interval(6.))*(nu0+interval(7.))*(nu0+interval(8.))*(nu0+interval(9.))*b2(dk0))/fac(6);
	x=x+((nu0+interval(6.))*(nu0+interval(7.))*(nu0+interval(8.))*(nu0+interval(9.))*b3(dk0))/fac(4);
	x=-interval(2.)*x;
	x=x+nu0*(nu0+interval(1.))*dk0*(dk0-interval(1.))*(b1(dk0)*b3(dk0)+power(b2(dk0),2)/interval(2.)+(interval(3.)*(dk0-interval(2.))*power(b1(dk0),2)*b2(dk0))/fac(3)+((dk0-interval(2.))*(dk0-interval(3.))*power(b1(dk0),4))/fac(4));
	x=x/(((nu0+interval(8.))*(nu0+interval(9.)))-((nu0+1)*(nu0+2)));
	return x;
}


//The following fuctions are used to compute the norm of the first iteration. 

interval E0(int dk0,IVector b,interval rho)//E0 for psi00+psi01+psi02+psi03+psi04 (Lemma 3.1)
{
	//Basically we compute abs(bi)*abs(delta^2 [psi0i/bi]-2sum cancelled terms) for each psi_0^i.
	interval nu0=nu_0(dk0);
	interval psi00,psi01,psi02,psi03,psi04;
	psi00=interval(1.)/power(rho-interval(1.),nu0)+interval(1.)/power(rho+interval(1.),nu0)-interval(2.)/power(rho,nu0);
	psi00=psi00*power(rho,nu0+interval(10.));
	for(int i=1;i<6;i++)
	{
		psi00=psi00-interval(2.)*binom(nu0+interval(2.)*i-interval(1.),2*i)*power(rho,10-(2*i));
	}
	psi00=abs(psi00);
	psi01=interval(1.)/power(rho-interval(1.),nu0+interval(2.))+interval(1.)/power(rho+interval(1.),nu0+interval(2.))-interval(2.)/power(rho,nu0+interval(2.));
	psi01=psi01*power(rho,nu0+interval(10.));
	for(int i=1;i<5;i++)
	{
		psi01=psi01-interval(2.)*binom(nu0+interval(2.)*(i+1)-interval(1.),2*i)*power(rho,10-2*(i+1));
	}
	psi01=abs(b[0])*abs(psi01);
	psi02=interval(1.)/power(rho-1,nu0+interval(4.))+interval(1.)/power(rho+1,nu0+interval(4.))-interval(2.)/power(rho,nu0+interval(4.));
	psi02=psi02*power(rho,nu0+interval(10.));
	for(int i=1;i<4;i++)
	{
		psi02=psi02-interval(2.)*binom(nu0+interval(2.)*(i+2)-interval(1.),2*i)*power(rho,10-2*(i+2));
	}
	psi02=abs(b[1])*abs(psi02);
	psi03=interval(1.)/power(rho-1,nu0+interval(6.))+interval(1.)/power(rho+1,nu0+interval(6.))-interval(2.)/power(rho,nu0+interval(6.));
	psi03=psi03*power(rho,nu0+interval(10.));
	for(int i=1;i<3;i++)
	{
		psi03=psi03-interval(2.)*binom(nu0+interval(2.)*(i+3)-interval(1.),2*i)*power(rho,10-2*(i+3));
	}
	psi03=abs(b[2])*abs(psi03);
	psi04=interval(1.)/power(rho-1,nu0+interval(8.))+interval(1.)/power(rho+1,nu0+interval(8.))-interval(2.)/power(rho,nu0+interval(8.));
	psi04=psi04*power(rho,nu0+interval(10.));
	for(int i=1;i<2;i++)
	{
		psi04=psi04-interval(2.)*binom(nu0+interval(2.)*(i+4)-interval(1.),2*i)*power(rho,10-2*(i+4));
	}
	psi04=abs(b[3])*abs(psi04);
	return psi00+psi01+psi02+psi03+psi04;
}

interval CN(int dk0,IVector b,int a,interval rho)//Funtion C_N(a,rho)
{
	interval c=interval(0.);
	for(int n=1;n<=4;n++)
	{
		c=c+abs(b[n-1])*power(rho,nu_0(dk0)*a-2*n);//Absolute values of psi_0^n evaluated at rho
	}
	return c;
}

interval f0(int dk0,IVector b,interval rho)//Computing the bound of the first iteration that comes from f - Taylor cancelled terms for psi00+psi01+psi02+psi03+psi04 (Lemma 3.2)
{
	interval nu0=nu_0(dk0);
	interval x,y,psi00,psi01,psi02,psi03,psi04;
	psi00=power(rho,-nu0);//Absolute values of psi0i evaluated at rho
	psi01=abs(b[0])*power(rho,-nu0-interval(2.));
	psi02=abs(b[1])*power(rho,-nu0-interval(4.));
	psi03=abs(b[2])*power(rho,-nu0-interval(6.));
	psi04=abs(b[3])*power(rho,-nu0-interval(8.));
	y=interval(2.)*(psi01*psi04)+interval(2.)*(psi02*(psi03+psi04))+power(psi03+psi04,2);//Function f(psi1) evaluated at 0
	x=binom(dk0,2)*(power(psi00,dk0-2)*y);
	y=interval(3.)*(power(psi01,2)*(psi03+psi04))+interval(3.)*(psi01*power(psi02+psi03+psi04,2))+power(psi02+psi03+psi04,3);
	x=x+binom(dk0,3)*(power(psi00,dk0-3)*y);
	y=interval(4.)*(power(psi01,3)*(psi02+psi03+psi04))+interval(6.)*(power(psi01,2)*power(psi02+psi03+psi04,2))+interval(4.)*(psi01*power(psi02+psi03+psi04,3))+power(psi02+psi03+psi04,4);
	x=x+binom(dk0,4)*(power(psi00,dk0-4)*y);
	for(int i=5;i<=dk0;i++)
	{
		x=x+binom(dk0,i)*power(psi00,dk0-i)*power(psi01+psi02+psi03+psi04,i);
	}
	return nu0*(nu0+1)*x*power(rho,nu0+interval(10.));
}

interval F0(int dk0,IVector b,interval rho,interval gamma)//Gives the first iteration bound (Corollary 3.4 and Lemma 3.7)
{
	return G(dk0,rho,gamma)*(E0(dk0,b,rho)+f0(dk0,b,rho));
}

//The following are used to compute the Lipschidz constant

interval A(int dk0,IVector b,interval rho)//Boundig the norm of the linear term A (Lemma 3.6)
{
	interval c=interval(0.),psi,nu0=nu_0(dk0);
	psi=CN(dk0,b,0,rho);
	for (int i=1;i<dk0;i++)
	{
		c=c+binom(dk0-interval(1.),i)*power(psi,dk0-i);
	}
	c=(nu0+2)*(nu0+1)*c;
	c=c+power(rho,2)*(power(interval(1.)-power(rho,-1),nu0+interval(2.))+power(interval(1.)+power(rho,-1),nu0+interval(2.))-interval(2.))-interval(2.)*binom(nu0+interval(2.),2);
	return abs(c);
}

interval R(int dk0,IVector b,interval rho,interval M0)//Bounding higher order terms (Lemma 3.8)
{
	interval psi,nu0=nu_0(dk0);
	psi=CN(dk0,b,0,rho);
	return (nu0+2)*(nu0+1)*(power(1+psi+M0/power(rho,interval(10.)),dk0-1)-power(1+psi,dk0-1));//multiplied by rho^2 from the inverse G
}

interval M0(int dk0,IVector b,interval rho,interval gamma)//Function M_0(rho,gamma)
{
	return 2*F0(dk0,b,rho,gamma);
}

interval Lip(int dk0,IVector b,interval rho,interval gamma)//Funtion that computes the Lipschidz constant for the existence (Lemma 3.8)
{
	return G(dk0,rho,gamma)*(A(dk0,b,rho)+R(dk0,b,rho,M0(dk0,b,rho,gamma)));
}

//<<<<<< Obtaining rho_1 for the existence of zeta_1 and explicit formula for psi^u-psi^s (Section 4) >>>>>>

//Norms of the inverse operators

interval c(interval gamma)//Definition of c(gamma) involved in the domains (see definition of D_{rho,gamma}^delta)
{
	interval x;
	x=power(gamma,2)+1;
	x=power(x,interval(0.5));
	x=gamma/x;
	return x;
}

interval Sd(interval rho,interval gamma,interval nu)// Bound of the first order inverse operator in a c-ascending domain aplied to a nu Banach space (Proposition interval(4)3)
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

interval eta1(int dk0,IVector b,interval rho,interval gamma)// Lemma 4.3
{
	interval r,x,nu0=nu_0(dk0);
	r=(rho+3*gamma/2);
	x=nu0;
	for(int n=1;n<=4;n++)
	{
		x=x+(2*n+nu0)*abs(b[n-1])*power(rho,-2*n);//Absolute values of dpsi_0^n evaluated at rho
	}
	x=x+(2*M0(dk0,b,rho,gamma))/(c(gamma)*power(rho,9));
	return x+((1+nu0)*M0(dk0,b,rho+gamma/2,gamma))/(power(r,10));
}

interval eta2(int dk0,IVector b,interval rho,interval gamma)//Corollary 4.4
{
	interval r,x,eta,nu0=nu_0(dk0);
	r=(rho+3*gamma/2);
	eta=eta1(dk0,b,rho,gamma);
	x=1+power(r,-1);
	x=Sd(rho,gamma,-interval(2*nu0+3))*eta*power((2*nu0)-eta,-2)*power(x,nu0+1);
	return x;
}

interval Gd(int dk0,IVector b,interval rho,interval gamma)// Bound of the second order inverse operator in a domain with rho aplied to a nu Banach space (Corollary interval(4.)5)
{
	interval nu0=nu_0(dk0);
	return eta1(dk0,b,rho,gamma)*eta2(dk0,b,rho,gamma)*(Sd(rho,gamma,10)+Sd(rho,gamma,13+2*nu0));
}

interval Ad(int dk0,IVector b,interval rho, interval gamma)//Funtion involved in the bound for the linear operator A^d (Lemma interval(4.)1)
{
	interval nu0=nu_0(dk0);
	interval c,r,M;
	r=rho+interval(3.)*gamma/interval(2.);
	M=M0(dk0,b,rho+gamma/interval(2.),gamma);
	c=interval(4.)*M*(nu0+1)*dk0*power(1+CN(dk0,b,0,r)+(M/power(r,10)),dk0-2);
	return c;
}

interval Lipd(int dk0,IVector b,interval rho,interval gamma)//Bound for the operator G^d[A^d*] (Lemma 4.7)
{
	return Gd(dk0,b,rho,gamma)*Ad(dk0,b,rho,gamma)*power(rho+interval(3.)*gamma/interval(2.),-10);
}

interval zeta11(int dk0,IVector b,interval rho,interval gamma)//Bound for zeta_1^1 (Lemma 4.7)
{
	return (eta1(dk0,b,rho,gamma)*Gd(dk0,b,rho,gamma)*Ad(dk0,b,rho,gamma))/(interval(1.)-Lipd(dk0,b,rho,gamma));
}

//<<<<<< Rigorous computations of psi^u-psi^s and dpsi^s (Section 5) >>>>>>

//Initial aproximation of the manifold and its derivative 

IVector psi0(int dk0,IVector b,IVector z)// psi0=psi00+psi01+psi02+psi03+psi04
{
	interval nu0=nu_0(dk0);
	IVector x(2);
	x=power(z,-nu0);
	for(int n=1;n<=4;n++)
	{
		x=x+b[n-1]*power(z,-2*n-nu0);//psi_0^n evaluated at z
	}
	return x;
}

IVector psi0(int dk0,IVector b,interval re,interval im)// psi0=psi00+psi01+psi02+psi03+psi04
{
	IVector z(2);
	z[0]=re;
	z[1]=im;
	return psi0(dk0,b,z);
}

IVector dpsi0(int dk0,IVector b,IVector z)// Derivative of psi0=psi00+psi01+psi02+psi03+psi04
{
	interval nu0=nu_0(dk0);
	IVector x(2);
	x=-nu0*power(z,-1-nu0);
	for(int n=1;n<=4;n++)
	{
		x=x-(2*n+nu0)*b[n-1]*power(z,-(2*n+1)-nu0);//dpsi_0^n evaluated at z
	}
	return x;
}

IVector dpsi0(int dk0,IVector b,interval re,interval im)// Derivative of psi0=psi00+psi01+psi02+psi03+psi04
{
	IVector z(2);
	z[0]=interval(re);
	z[1]=im;
	return dpsi0(dk0,b,z);
}

//<<<<<< Iterative methods for psi^u, psi^s (Section interval(5.)1) and derivative psi^s (Section interval(5.)2) >>>>>>

//The following functions are defined to obtain the generating functions F and dF for the iterative methods

int omega(int dk0)
{
		return 10*dk0;
}

IVector T0(int dk0,IVector b,IVector z)//T0 for psi00+psi01+psi02+psi03+psi04 at point z (Lemma interval(5.)1)
{
	int Omega=omega(dk0);
	interval nu0=nu_0(dk0);
	IVector psi00(2),psi01(2),psi02(2),psi03(2),psi04(2);
	psi00[0]=interval(0.);
	psi00[1]=psi00[0];
	psi01=psi00;
	psi02=psi00;
	psi03=psi00;
	psi04=psi00;
	for(int i=6;i<=pNn(dk0,Omega,0);i++)
	{
		psi00=psi00+binom(nu0+interval(2.)*i-1,2*i)*power(z,-nu0-interval(2.)*i);
	}
	for(int i=5;i<=pNn(dk0,Omega,1);i++)
	{
		psi01=psi01+binom(nu0+interval(2.)*(i+1)-interval(1.),2*i)*power(z,-nu0-interval(2.)*(i+1));
	}
	psi01=b[0]*psi01;
	for(int i=4;i<=pNn(dk0,Omega,2);i++)
	{
		psi02=psi02+binom(nu0+interval(2.)*(i+2)-interval(1.),2*i)*power(z,-nu0-interval(2.)*(i+2));
	}
	psi02=b[1]*psi02;
	for(int i=3;i<=pNn(dk0,Omega,3);i++)
	{
		psi03=psi03+binom(nu0+interval(2.)*(i+3)-interval(1.),2*i)*power(z,-nu0-interval(2.)*(i+3));
	}
	psi03=b[2]*psi03;
	for(int i=2;i<=pNn(dk0,Omega,4);i++)
	{
		psi04=psi04+binom(nu0+interval(2.)*(i+4)-interval(1.),2*i)*power(z,-nu0-interval(2.)*(i+4));
	}
	psi04=b[3]*psi04;
	return -2*(psi00+psi01+psi02+psi03+psi04);
}

IVector T0b(int dk0,IVector b,IVector z)//Rest from T0 for psi00+psi01+psi02+psi03+psi04 at point z. (Lemma 51 and using Lemma 33)
{
	int Omega=omega(dk0);
	interval psi00,psi01,psi02,psi03,psi04,c,nu0=nu_0(dk0);
	IVector x(2);
	c=mod(z);
	psi00=nu0*((2*pNn(dk0,Omega,0)-1)/fac(2*(pNn(dk0,Omega,0)+1)))*(power(c-1,-(2*pNn(dk0,Omega,0)+nu0))+power(c+1,-(2*pNn(dk0,Omega,0)+nu0))-2*power(c,-(2*pNn(dk0,Omega,0)+nu0)));
	psi01=abs(b[0])*binom(2*pNn(dk0,Omega,1)+1+nu0,2*pNn(dk0,Omega,1))*(power(c-1,-(2+2*pNn(dk0,Omega,1)+nu0))+power(c+1,-(2+2*pNn(dk0,Omega,1)+nu0))-2*power(c,-(2+2*pNn(dk0,Omega,1)+nu0)));
	psi02=abs(b[1])*binom(2*2+2*pNn(dk0,Omega,2)-1+nu0,2*pNn(dk0,Omega,2))*(power(c-1,-(2*2+2*pNn(dk0,Omega,2)+nu0))+power(c+1,-(2*2+2*pNn(dk0,Omega,2)+nu0))-2*power(c,-(2*2+2*pNn(dk0,Omega,2)+nu0)));
	psi03=abs(b[2])*binom(2*3+2*pNn(dk0,Omega,3)-1+nu0,2*pNn(dk0,Omega,3))*(power(c-1,-(2*3+2*pNn(dk0,Omega,3)+nu0))+power(c+1,-(2*3+2*pNn(dk0,Omega,3)+nu0))-2*power(c,-(2*3+2*pNn(dk0,Omega,3)+nu0)));
	psi04=abs(b[3])*binom(2*4+2*pNn(dk0,Omega,4)-1+nu0,2*pNn(dk0,Omega,4))*(power(c-1,-(2*4+2*pNn(dk0,Omega,4)+nu0))+power(c+1,-(2*4+2*pNn(dk0,Omega,4)+nu0))-2*power(c,-(2*4+2*pNn(dk0,Omega,4)+nu0)));
	x[0]=interval(-1.,1.)*(psi00+psi01+psi02+psi03+psi04);
	x[1]=x[0];
	return x;
}

IVector T1(int dk0,IVector b,IVector psi,IVector z)//Taylor for psi00+psi01+psi02+psi03+psi04 at point z (Lemma interval(5.)1)
{
	interval nu0=nu_0(dk0);
	IVector x(2),y(2),psi00(2),psi01(2),psi02(2),psi03(2),psi04(2);
	psi00=power(z,-nu0);
	psi01=b[0]*power(z,-nu0-interval(2.));
	psi02=b[1]*power(z,-nu0-interval(4.));
	psi03=b[2]*power(z,-nu0-interval(6.));
	psi04=b[3]*power(z,-nu0-interval(8.));
	x=dk0*cprod(power(psi00,dk0-1),psi);
	y=interval(2.)*cprod(psi01,psi04+psi)+interval(2.)*cprod(psi02,psi03+psi04+psi)+power(psi03+psi04+psi,2);
	x=x+binom(dk0,2)*cprod(power(psi00,dk0-2),y);
	y=interval(3.)*cprod(power(psi01,2),psi03+psi04+psi)+interval(3.)*cprod(psi01,power(psi02+psi03+psi04+psi,2))+power(psi02+psi03+psi04+psi,3);
	x=x+binom(dk0,3)*cprod(power(psi00,dk0-3),y);
	y=interval(4.)*cprod(power(psi01,3),psi02+psi03+psi04+psi)+interval(6.)*cprod(power(psi01,2),power(psi02+psi03+psi04+psi,2))+interval(4.)*cprod(psi01,power(psi02+psi03+psi04+psi,3))+power(psi02+psi03+psi04+psi,4);
	x=x+binom(dk0,4)*cprod(power(psi00,dk0-4),y);
	for(int i=5;i<=dk0;i++)
	{
		x=x+binom(dk0,i)*cprod(power(psi00,dk0-i),power(psi01+psi02+psi03+psi04+psi,i));
	}
	return nu0*(nu0+1)*x;
}

IVector F(int dk0,IVector b,IVector psi,IVector z)// Funtion involved in inner equation used to compute Wu and Ws (see Section interval(5.)interval(1.)2)
{
	return T0(dk0,b,z)+T1(dk0,b,psi,z)+T0b(dk0,b,z);
}

IVector dF(int dk0,IVector b,IVector psi,IVector z)// Funtion involved in inner equation used to compute Wu and Ws (see Section interval(5.)2)
{
	interval nu0=nu_0(dk0);
	return (nu0+interval(1.))*(nu0+interval(2.))*power(psi0(dk0,b,z)+psi,dk0-1);
}

// I treat Fu and Fs separately only for the sake of documentation. The function is the same, but the justification of the formulae is written out differently.

IVector Fu(int dk0,IVector b,IVector x,IVector z)//Iterative function unstable
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

IVector Fs(int dk0,IVector b,IVector x,IVector z)//Iterative function stable
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

IVector dFs(int dk0,IVector b,IVector dx,IVector z)//Iterative function derivative stable
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

//Initial conditions for the iterative methods (Sections interval(5.)1 and interval(5.)2)

IVector psi10(int dk0,IVector b,interval n,interval rho,interval gamma)//Initial condition for psi1 contain the error given by the theorem of existence, the M0 is the bound given by the first iteration and the space is nu_N=10
{
	interval varrho=rho+2*gamma;
	IVector z({n,-varrho}),psi10(2);
	psi10[0]=interval(-1.,1.)*M0(dk0,b,rho,gamma)/(power(mod(z),10+nu_0(dk0)));
	psi10[1]=psi10[0];
	return psi10;
}

IVector d0_u(int dk0,IVector b,interval n,interval rho,interval gamma)//Initial condition for d^u
{
	return psi10(dk0,b,n,rho,gamma) - psi10(dk0,b,n-1,rho,gamma);
}

IVector d0_s(int dk0,IVector b,interval n,interval rho,interval gamma)//Initial condition for d^s 
{
	return psi10(dk0,b,n,rho,gamma) - psi10(dk0,b,n+1,rho,gamma);
}

IVector x0_u(int dk0,IVector b,interval n,interval rho,interval gamma)// Initial conditions unstable combined
{
	return combine(d0_u(dk0,b,n,rho,gamma),psi10(dk0,b,n,rho,gamma));
}

IVector x0_s(int dk0,IVector b,interval n,interval rho,interval gamma)// Initial conditions stable combined
{
	return combine(d0_s(dk0,b,n,rho,gamma),psi10(dk0,b,n,rho,gamma));
}

IVector dpsi10(int dk0,IVector b,interval n,interval rho,interval gamma)//Initial condition for derivative psi1 contain the error given by the theorem of existence, the space is (10+nu_0)
{
	interval M1,nu0=nu_0(dk0),varrho=rho+2*gamma;
	IVector z({n,-varrho}),dpsi10(2);
	M1=2*M0(dk0,b,rho,gamma)/c(gamma);
	M1=M1+((10+nu0)*M0(dk0,b,rho+gamma/2,gamma))/(rho+interval(3.)*gamma/interval(2.));//Here we use that we are computing the derivative at im z= -(rho+2gamma)
	M1=M1/(power(mod(z),10+nu0));
	dpsi10[0]=interval(-1.,1.)*M1;
	dpsi10[1]=dpsi10[0];
	return dpsi0(dk0,b,z)+dpsi10;
}

IVector dd0_s(int dk0,IVector b,interval n,interval rho,interval gamma)//Initial condition for derivative d^s 
{
	return dpsi10(dk0,b,n,rho,gamma) - dpsi10(dk0,b,n+1,rho,gamma);
}

IVector dx0_s(int dk0,IVector b,interval n,interval rho,interval gamma)// Initial conditions derivative stable combined
{
	return dcombine(d0_s(dk0,b,n,rho,gamma),psi10(dk0,b,n,rho,gamma),dd0_s(dk0,b,n,rho,gamma),dpsi10(dk0,b,n,rho,gamma));//Using Cauchy estimates with r=1/2
}

//Definition of the iterative methods

IVector Wu(int dk0,IVector b,interval re,int L,interval rho,interval gamma)// Psi1 unstable at points re-1-irho and re-irho interated from -n
{
	interval varrho=-rho-2*gamma;//Height where to compute the difference (where we have defined the interval re(z) in [0,1])
	IVector z(2),x=x0_u(dk0,b,interval(-L+re),rho,gamma),psi1(2),psi2(2),d(2);
	z[0]=interval(-L+re);
	z[1]=varrho;
	for(int i=0;i<L-1;i++)
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

IVector dWs(int dk0,IVector b,interval re,int L,interval rho,interval gamma)// psi1, dpsi1 stable at points re-irho and re+1-irho interated from n
{
	interval varrho=-rho-2*gamma;//Height where to compute the difference (where we have defined the interval re(z) in [0,1])
	IVector z(2),x=dx0_s(dk0,b,interval(L+re),rho,gamma),psi1(2),psi2(2),d(2),dpsi1(2),dpsi2(2),dd(2);
	z[0]=interval(L+re);
	z[1]=varrho;
	for(int i=0;i<L-1;i++)
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

//<<<<<< Definition of the Simpson sum and its error (Section interval(2.)2) >>>>>>

IVector SSum(int dk0,IVector b,int L,interval rho,interval gamma,int P)// Simpson sum of p2(z)*e^(2pi*z) between -irho and 1-irho
{
	interval s,varrho,nu0=nu_0(dk0);
	IVector err(2),SS(2),x(4),y(2),dx(8),d(2),dd(2),psi1sL(2),psi1uL(2),psi1sR(2),psi1uR(2),dpsisL(2),dpsisR(2),z(2),deltaL(2),deltaR(2);
	varrho=rho+2*gamma;
	SS[0]=interval(0.);
	SS[1]=interval(0.);
	for(int i=0;i<P+1;i++)
	{
		s=interval(i)/interval(P);
		x=Wu(dk0,b,s,L,rho,gamma);//Computing psi_1^u at -irho+s-1 and -irho+s
		split(x,psi1uL,psi1uR);
		dx=dWs(dk0,b,s-1,L,rho,gamma);//Computing psi_1^s and dpsi^s at -irho+s-1 and -irho+s
		dsplit(dx,psi1sR,psi1sL,dpsisR,dpsisL);
		deltaL=psi1uL-psi1sL;
		deltaR=psi1uR-psi1sR;
		z[1]=-varrho;
		z[0]=interval(s-1);
		err=power(mod(z),-11-nu0)*zeta11(dk0,b,rho,gamma)*interval(-1.,1.);//Observe that this error decreses when considering a bigger rho, which means smaller interval in the end
		y=cprod(deltaR,dpsisL+err);
		z[0]=interval(s);
		err=power(mod(z),-11-nu0)*zeta11(dk0,b,rho,gamma)*interval(-1.,1.);
		y=y-cprod(deltaL,dpsisR+err);
		z[0]=interval(0.);
		z[1]=2*interval::pi()*s;
		y=cprod(y,exp(z));
		if (i==0 || i==P){
			y=y;
		}
		else{
			if (pow(-1,i)==1){
				y=2*y;
			}
			else{
				y=4*y;
			}
		}
		SS=SS+y;
	}
	SS=SS/interval(3*P);
	return SS;
}

interval ESS(int dk0,IVector b,interval rho,interval gamma,int P)//Error for the Simpson sum (Lemma 2.3)
{
	interval r,a,x,y,nu0=nu_0(dk0);
	a=c(gamma);
	r=rho+2*gamma-a/2;
	y=interval(0.);
	for (int i=0;i<5;i++)
	{
		y=y+(fac(4)/fac(i))*power(interval::pi(),i)*power(a,i-4);
	}
	x=eta1(dk0,b,rho,gamma)/power(r,nu0+1);
	x=x+(zeta11(dk0,b,rho,gamma)/power(r,11+nu0));
	x=(x*M0(dk0,b,rho+gamma/2,gamma))/power(r,10+nu0);
	x=(x*pow(2,4)*y)/(45*pow(P,4));
	return x;
}

//<<<<<< Computer assited proof >>>>>>

// Remark 1: The method is really slow if m increases since it has to compute two times de difference for each step (see function SSum). My recomendation is to check 
//           first the sum with m small (say 10) and see that we obtain a good interval for the L chosen (observe that increasing the L also makes the process slower).
//           Then, once we obtain a good interval, repite the process including the error with a good m for it 
// Remark 2: Observe that the imaginary part where we compute the distance is -(rho+2gamma), this means that the error ESS improves when gamma grows but
//			 on the other hand, the interval for the Simpson sum gets bigger since we are computing it at imaginary part bigger (here the difficulties with the number of digids). 
// Summarising: If we don't care about the time, take gamma as small and L, rho as big as wished to get a good Simpson sum and m as big as needed to have a good
// 				controll of the error. Observe that rho cannot be too big since then, we cannot compute the difference.

void test()//Computer assisted proof for the Stokes constant p_2^{-1}(n)
{
	int dk0=2,P=10,L=1000;//For dk0=2, we can prove it is different from zero with L=1000 (really on the limit)
	interval rho=interval(8.215),gamma=interval(0.949327);//rho=8.215 for dk_0=2
	interval nu0=nu_0(dk0),fk0=interval(1.),err,varrho=rho+2*gamma;
	IVector b(4),SS(2),b0(2);
	b[0]=b1(dk0);
	b[1]=b2(dk0);
	b[2]=b3(dk0);
	b[3]=b4(dk0);
	cout << "-- Polynomial I={d_{k_0}} -- " << endl;
	cout << "dk0 = " << dk0 << endl;
	cout << "P = " << P << endl;
	cout << "L = " << L << endl;
	cout << "rho = " << rho << endl;
	cout << "varrho = " << rho+2*gamma << endl;
	cout << "Lipschidz <=interval(0.5) " << Lip(dk0,b,rho,gamma) << " <=interval(0.5) " << endl;//It has to be less than interval(0.5) to have existence of psi_interval(1.) (see Section interval(3.)3)
	cout << "gamma = " << gamma << endl;
	cout << "c(gamma) =" << c(gamma) << " <=0.688493 " << endl;//Condition on gamma (see Proposition interval(4.)3). gamma=0.945 close to the optimal 
	cout << "Lipschidz^d =" << Lipd(dk0,b,rho,gamma) << " < 1 " << endl;//It has to be less than 1 to have existence of zeta_1 (see proof of Lemma interval(4.)6). 
	SS=exp(interval(2.)*interval::pi()*varrho)*SSum(dk0,b,L,rho,gamma,P);//Simpson sum without error. Here is where we increase the rho to obtain an smaller interval. 
	cout << "Simpson sum  = (" << SS <<")" << endl;
	err=exp(interval(2.)*interval::pi()*varrho)*ESS(dk0,b,rho,gamma,P);//Error Simpson sum. Here we see which m do we need.
	cout << "ESS " << err << endl;
	//SS[0]=SS[0]+err*interval(-1.,1.);
	//SS[1]=SS[1]+err*interval(-1.,1.);
	cout << "p_2^-1 = " << SS << endl;//Stokes constant
	b0[0]=0;
	b0[1]=(interval(2.)*interval::pi()*(dk0-interval(2.)))/(dk0-interval(1.));
	b0=power((nu0*(nu0+1))/abs(fk0),nu0)*exp(b0);
	cout << "Splitting constant = " << cprod(b0,SS) << endl; // Here we use the formula given in Theorem interval(2.)interval(4.)
}

void test3()//Loop for values of Stokes constant
{
	int dk0=21,m=200,L=2500;
	interval gamma=0.945,rho=5.25;
	interval nu0=interval(2.)/(dk0-interval(1.)),fk0=interval(interval(1.)),err;
	IVector b(4),SS(2),b0(2);
	b[0]=b1(dk0);
	b[1]=b2(dk0);
	b[2]=b3(dk0);
	b[3]=b4(dk0);
	for (int i=1;i<=29;i++)
	{
	cout << "dk0 = " << dk0 << endl;
	cout << "m = " << m << endl;
	cout << "L = " << L << endl;
	cout << "rho = " << rho << endl;
	cout << "varrho = " << rho+interval(2.)*gamma << endl;
	cout << "Lipschidz <=interval(0.5) " << Lip(dk0,b,rho,gamma) << " <=interval(0.5) " << endl;//It has to be less than interval(0.5) to have existence of psi_interval(1.) (see Section interval(3.)3)
	cout << "gamma = " << gamma << endl;
	cout << "c(gamma) <= 0.688493 " << c(gamma) << " <=0.688493 " << endl;//Condition on gamma (see Proposition interval(4.)3). gamma=0.945 close to the optimal 
	cout << "Lipschidz^d < 1 " << Lipd(dk0,b,rho,gamma) << " < 1 " << endl;//It has to be less than 1 to have existence of zeta_1 (see proof of Lemma interval(4.)6). 
	SS=SSum(dk0,b,L,rho,gamma,m);//Simpson sum without error. Here is where we increase the rho to obtain an smaller interval. 
	cout << "Simpson sum  = (" << SS <<")" << endl;
	err=exp(interval(2.)*interval::pi()*rho)*ESS(dk0,b,rho,gamma,m);//Error Simpson sum. Here we see which m do we need.
	cout << "ESS " << err << endl;
	SS[0]=SS[0]+err*interval(-1.,1.);
	SS[1]=SS[1]+err*interval(-1.,1.);
	cout << "p_2^-1 = " << SS << endl;//Stokes constant
	b0[0]=0;
	b0[1]=(interval(2.)*interval::pi()*(dk0-interval(2.)))/(dk0-interval(1.));
	b0=power((nu0*(nu0+1))/abs(fk0),nu0)*exp(b0);
	cout << "Splitting constant = " << cprod(b0,SS) << endl; // Here we use the formula given in Theorem interval(2.)interval(4.)
	dk0=dk0+1;
	}
}

void test1()//Computing the convergence of the exponent 
{
	int dk0=200,m=10,L=2500;
	interval gamma=0.945,rho=4.;
	interval nu0=interval(2.)/(dk0-interval(1.)),fk0=interval(interval(1.)),x,x1,y;
	IVector b(4),SS(2),b0(2);
	b[0]=b1(dk0);
	b[1]=b2(dk0);
	b[2]=b3(dk0);
	b[3]=b4(dk0);
	SS=SSum(dk0,b,L,rho,gamma,m);//Simpson sum without error. Here is where we increase the rho to obtain an smaller interval. 
	b0[0]=0;
	b0[1]=(interval(2.)*interval::pi()*(dk0-interval(2.)))/(dk0-interval(1.));
	b0=exp(b0);
	SS=cprod(b0,SS);
	x=log(SS[1]);
	for (int i=0;i<=100;i++)
	{
		dk0=dk0+1;
		SS=SSum(dk0,b,L,rho,gamma,m);//Simpson sum without error. Here is where we increase the rho to obtain an smaller interval. 
		b0[0]=0;
		b0[1]=(interval(2.)*interval::pi()*(dk0-interval(2.)))/(dk0-interval(1.));
		b0=exp(b0);
		SS=cprod(b0,SS);
		x1=log(SS[1]);
		//cout << "[" << log(dk0) << "," << log(lb) << "]," << endl;
		//cout << "[" << log(dk0) << "," << log(rb) << "]," << endl;
		y=(x1-x)/(log(dk0-1)-log(dk0-2));
		cout << "[" << dk0 << "," << mid(y) << "]," << endl;
		x=x1;
	}
}

void Lipschidz()//This test can be used to find the minimal rho_0 (Proposition 2.1)
{
	int dk0=950;
	interval rho=interval(4.),gamma=interval(0.949327);
	IVector b(4);
	b[0]=b1(dk0);
	b[1]=b2(dk0);
	b[2]=b3(dk0);
	b[3]=b4(dk0);
	for (int i=0;i<=10;i++)
	{
		rho=rho-interval(1.)/interval(100);
		cout << "dk0 = " << dk0 << " rho = " << rho << endl;
		cout << "Lipschidz = " << Lip(dk0,b,rho,gamma) << " <=0.5 " << endl;//It has to be less than 0.5 to have existence of psi_1.
	}
}

void testing()
{
	int dk0=2,L=3000,P=5000;
	interval varrho,rho=interval(8.215),gamma=interval(0.949327);//Existence for dk0=2 at rho=8.215
	IVector b(4),z(2);
	varrho=rho+2*gamma;
	b[0]=b1(dk0);
	b[1]=b2(dk0);
	b[2]=b3(dk0);
	b[3]=b4(dk0);
	z[0]=interval(0.);
	z[1]=-varrho;
	cout << "SSum =" << SSum(dk0,b,L,rho,gamma,1) << endl;
}

int main()
{
  	cout.precision(10);
	try
	{
		std::chrono::time_point<std::chrono::system_clock> start, end;
		start = std::chrono::system_clock::now();
		test();
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
