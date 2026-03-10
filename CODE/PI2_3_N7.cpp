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

//----------------------------------- Polynomial case with two monomials (Section 5.5.3) ----------------------------------------------------------

//This code gives an interval for the Stokes constant p_2^{-1}(a_2) for degree d_{k_0}=2, d_{k_1}=3 and every value of a_2 (coefficient of the second monomial to be introduced in the Stokes()).
//It uses a first approximation for phi_0 with 7 terms, that is, psi_0=psi_0^0+psi_0^1+psi_0^2+psi_0^3+psi_0^4+psi_0^5+psi_0^6+psi_0^7

//<<<<<< Implementation of some general functions used in the program >>>>>>

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

//<<<<<< Obtaining rho_0 for the existence of solutions of the inner equation (Section 5.2) >>>>>>

//Norms of the inverse operators

int pNn(int N,int n)//Definition of p_N^n
{
	interval x;
	x=interval(N-n);
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
		x=x*Pi/interval(2.);
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
//Here for the polynomial case kappa=2./(dk0-1)=2 and nu_N=18
{
	interval x;
	x=S(rho,gamma,interval(3)+4);
	x=x*(S(rho,gamma,interval(15))+S(rho,gamma,interval(22)));
	return x;
}

//Here we define the coefficients b_{j*d_{k_0}}(d_{k_0}) for j={1,2,3,4,5,6,7} of the functions psi_0^j

interval b1(interval a2)//Coefficient of psi_0^1
{
	return interval((a2-interval(10.))/interval(8.));
}

interval b2(interval a2)//Coefficient of psi_0^2
{
	interval x;
	x=interval(6.)*power(b1(a2),2)+interval(3.)*a2*b1(a2);
	x=x-interval(2.)*(interval(7.)+binom(7,4)*b1(a2));
	return x/interval(30.);
}

interval b3(interval a2)//Coefficient of psi_0^3
{
	interval x;
	x=interval(12.)*b1(a2)*b2(a2)+a2*(interval(3.)*b2(a2)+interval(3.)*power(b1(a2),2));
	x=x-interval(2.)*(interval(9.)+binom(9,6)*b1(a2)+binom(9,4)*b2(a2));
	x=x/interval(60.);
	return x;
}

interval b4(interval a2)//Coefficient of psi_0^4
{
	interval x;
	x=interval(12.)*b1(a2)*b3(a2)+interval(6.)*power(b2(a2),2)+a2*(interval(3.)*b3(a2)+interval(6.)*b1(a2)*b2(a2)+power(b1(a2),3));
	x=x-interval(2.)*(interval(11.)+binom(11,8)*b1(a2)+binom(11,6)*b2(a2)+binom(11,4)*b3(a2));
	x=x/interval(98.);
	return x;
}

interval b5(interval a2)//Coefficient of psi_0^5
{
	interval x;
	x=interval(12.)*b1(a2)*b4(a2)+interval(12.)*b2(a2)*b3(a2)+a2*(interval(3.)*b4(a2)+interval(3.)*(2*b1(a2)*b3(a2)+power(b2(a2),2))+3*power(b1(a2),2)*b2(a2));
	x=x-interval(2.)*(interval(13.)+binom(13,10)*b1(a2)+binom(13,8)*b2(a2)+binom(13,6)*b3(a2)+binom(13,4)*b4(a2));
	x=x/power(interval(12.),2);
	return x;
}

interval b6(interval a2)//Coefficient of psi_0^6
{
	interval x;
	x=interval(6.)*(interval(2.)*b1(a2)*b5(a2)+interval(2.)*b2(a2)*b4(a2)+power(b3(a2),2));
	x=x+a2*(interval(3.)*b5(a2)+interval(3.)*(2*b1(a2)*b4(a2)+2*b2(a2)*b3(a2))+3*power(b1(a2),2)*b3(a2)+3*b1(a2)*power(b2(a2),2));
	x=x-interval(2.)*(interval(15.)+binom(15,12)*b1(a2)+binom(15,10)*b2(a2)+binom(15,8)*b3(a2)+binom(15,6)*b4(a2)+binom(15,4)*b5(a2));
	x=x/((interval(14*15))-interval(12.));
	return x;
}

interval b7(interval a2)//Coefficient of psi_0^7
{
	interval x;
	x=interval(6.)*(interval(2.)*b1(a2)*b6(a2)+interval(2.)*b2(a2)*b5(a2)+2*b3(a2)*b4(a2));
	x=x+a2*(interval(3.)*b6(a2)+interval(3.)*(2*b1(a2)*b5(a2)+2*b2(a2)*b4(a2)+power(b3(a2),2))+3*power(b1(a2),2)*b4(a2)+3*b1(a2)*2*b2(a2)*b3(a2)+power(b2(a2),3));
	x=x-interval(2.)*(interval(17.)+binom(17,14)*b1(a2)+binom(17,12)*b2(a2)+binom(17,10)*b3(a2)+binom(17,8)*b4(a2)+binom(17,6)*b5(a2)+binom(17,4)*b6(a2));
	x=x/((interval(16*17))-interval(12.));
	return x;
}

//The following fuctions are used to compute the norm of the first iteration. 

interval E0(IVector b,interval rho)//The sum of E_N^n for all 0<=n<=N (Lemma 5.2.2)
{
	interval psi00,psi01,psi02,psi03,psi04,psi05,psi06,psi07;
	psi00=interval(1.)/power(rho-1.,2)+1./power(rho+interval(1.),2)-interval(2.)/power(rho,2);
	psi00=psi00*power(rho,2+interval(16.));
	for(int i=1;i<9;i++)
	{
		psi00=psi00-interval(2.)*binom(2+interval(2.)*i-interval(1.),2*i)*power(rho,16-(2*i));
	}
	psi00=abs(psi00);
	psi01=interval(1.)/power(rho-interval(1.),2+interval(2.))+interval(1.)/power(rho+interval(1.),2+interval(2.))-interval(2.)/power(rho,2+interval(2.));
	psi01=psi01*power(rho,2+interval(16.));
	for(int i=1;i<8;i++)
	{
		psi01=psi01-interval(2.)*binom(2+interval(2.)*(i+1)-interval(1.),2*i)*power(rho,16-2*(i+1));
	}
	psi01=abs(b[0])*abs(psi01);
	psi02=interval(1.)/power(rho-1,2+interval(4.))+interval(1.)/power(rho+1,2+interval(4.))-interval(2.)/power(rho,2+interval(4.));
	psi02=psi02*power(rho,2+interval(16.));
	for(int i=1;i<7;i++)
	{
		psi02=psi02-interval(2.)*binom(2+interval(2.)*(i+2)-interval(1.),2*i)*power(rho,16-2*(i+2));
	}
	psi02=abs(b[1])*abs(psi02);
	psi03=interval(1.)/power(rho-1,2+interval(6.))+interval(1.)/power(rho+1,2+interval(6.))-interval(2.)/power(rho,2+interval(6.));
	psi03=psi03*power(rho,2+interval(16.));
	for(int i=1;i<6;i++)
	{
		psi03=psi03-interval(2.)*binom(2+interval(2.)*(i+3)-interval(1.),2*i)*power(rho,16-2*(i+3));
	}
	psi03=abs(b[2])*abs(psi03);
	psi04=interval(1.)/power(rho-1,2+interval(8.))+interval(1.)/power(rho+1,2+interval(8.))-interval(2.)/power(rho,2+interval(8.));
	psi04=psi04*power(rho,2+interval(16.));
	for(int i=1;i<5;i++)
	{
		psi04=psi04-interval(2.)*binom(2+interval(2.)*(i+4)-interval(1.),2*i)*power(rho,16-2*(i+4));
	}
	psi04=abs(b[3])*abs(psi04);
	psi05=interval(1.)/power(rho-1,2+interval(10.))+interval(1.)/power(rho+1,2+interval(10.))-interval(2.)/power(rho,2+interval(10.));
	psi05=psi05*power(rho,2+interval(16.));
	for(int i=1;i<4;i++)
	{
		psi05=psi05-interval(2.)*binom(2+interval(2.)*(i+5)-interval(1.),2*i)*power(rho,16-2*(i+5));
	}
	psi05=abs(b[4])*abs(psi05);
	psi06=interval(1.)/power(rho-1,2+interval(12.))+interval(1.)/power(rho+1,2+interval(12.))-interval(2.)/power(rho,2+interval(12.));
	psi06=psi06*power(rho,2+interval(16.));
	for(int i=1;i<3;i++)
	{
		psi06=psi06-interval(2.)*binom(2+interval(2.)*(i+6)-interval(1.),2*i)*power(rho,16-2*(i+6));
	}
	psi06=abs(b[5])*abs(psi06);
	psi07=interval(1.)/power(rho-1,2+interval(14.))+interval(1.)/power(rho+1,2+interval(14.))-interval(2.)/power(rho,2+interval(14.));
	psi07=psi07*power(rho,2+interval(16.));
	for(int i=1;i<2;i++)
	{
		psi07=psi07-interval(2.)*binom(2+interval(2.)*(i+7)-interval(1.),2*i)*power(rho,16-2*(i+7));
	}
	psi07=abs(b[6])*abs(psi07);
	return psi00+psi01+psi02+psi03+psi04+psi05+psi06+psi07;
}

interval CN(IVector b,int a,interval rho)//Definition of C_N(a,rho) (see Lemma 5.2.4)
{
	interval c=interval(0.);
	for(int n=1;n<=7;n++)
	{
		c=c+abs(b[n-1])*power(rho,2*(a-n));//Absolute values of psi_0^n evaluated at rho
	}
	return c;
}

interval f0(interval a2,IVector b,interval rho)//The sum of F_N^k (Lemma 5.2.5)
{
	interval x,y,psi00,psi01,psi02,psi03,psi04,psi05,psi06,psi07;
	psi00=power(rho,-2);//Absolute values of psi0i evaluated at rho
	psi01=abs(b[0])*power(rho,-2-interval(2.));
	psi02=abs(b[1])*power(rho,-2-interval(4.));
	psi03=abs(b[2])*power(rho,-2-interval(6.));
	psi04=abs(b[3])*power(rho,-2-interval(8.));
	psi05=abs(b[4])*power(rho,-2-interval(10.));
	psi06=abs(b[5])*power(rho,-2-interval(12.));
	psi07=abs(b[6])*power(rho,-2-interval(14.));
	x=interval(2.)*(psi01*psi07)+interval(2.)*(psi02*(psi06+psi07))+2*psi03*(psi05+psi06+psi07)+power(psi05+psi06+psi07,2);//Function f(psi1) evaluated at 0
	x=interval(6.)*x;
	y=interval(3.)*power(psi00,2)*psi07+interval(3.)*psi00*(interval(2.)*psi01*(psi06+psi07)+interval(2.)*psi02*(psi05+psi06+psi07)+2*psi03*(psi04+psi05+psi06+psi07)+power(psi04+psi05+psi06+psi07,2));
	y=y+interval(3.)*power(psi01,2)*(psi05+psi06+psi07)+interval(3.)*psi01*(2*psi02*(psi04+psi05+psi06+psi07)+power(psi03+psi04+psi05+psi06+psi07,2));
	y=y+3*power(psi02,2)*(psi03+psi04+psi05+psi06+psi07)+3*psi02*power(psi03+psi04+psi05+psi06+psi07,2)+power(psi03+psi04+psi05+psi06+psi07,3);
	x=x+a2*y;
	return x*power(rho,interval(18.));
}

interval F0(interval a2,IVector b,interval rho,interval gamma)//Computation of G*E_0^b (Corollaries 5.2.7 and 5.2.9)
{
	return G(rho,gamma)*(E0(b,rho)+f0(a2,b,rho));
}

//The following are used to compute the Lipschidz constant (Lemma 5.2.11)

interval A(interval a2,IVector b,interval rho)//Boundig the norm of the linear term A (Lemma 5.2.10)
{
	interval c=interval(0.),psi;//Absolute values of psi0i evaluated at rho
	psi=CN(b,0,rho);
	c=interval(12.)*psi;
	c=c+((interval(3.)*a2)/power(rho,2))*power(1.+psi,2);
	c=c+power(rho,2)*(power(interval(1.)-power(rho,-1),2+interval(2.))+power(interval(1.)+power(rho,-1),2+interval(2.))-interval(2.))-interval(2.)*binom(2+interval(2.),2);
	return abs(c);//multiplied by rho^2 from the inverse G
}

interval R(interval a2,IVector b,interval rho,interval M0)//Bounding higher order terms R^b (Lemma 5.2.11)
{
	interval c,psi;//Absolute values of psi0i evaluated at rho
	psi=CN(b,0,rho);
	c=interval(12.)*M0/power(rho,interval(16.));
	c=c+interval(3.)*a2*(power(interval(1.)+psi+M0/power(rho,interval(16.)),2)-power(interval(1.)+psi,2));
	return abs(c);//multiplied by rho^2 from the inverse G
}

interval M0(interval a2,IVector b,interval rho,interval gamma)//Function M_0(rho,gamma)
{
	return 2*F0(a2,b,rho,gamma);
}

interval Lip(interval a2,IVector b,interval rho,interval gamma)//The Lipschidz constant for the existence (Lemma 5.2.11)
{
	return G(rho,gamma)*(A(a2,b,rho)+R(a2,b,rho,M0(a2,b,rho,gamma)));
}

//<<<<<< Obtaining rho_1 for the existence of zeta_1 and explicit formula for psi^u-psi^s (Section 5.3) >>>>>>

//Norms of the inverse operators

interval c(interval gamma)//Definition of c(gamma) involved in the domains (see definition of the domain  D_{rho,gamma}^delta)
{
	interval x;
	x=power(gamma,2)+1;
	x=power(x,interval(1.)/interval(2.));
	x=gamma/x;
	return x;
}

interval Sd(interval rho,interval gamma,interval nu)// Bound of the first order inverse operator S^d  (Proposition 5.3.2)
{
	interval s,r,x1,x2,Pi=interval::pi();
	r=rho+(interval(3.)*gamma)/interval(2.);
	if (leftBound(nu)>0){
		s=c(gamma)-interval(3.);
		x1=power(interval(1.)-c(gamma)/interval(2.),2)/power(r,2);
		x1=power(1+x1,nu/interval(2.));
		x1=(interval(1.)/(interval(2.)*Pi*power(c(gamma),2)*r))*x1;
		x2=(interval(3.)-2*c(gamma))/(interval(4.)*power(r,2));
		x2=power(1.-x2,-nu/interval(2.));
		x2=(interval(1.)+interval(1.)/c(gamma))*B(nu)*x2;
		x1=x1+x2;
	}
	else{
		s=interval(5.)-c(gamma);
		if (leftBound(nu)<-1){
			x1=interval(1.)+(nu+1)/(2*Pi*r);
			x1=1+interval(1.)/x1;
			x1=1+x1/(interval(2.)*Pi*c(gamma)*r);
			x2=power(1-c(gamma)/interval(2.),2)+power(r,2);
			x2=1+(interval(5.)-interval(2.)*c(gamma))/(interval(4.)*x2);
			x2=power(x2,-nu/interval(2.));
			x1=x1*x2;
		}
		else{
			x1=power(1-c(gamma)/interval(2.),2)/power(r,2);
			x1=power(1+x1,-nu/interval(2.));
			x1=x1/(-nu*c(gamma));
			x2=(1+interval(1.)/c(gamma))/(interval(2.)*Pi*power(c(gamma),2)*r);
			x1=x1+x2;
		}
	}
	x2=power(interval(2.)-c(gamma)/interval(2.),2)+power(r,2);
	x2=1+s/x2;
	x2=power(x2,-nu/interval(2.));
	x1=interval(1.)/r+x1*x2;
	return x1;
}

//In order to define the bound for the last inverse, we need to define the the following two bounds related to eta_1^d and eta_2^d

interval eta1(interval a2,IVector b,interval rho,interval gamma)//Computation of eta_1^b (Lemma 5.3.3)
{
	interval r,x;
	r=(rho+interval(3.)*gamma/interval(2.));
	x=interval(2.);
	for(int n=1;n<=7;n++)
	{
		x=x+2*(n+1)*abs(b[n-1])*power(rho,-2*n);//Absolute values of dpsi_0^n evaluated at rho
	}
	x=x+(2*M0(a2,b,rho,gamma))/(c(gamma)*power(rho,15));
	return x+(3*M0(a2,b,rho+gamma/2,gamma))/(power(r,16));
}

interval eta2(interval a2,IVector b,interval rho,interval gamma)//Computation of eta_2^b (Corollary 5.3.6)
{
	interval r,x,eta,nu0=interval(2.);
	r=rho+3*gamma/2;
	eta=eta1(a2,b,rho,gamma);
	x=1+power(r,-1);
	x=Sd(rho,gamma,-interval(2*nu0+3))*eta*power((2*nu0)-eta,-2)*power(x,nu0+1);
	return x;
}

interval Gd(interval a2,IVector b,interval rho,interval gamma)// Bound of the second order inverse operator G^d (Corollary 5.3.7)
//Here for the polynomial case kappa=2 and nuN=18
{
	return eta1(a2,b,rho,gamma)*eta2(a2,b,rho,gamma)*(Sd(rho,gamma,interval(16))+Sd(rho,gamma,interval(23)));
}

interval Ad(interval a2,IVector b,interval rho, interval gamma)//Computation of A^{d,b} (Lemma 5.3.1)
{
	interval c,r,M;
	r=rho+interval(3.)*gamma/interval(2.);
	M=M0(a2,b,rho+gamma/interval(2.),gamma);
	c=interval(2.)*M*(interval(12.)+((interval(6.)*a2)/power(r,interval(2.)))*(1+CN(b,0,r)+(M/power(r,16))));
	return c;
}

interval Lipd(interval a2,IVector b,interval rho,interval gamma)//Bound for the operator G^d[A^d*] (Lemma 5.3.8)
{
	return Gd(a2,b,rho,gamma)*Ad(a2,b,rho,gamma)*power(rho+interval(3.)*gamma/interval(2.),-16);
}

interval zeta11(interval a2,IVector b,interval rho,interval gamma)//Computation of zeta_1^{1,b} (Lemma 5.3.8)
{
	return (eta1(a2,b,rho,gamma)*Gd(a2,b,rho,gamma)*Ad(a2,b,rho,gamma))/(1-Lipd(a2,b,rho,gamma));
}

//<<<<<< Rigorous computations of psi^u-psi^s and dpsi^s (Section 5.4) >>>>>>

//Initial aproximation of the manifold and its derivative

IVector psi0(IVector b,IVector z)//Computation of psi_0=psi_0^0+psi_0^1+psi_0^2+psi_0^3+psi_0^4+psi_0^5+psi_0^6+psi_0^7
{
	IVector x(2);
	x=power(z,-interval(2.));
	for(int n=1;n<=7;n++)
	{
		x=x+b[n-1]*power(z,-2*n-interval(2.));//psi_0^n evaluated at z
	}
	return x;
}

IVector psi0(IVector b,interval re,interval im)//Computation of psi0(re,im)
{
	IVector z(2);
	z[0]=re;
	z[1]=im;
	return psi0(b,z);
}

IVector dpsi0(IVector b,IVector z)//Derivative of psi_0=psi_0^0+psi_0^1+psi_0^2+psi_0^3+psi_0^4+psi_0^5+psi_0^6+psi_0^7
{
	IVector x(2);
	x=-2*power(z,-1-2);
	for(int n=1;n<=7;n++)
	{
		x=x-(2*(n+1))*b[n-1]*power(z,-(2*n+1)-2);//dpsi_0^n evaluated at z
	}
	return x;
}

IVector dpsi0(IVector b,interval re,interval im)// Derivative of psi_0(re,im)
{
	IVector z(2);
	z[0]=interval(re);
	z[1]=im;
	return dpsi0(b,z);
}

//<<<<<< Iterative methods for psi^u, psi^s (Section 5.5.1) and derivative psi^s (Section 5.5.2) >>>>>>

//The following functions are defined to obtain the generating functions F and dF for the iterative methods

IVector T0(IVector b,IVector z)//T_0^Omega(z) for psi00+psi01+psi02+psi03+psi04+psi05+psi06+psi07 (Lemma 5.4.4)
{
	int Omega=20;
	interval nu0=interval(2.);
	IVector psi00(2),psi01(2),psi02(2),psi03(2),psi04(2),psi05(2),psi06(2),psi07(2);
	psi00[0]=interval(0.);
	psi00[1]=psi00[0];
	psi01=psi00;
	psi02=psi00;
	psi03=psi00;
	psi04=psi00;
	psi05=psi00;
	psi06=psi00;
	psi07=psi00;
	for(int i=9;i<=pNn(Omega,0);i++)
	{
		psi00=psi00+binom(nu0+interval(2.)*i-1,2*i)*power(z,-nu0-interval(2.)*i);
	}
	for(int i=8;i<=pNn(Omega,1);i++)
	{
		psi01=psi01+binom(nu0+interval(2.)*(i+1)-interval(1.),2*i)*power(z,-nu0-interval(2.)*(i+1));
	}
	psi01=b[0]*psi01;
	for(int i=7;i<=pNn(Omega,2);i++)
	{
		psi02=psi02+binom(nu0+interval(2.)*(i+2)-interval(1.),2*i)*power(z,-nu0-interval(2.)*(i+2));
	}
	psi02=b[1]*psi02;
	for(int i=6;i<=pNn(Omega,3);i++)
	{
		psi03=psi03+binom(nu0+interval(2.)*(i+3)-interval(1.),2*i)*power(z,-nu0-interval(2.)*(i+3));
	}
	psi03=b[2]*psi03;
	for(int i=5;i<=pNn(Omega,4);i++)
	{
		psi04=psi04+binom(nu0+interval(2.)*(i+4)-interval(1.),2*i)*power(z,-nu0-interval(2.)*(i+4));
	}
	psi04=b[3]*psi04;
	for(int i=4;i<=pNn(Omega,5);i++)
	{
		psi05=psi05+binom(nu0+interval(2.)*(i+5)-interval(1.),2*i)*power(z,-nu0-interval(2.)*(i+5));
	}
	psi05=b[4]*psi05;
	for(int i=3;i<=pNn(Omega,6);i++)
	{
		psi06=psi06+binom(nu0+interval(2.)*(i+6)-interval(1.),2*i)*power(z,-nu0-interval(2.)*(i+6));
	}
	psi06=b[5]*psi06;
	for(int i=2;i<=pNn(Omega,7);i++)
	{
		psi07=psi07+binom(nu0+interval(2.)*(i+7)-interval(1.),2*i)*power(z,-nu0-interval(2.)*(i+7));
	}
	psi07=b[6]*psi07;
	return -2*(psi00+psi01+psi02+psi03+psi04+psi05+psi06+psi07);
}

IVector T0b(IVector b,IVector z)///Sum of T_0^{n,b}(z) (Lemma 5.4.4 and Remark 5.4.5)
{
	int Omega=20;
	interval psi00,psi01,psi02,psi03,psi04,psi05,psi06,psi07,c,nu0=interval(2.);
	IVector x(2);
	c=mod(z);
	psi00=nu0*((2*pNn(Omega,0)-1)/fac(2*(pNn(Omega,0)+1)))*(power(c-1,-(2*pNn(Omega,0)+nu0))+power(c+1,-(2*pNn(Omega,0)+nu0))-2*power(c,-(2*pNn(Omega,0)+nu0)));
	psi01=abs(b[0])*binom(2*pNn(Omega,1)+1+nu0,2*pNn(Omega,1))*(power(c-1,-(2+2*pNn(Omega,1)+nu0))+power(c+1,-(2+2*pNn(Omega,1)+nu0))-2*power(c,-(2+2*pNn(Omega,1)+nu0)));
	psi02=abs(b[1])*binom(2*2+2*pNn(Omega,2)-1+nu0,2*pNn(Omega,2))*(power(c-1,-(2*2+2*pNn(Omega,2)+nu0))+power(c+1,-(2*2+2*pNn(Omega,2)+nu0))-2*power(c,-(2*2+2*pNn(Omega,2)+nu0)));
	psi03=abs(b[2])*binom(2*3+2*pNn(Omega,3)-1+nu0,2*pNn(Omega,3))*(power(c-1,-(2*3+2*pNn(Omega,3)+nu0))+power(c+1,-(2*3+2*pNn(Omega,3)+nu0))-2*power(c,-(2*3+2*pNn(Omega,3)+nu0)));
	psi04=abs(b[3])*binom(2*4+2*pNn(Omega,4)-1+nu0,2*pNn(Omega,4))*(power(c-1,-(2*4+2*pNn(Omega,4)+nu0))+power(c+1,-(2*4+2*pNn(Omega,4)+nu0))-2*power(c,-(2*4+2*pNn(Omega,4)+nu0)));
	psi05=abs(b[4])*binom(2*5+2*pNn(Omega,5)-1+nu0,2*pNn(Omega,5))*(power(c-1,-(2*5+2*pNn(Omega,5)+nu0))+power(c+1,-(2*5+2*pNn(Omega,5)+nu0))-2*power(c,-(2*5+2*pNn(Omega,5)+nu0)));
	psi06=abs(b[5])*binom(2*6+2*pNn(Omega,6)-1+nu0,2*pNn(Omega,6))*(power(c-1,-(2*6+2*pNn(Omega,6)+nu0))+power(c+1,-(2*6+2*pNn(Omega,6)+nu0))-2*power(c,-(2*6+2*pNn(Omega,6)+nu0)));
	psi07=abs(b[6])*binom(2*7+2*pNn(Omega,7)-1+nu0,2*pNn(Omega,7))*(power(c-1,-(2*7+2*pNn(Omega,7)+nu0))+power(c+1,-(2*7+2*pNn(Omega,7)+nu0))-2*power(c,-(2*7+2*pNn(Omega,7)+nu0)));
	x[0]=interval(-1.,1.)*(psi00+psi01+psi02+psi03+psi04+psi05+psi06+psi07);
	x[1]=x[0];
	return x;
}

IVector T1(interval a2,IVector b,IVector psi,IVector z)//T_1^Omega for psi00+psi01+psi02+psi03+psi04+psi05 at point z (Lemma 5.4.4)
{
	IVector x(2),y(2),psi00(2),psi01(2),psi02(2),psi03(2),psi04(2),psi05(2),psi06(2),psi07(2);
	psi00=power(z,-2);
	psi01=b[0]*power(z,-2-interval(2.));
	psi02=b[1]*power(z,-2-interval(4.));
	psi03=b[2]*power(z,-2-interval(6.));
	psi04=b[3]*power(z,-2-interval(8.));
	psi05=b[4]*power(z,-2-interval(10.));
	psi06=b[5]*power(z,-2-interval(12.));
	psi07=b[6]*power(z,-2-interval(14.));
	x=interval(2.)*cprod(psi00,psi)+interval(2.)*cprod(psi01,psi07+psi)+interval(2.)*cprod(psi02,psi06+psi07+psi)+2*cprod(psi03,psi05+psi06+psi07+psi)+power(psi04+psi05+psi06+psi07+psi,2);
	x=interval(6.)*x;
	y=interval(3.)*cprod(power(psi00,2),psi07+psi)+interval(3.)*cprod(psi00,interval(2.)*cprod(psi01,psi06+psi07+psi)+2*cprod(psi02,psi05+psi06+psi07+psi)+2*cprod(psi03,psi04+psi05+psi06+psi07+psi)+power(psi04+psi05+psi06+psi07+psi,2));
	y=y+interval(3.)*cprod(power(psi01,2),psi05+psi06+psi07+psi)+interval(3.)*cprod(psi01,2*cprod(psi02,psi04+psi05+psi06+psi07+psi)+power(psi03+psi04+psi05+psi06+psi07+psi,2));
	y=y+3*cprod(power(psi02,2),psi03+psi04+psi05+psi06+psi07+psi)+3*cprod(psi02,power(psi03+psi04+psi05+psi06+psi07+psi,2))+power(psi03+psi04+psi05+psi06+psi07+psi,3);
	x=x+a2*y;
	return x;
}

IVector F(interval a2,IVector b,IVector psi,IVector z)//Right hand side of the iterative method for psi_1^u and psi_1^s
{
	return T0(b,z)+T1(a2,b,psi,z)+T0b(b,z);
}

IVector dF(interval a2,IVector b,IVector psi,IVector z)//Right hand side of the iterative method for the derivative of psi_1^s
{
	IVector Psi0=psi0(b,z);
	return interval(12.)*(Psi0+psi)+interval(3.)*a2*power(Psi0+psi,2);
}

// We treat Fu and Fs separately only for the sake of documentation. The function is the same, but the justification of the formulae is written out differently.

IVector Fu(interval a2,IVector b,IVector x,IVector z)//Iterative function for psi_1^u
{
	// below the convention is
	// d=d(n)=psi(n)-psi(n-1)
	// psi=psi(n)
	IVector d(2),psi(2);
	split(x,d,psi);
	// below we comute d(n+1):
	d = F(a2,b,psi,z)+d;
	// Note that 
	// d(n+1) = psi(n+1)-psi(n)
	// psi(n+1) = psi(n)+d(n+1):
	psi = psi+d;
	return combine(d,psi);
}

IVector Fs(interval a2,IVector b,IVector x,IVector z)//Iterative function for psi_1^s
{
	// below the convention is:
	// d=d(n)=psi(n)-psi(n+1)
	// psi=psi(n)
	IVector d(2),psi(2);
	split(x,d,psi);
	// below we compute d(n-1):
	d = F(a2,b,psi,z)+d;
	// Note that
	// d(n-1) = psi(n-1) - psi(n)
	// psi(n-1) = psi(n) + d(n-1)
	psi = psi + d;
	return combine(d,psi);
}

IVector dFs(interval a2,IVector b,IVector dx,IVector z)//Iterative function for the derivative of psi_1^s
{
	// below the convention is:
	// d=d(n)=psi(n)-psi(n+1)
	// psi=psi(n)
	// dd=dd(n)=dpsi(n)-dpsi(n+1)
	// dpsi=dpsi(n)
	IVector d(2),psi(2),dd(2),dpsi(2);
	dsplit(dx,d,psi,dd,dpsi);
	// below we compute d(n-1) and dd(n-1):
	d = F(a2,b,psi,z)+d;
	dd = cprod(dF(a2,b,psi,z),dpsi)+dd;
	// Note that:
	// d(n-1) = psi(n-1) - psi(n) (same for dd(n-1)
	// psi(n-1) = psi(n) + d(n-1) (same for dpsi(n-1)
	psi = psi + d;
	dpsi = dpsi + dd;
	return dcombine(d,psi,dd,dpsi);
}

//Initial conditions for the iterative methods

IVector psi10(interval a2,IVector b,interval n,interval rho,interval gamma)//Initial condition for psi_1^{u,s} (using Proposition 5.1.2)
{
	interval varrho=rho+2*gamma;
	IVector z({n,-varrho}),psi10(2);
	psi10[0]=interval(-1.,1.)*M0(a2,b,rho,gamma)/(power(mod(z),18));
	psi10[1]=psi10[0];
	return psi10;
}

IVector d0_u(interval a2,IVector b,interval n,interval rho,interval gamma)//Initial condition for d^u
{
	return psi10(a2,b,n,rho,gamma) - psi10(a2,b,n-1,rho,gamma);
}

IVector d0_s(interval a2,IVector b,interval n,interval rho,interval gamma)//Initial condition for d^s 
{
	return psi10(a2,b,n,rho,gamma) - psi10(a2,b,n+1,rho,gamma);
}

IVector x0_u(interval a2,IVector b,interval n,interval rho,interval gamma)// Initial conditions unstable combined
{
	return combine(d0_u(a2,b,n,rho,gamma),psi10(a2,b,n,rho,gamma));
}

IVector x0_s(interval a2,IVector b,interval n,interval rho,interval gamma)// Initial conditions stable combined
{
	return combine(d0_s(a2,b,n,rho,gamma),psi10(a2,b,n,rho,gamma));
}

IVector dpsi10(interval a2,IVector b,interval n,interval rho,interval gamma)//Initial condition for the derivative of psi_1^s (see Remark 5.3.4)
{
	interval M1,nu0=interval(2.),varrho=rho+2*gamma;
	IVector z({n,-varrho}),dpsi10(2);
	M1=2*M0(a2,b,rho,gamma)/c(gamma);
	M1=M1+((16+nu0)*M0(a2,b,rho+gamma/2,gamma))/(rho+interval(3.)*gamma/interval(2.));//Here we use that we are computing the derivative at im z= -i(rho+2gamma)
	M1=M1/(power(mod(z),16+nu0));
	dpsi10[0]=interval(-1.,1.)*M1;
	dpsi10[1]=dpsi10[0];
	return dpsi0(b,z)+dpsi10;
}

IVector dd0_s(interval a2,IVector b,interval n,interval rho,interval gamma)//Initial condition for derivative d^s 
{
	return dpsi10(a2,b,n,rho,gamma) - dpsi10(a2,b,n+1,rho,gamma);
}

IVector dx0_s(interval a2,IVector b,interval n,interval rho,interval gamma)// Initial conditions derivative stable combined
{
	return dcombine(d0_s(a2,b,n,rho,gamma),psi10(a2,b,n,rho,gamma),dd0_s(a2,b,n,rho,gamma),dpsi10(a2,b,n,rho,gamma));//Using Cauchy estimates with r=1/2
}

//Definition of the iterative methods

IVector Wu(interval a2,IVector b,interval re,int L,interval rho,interval gamma)// psi_1^u at re-1-ivarrho and re-ivarrho interated from -L
{
	interval varrho=-rho-2*gamma;//Height where to compute the difference (where we have defined the interval re(z) in [0,1])
	IVector z(2),x=x0_u(a2,b,interval(-L+re),rho,gamma),psi1(2),psi2(2),d(2);
	z[0]=interval(-L+re);
	z[1]=varrho;
	for(int i=0;i<L-1;i++)
	{
		x=Fu(a2,b,x,z);
		z[0]=z[0]+1;
	}
	split(x,d,psi1);
	x=Fu(a2,b,x,z);
	split(x,d,psi2);
	x=combine(psi1,psi2);
	return x;
}

IVector dWs(interval a2,IVector b,interval re,int L,interval rho,interval gamma)// psi_1^s and its derivative at re-ivarrho and re+1-ivarrho interated from L
{
	interval varrho=-rho-2*gamma;//Height where to compute the difference (where we have defined the interval re(z) in [0,1])
	IVector z(2),x=dx0_s(a2,b,interval(L+re),rho,gamma),psi1(2),psi2(2),d(2),dpsi1(2),dpsi2(2),dd(2);
	z[0]=interval(L+re);
	z[1]=varrho;
	for(int i=0;i<L-1;i++)
	{
		x=dFs(a2,b,x,z);
		z[0]=z[0]-1;
	}
	dsplit(x,d,psi1,dd,dpsi1);
	x=dFs(a2,b,x,z);
	dsplit(x,d,psi2,dd,dpsi2);
	x=dcombine(psi1,psi2,dpsi1,dpsi2);
	return x;
}

//<<<<<< Definition of the Simpson sum and its error (Section 5.1.2) >>>>>>

IVector SSum(interval a2,IVector b,int L,interval rho,interval gamma,int P)// Simpson sum of p2(z)*e^(2pi*z) between -ivarrho and 1-ivarrho 
{
	interval s,varrho,nu0=interval(2.);
	IVector err(2),SS(2),x(4),y(2),dx(8),d(2),dd(2),psi1sL(2),psi1uL(2),psi1sR(2),psi1uR(2),dpsisL(2),dpsisR(2),z(2),deltaL(2),deltaR(2);
	varrho=rho+2*gamma;
	SS[0]=interval(0.);
	SS[1]=interval(0.);
	for(int i=0;i<P+1;i++)
	{
		s=interval(i)/interval(P);
		x=Wu(a2,b,s,L,rho,gamma);//Computing psi_1^u at -ivarrho+s-1 and -ivarrho+s
		split(x,psi1uL,psi1uR);
		dx=dWs(a2,b,s-1,L,rho,gamma);//Computing psi_1^s and its derivative at -ivarrho+s-1 and -ivarrho+s
		dsplit(dx,psi1sR,psi1sL,dpsisR,dpsisL);
		deltaL=psi1uL-psi1sL;
		deltaR=psi1uR-psi1sR;
		z[1]=-varrho;
		z[0]=interval(s-1);
		err=power(mod(z),-17-nu0)*zeta11(a2,b,rho,gamma)*interval(-1.,1.);//Observe that this error decreses when considering a bigger rho, which means smaller interval in the end
		y=cprod(deltaR,dpsisL+err);
		z[0]=interval(s);
		err=power(mod(z),-17-nu0)*zeta11(a2,b,rho,gamma)*interval(-1.,1.);
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

interval ESS(interval a2,IVector b,interval rho,interval gamma,int P)//Error for the Simpson sum (Lemma 5.1.6)
{
	interval r,a,x,y,nu0=interval(2.);
	a=c(gamma);
	r=rho+2*gamma-a/2;
	y=interval(0.);
	for (int i=0;i<5;i++)
	{
		y=y+(fac(4)/fac(i))*power(interval::pi(),i)*power(a,i-4);
	}
	x=eta1(a2,b,rho,gamma)/power(r,nu0+1);
	x=x+(zeta11(a2,b,rho,gamma)/power(r,17+nu0));
	x=(x*M0(a2,b,rho+gamma/2,gamma))/power(r,16+nu0);
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

void Stokes()//Computer assisted proof for the Stokes constant p_2^{-1}(a_2)
{
	int P=10000,L=200;
	interval rho=interval(9.155),gamma=interval(0.949327);
	interval a2=interval(6.02),fk0=interval(1.),err,varrho=rho+2*gamma;
	IVector b(7),SS(2);
	b[0]=b1(a2);
	b[1]=b2(a2);
	b[2]=b3(a2);
	b[3]=b4(a2);
	b[4]=b5(a2);
	b[5]=b6(a2);
	b[6]=b7(a2);
	cout << "-- Polynomial I={2,3} -- " << endl;
	cout << "a_2 = " << a2 << endl;
	cout << "P = " << P << endl;
	cout << "L = " << L << endl;
	cout << "rho = " << rho << endl;
	cout << "varrho = " << rho+2*gamma << endl;
	cout << "Lipschidz =" << Lip(a2,b,rho,gamma) << " <=interval(0.5) " << endl;//It has to be less than 1/2 to have existence of psi_1 (see Section 5.2.3)
	cout << "gamma = " << gamma << endl;
	cout << "c(gamma) =" << c(gamma) << " <=0.688493 " << endl;//Condition on gamma (see Proposition 5.3.2).
	cout << "Lipschidz^d =" << Lipd(a2,b,rho,gamma) << " < 1 " << endl;//It has to be less than 1 to have existence of zeta_1 (see proof of Lemma 5.3.8). 
	SS=exp(interval(2.)*interval::pi()*varrho)*SSum(a2,b,L,rho,gamma,P);//Simpson sum without error.  
	cout << "Simpson sum  = (" << SS <<")" << endl;
	err=exp(interval(2.)*interval::pi()*varrho)*ESS(a2,b,rho,gamma,P);//Error Simpson sum. Here we see which P we need.
	cout << "ESS " << err << endl;
	SS[0]=SS[0]+err*interval(-1.,1.);
	SS[1]=SS[1]+err*interval(-1.,1.);
	cout << "p_2^-1 = " << SS << endl;//Stokes constant
	cout << "Splitting constant = " << power(interval(6.)/abs(fk0),2)*SS << endl;// Here we use Definition 1.2.13
}

void Lipschidz()//This test can be used to find the minimal rho_0(d_{k_0}) (Section 5.2.3)
{
	interval a2=interval(3.);
	interval rho=interval(8.7),gamma=gamma=interval(0.949327);
	IVector b(7);
	b[0]=b1(a2);
	b[1]=b2(a2);
	b[2]=b3(a2);
	b[3]=b4(a2);
	b[4]=b5(a2);
	b[5]=b6(a2);
	b[6]=b7(a2);
	for (int i=0;i<=10;i++)
	{
		rho=rho-interval(1.)/interval(100);
		cout << "a_2 = " << a2 << " rho = " << rho << endl;
		cout << "Lipschidz = " << Lip(a2,b,rho,gamma) << " <=0.5 " << endl;//It has to be less than 0.5 to have existence of psi_1.
	}
}

int main()
{
  	cout.precision(10);
	try
	{
		std::chrono::time_point<std::chrono::system_clock> start, end;
		start = std::chrono::system_clock::now();
		//Lipschidz();
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