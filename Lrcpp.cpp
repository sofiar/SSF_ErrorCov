#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <numeric>
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
// [[Rcpp::export]]
double LL(double Alpha,double Beta,int Ng, const NumericVector &XXcov,NumericVector Ral,NumericMatrix Quienesg,
       int Npg, bool LOG)
{
  NumericVector Xcov(XXcov);
  NumericVector ral(Ral);
  NumericMatrix quienesg(Quienesg);
  double alph(Alpha);
  double b(Beta);
  int ng(Ng);
  int npg(Npg);

  
  double res=0; double pp;
  
   for(int g=0; g<=(ng-1); g++)
    {
    NumericVector Xg(npg);
    NumericVector ralg(npg);

    Xg=Xcov[quienesg(_,g)];
    ralg=ral[quienesg(_,g)];
    pp=(alph*Xg[0]+b*ralg[0])-log(sum(exp(alph*Xg+b*ralg)));

    
    res=res+pp;
    }  
   
   if (!LOG)
   {res=exp(res);}
   
  return res;
}

// [[Rcpp::export]]
NumericVector sampleX(double Alpha,double Beta,int Ng,
                      const std::vector<double>& XXcov,
                      const NumericVector &Ral,const NumericVector &Xobs,const NumericMatrix &Quienesg,
                      int Npg,int NNtot,const NumericMatrix &QQ)
{
  //NumericVector Xcov = XXcov;
  vector<double> Xcov=XXcov;
  NumericVector ral= Ral;
  NumericMatrix quienesg = Quienesg;
  NumericMatrix Q=QQ;
  double alph(Alpha);
  double b(Beta);
  int ng(Ng);
  int npg(Npg);
  int Ntot(NNtot);
  //vector<double> Xcan(Ntot);
  //vector<double> Xcan(Ntot);
  
  NumericVector Xcan(Ntot);
  NumericVector Xnew(Ntot);
  NumericVector XNcov(Ntot);

  XNcov=Xcov;
  double a0; double a1; double p1;double xx;
  
  for(int j=0; j<=(Ntot-1); j++)
  {
    Xnew=XNcov;
    //if X_j=0
    Xcan=XNcov;
    Xcan[j]=0;
    a0=LL(alph,b, Ng-1 ,Xcan,ral,Quienesg ,Npg+1,FALSE)*Q(0,Xobs[j]);
    

    //if X_j=1 
    Xcan=XNcov;
    Xcan[j]=1;
    a1=LL(alph,b, Ng-1 ,Xcan,ral,Quienesg,Npg+1,FALSE)*Q(1,Xobs[j]);
    
    p1=a1/(a0+a1);
    xx=R::rbinom(1,p1);
    
    Xnew[j]=xx;
    XNcov=Xnew;  
    
      
    
   
  }
  return XNcov;
}




