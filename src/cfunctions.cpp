# include <numeric>
# include <algorithm>
# include <Rcpp.h>

// [[Rcpp::depends(Rcpp)]]

using namespace Rcpp ;

// [[Rcpp::export]]
double fact(int i)
{
  return (i>1 ? i*fact(i-1) : 1.0);
}




// [[Rcpp::export]]
double innerProduct(NumericVector x, NumericVector y) {
  return std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
}


// [[Rcpp::export]]
double likelihoodGP2(double sumlimit, double lambda, double alpha1, double alpha2,
                                 double alpha3, double eta, int T, int seas1, int seas2,
                                 std::vector<int> data)
{
  double U = 1/( 1 - alpha1-alpha2-alpha3);

  double beta1 = lambda * U *alpha1;
  double beta2 = lambda * U *alpha2;
  double beta3 = lambda * U *alpha3;
  double zeta = lambda*U*(1-2*alpha1-alpha3);
  double pyz;
  double gryz;
  double nain;
  double help;
  double dwarf;
  double zwerg;

  int t;
  double mlef = 0.0;
  for (t = seas2 + 1; t <= T;t++ ) {
    int x = data[t-1];
    int y = data[t-seas1-1];
    int z = data[t-seas2-1];

    int minxyz  = std::min(x, y+z);
    int r;
    double nain = 0.0;
    for (r=0; r <= minxyz;r++) {
      int s, v, w;
      zwerg=0.0;
      for ( s=0; s <= y; s++)
      {
        for ( v=0; v <= y; v++)
        {
          for ( w=0; w <= y; w++)
          {
            if (r-s-v >= 0 && z-r+v-w >= 0 && y-s-v-w >= 0){
              zwerg = zwerg  + beta3*pow(beta3+s*eta,s-1)/fact(s)*exp(-beta3-s*eta) * beta1*pow(beta1+v*eta,v-1)/fact(v)*exp(-beta1-v*eta) * beta1*pow(beta1+w*eta,w-1)/fact(w)*exp(-beta1-w*eta) * beta2*pow(beta2+(r-s-v)*eta,r-s-v-1)/fact(r-s-v)*exp(-beta2-(r-s-v)*eta) * lambda*pow(lambda+(z-r+v-w)*eta,z-r+v-w-1)/fact(z-r+v-w)*exp(-lambda-(z-r+v-w)*eta) * zeta*pow(zeta+(y-s-v-w)*eta,y-s-v-w-1)/fact(y-s-v-w)*exp(-zeta-(y-s-v-w)*eta);
            } //end if
          } //for w
        } //for v
      } // for s

      help = (beta1+beta3) * pow(lambda*U*(1-alpha1-alpha3),2) * exp(-(beta1+beta3)-2*(lambda*U*(1-alpha1-alpha3))-y*eta-z*eta);
      int minyz  = std::min(y, z);

      int j;
      dwarf = 0.0;
      for (j=0;j <=minyz; j++) {
        dwarf = dwarf + pow(lambda*U*(1-alpha1-alpha3) + eta*(y-j),y-j-1) / fact(y-j) * pow(lambda*U*(1-alpha1-alpha3) + eta*(z-j),z-j-1) /fact(z-j) * pow(lambda*U*(alpha1+alpha3) + eta*j,j-1) /fact(j) * exp(j*eta);
      }

      pyz = dwarf*help;
      gryz = zwerg*pow(pyz,-1);

      nain = nain + gryz*lambda*pow(lambda+(x-r)*eta,x-r-1)*exp(-lambda-(x-r)*eta) /fact(x-r);
    } //end r

    mlef = mlef + log(nain);
  } //end t

  return (-mlef);
} //end of function


// [[Rcpp::export]]
double likelihoodGP1(double sumlimit, double lambda, double alpha, double eta, int T, int seas, std::vector<int> data)
{

  double psi = eta*(1-alpha)/lambda;
  double gry;

  int t;

  double mlef = 0.0;
  for (t = seas + 1; t <= T; t++) {

    int x = data[t-1];
    int y = data[t-seas-1];


    int minxy = std::min(x,y);

    int r;
    double nain = 0.0;
    for (r=0;r <=minxy; r++) {
      if (y >= r) {
        nain = nain + fact(y)/fact(y-r)/fact(r) * alpha *(1-alpha) * pow(alpha+psi*r,r-1)*pow(1-alpha+psi*(y-r),y-r-1)/pow(1+psi*y,y-1) *
          lambda*pow(lambda+(x-r)*eta,x-r-1)*exp(-lambda-(x-r)*eta)/fact(x-r);
      } //end if

    } //end r
    mlef = mlef + log(nain);

  } // end t

  return(-mlef);
} // end of function


// [[Rcpp::export]]
std::vector<int> simGP2(double sumlimit, double lambda, double alpha1, double alpha2,
                        double alpha3, double eta, int T,int N, int seas1, int seas2,
                        std::vector<int> data, std::vector<double> uniform, std::vector<int> innovations)
{
  double U = 1/( 1 - alpha1-alpha2-alpha3);

  double beta1 = lambda * U *alpha1;
  double beta2 = lambda * U *alpha2;
  double beta3 = lambda * U *alpha3;
  double zeta = lambda*U*(1-2*alpha1-alpha3);
  double pyz;
  double gryz;
  double nain;
  double help;
  double dwarf;
  double zwerg;

  int t;
  for (t = N + 1; t <= T; t++ ) {
    int y = data[t-seas1-1];
    int z = data[t-seas2-1];

    double unif = uniform[t-1];

    int r = 0.0;
    double nain = 0.0;
    while (nain <= unif) {
      int s, v, w;
      zwerg=0.0;
      for ( s=0; s <= sumlimit; s++) //y
      {
        for ( v=0; v <= sumlimit; v++)
        {
          for ( w=0; w <= sumlimit; w++)
          {
            if (r-s-v >= 0 && z-r+v-w >= 0 && y-s-v-w >= 0){
              zwerg = zwerg  + beta3*pow(beta3+s*eta,s-1)/fact(s)*exp(-beta3-s*eta) * beta1*pow(beta1+v*eta,v-1)/fact(v)*exp(-beta1-v*eta) * beta1*pow(beta1+w*eta,w-1)/fact(w)*exp(-beta1-w*eta) * beta2*pow(beta2+(r-s-v)*eta,r-s-v-1)/fact(r-s-v)*exp(-beta2-(r-s-v)*eta) * lambda*pow(lambda+(z-r+v-w)*eta,z-r+v-w-1)/fact(z-r+v-w)*exp(-lambda-(z-r+v-w)*eta) * zeta*pow(zeta+(y-s-v-w)*eta,y-s-v-w-1)/fact(y-s-v-w)*exp(-zeta-(y-s-v-w)*eta);
            } //end if
          } //for w
        } //for v
      } // for s

      help = (beta1+beta3) * pow(lambda*U*(1-alpha1-alpha3),2) * exp(-(beta1+beta3)-2*(lambda*U*(1-alpha1-alpha3))-y*eta-z*eta);
      int minyz  = std::min(y, z);

      int j;
      dwarf = 0.0;
      for (j=0;j <=minyz; j++) {
        dwarf = dwarf + pow(lambda*U*(1-alpha1-alpha3) + eta*(y-j),y-j-1) / fact(y-j) * pow(lambda*U*(1-alpha1-alpha3) + eta*(z-j),z-j-1) /fact(z-j) * pow(lambda*U*(alpha1+alpha3) + eta*j,j-1) /fact(j) * exp(j*eta);
      }

      pyz = dwarf*help;
      gryz = zwerg*pow(pyz,-1);

      nain = nain + gryz;
      r = r+1;
    } //end r

    data[t-1] = r-1 + innovations[t-1];
  } //end t

  return (data);
} //end of function


// [[Rcpp::export]]
std::vector<int> simGP1(double sumlimit, double lambda, double alpha, double eta, int T, int N,
                        int seas, std::vector<int> data,  std::vector<double> uniform, std::vector<int> innovations)
{

  double psi = eta*(1-alpha)/lambda;
  double gry;

  int t;

  double mlef = 0.0;
  for (t = N + 1; t <= T; t++) {


    int y = data[t-seas-1];

    double unif = uniform[t-1];

    int r = 0.0;
    double nain = 0.0;
    while (nain <= unif) {
      if (y >= r) {
        nain = nain + fact(y)/fact(y-r)/fact(r) * alpha *(1-alpha) * pow(alpha+psi*r,r-1) *pow(1-alpha+psi*(y-r),y-r-1)/pow(1+psi*y,y-1);
      } //end if
      r = r+1;
    } //end r
    data[t-1] = r-1 + innovations[t-1];

  } // end t

  return(data);
} // end of function


// [[Rcpp::export]]
double dGP2h(int x, int y, int z, double lambda, double alpha1, double alpha2, double alpha3, double eta)
{
  double U = 1/( 1 - alpha1-alpha2-alpha3);

  double beta1 = lambda * U *alpha1;
  double beta2 = lambda * U *alpha2;
  double beta3 = lambda * U *alpha3;
  double zeta = lambda*U*(1-2*alpha1-alpha3);
  double pyz;
  double gryz;
  double help;
  double dwarf;
  double zwerg;


  int minxyz  = std::min(x, y+z);
  int r;
  double nain = 0.0;
  for (r=0; r <= minxyz;r++) {
    int s, v, w;
    zwerg=0.0;
    for ( s=0; s <= y; s++)
    {
      for ( v=0; v <= y; v++)
      {
        for ( w=0; w <= y; w++)
        {
          if (r-s-v >= 0 && z-r+v-w >= 0 && y-s-v-w >= 0){
            zwerg = zwerg  + beta3*pow(beta3+s*eta,s-1)/fact(s)*exp(-beta3-s*eta) * beta1*pow(beta1+v*eta,v-1)/fact(v)*exp(-beta1-v*eta) * beta1*pow(beta1+w*eta,w-1)/fact(w)*exp(-beta1-w*eta) * beta2*pow(beta2+(r-s-v)*eta,r-s-v-1)/fact(r-s-v)*exp(-beta2-(r-s-v)*eta) * lambda*pow(lambda+(z-r+v-w)*eta,z-r+v-w-1)/fact(z-r+v-w)*exp(-lambda-(z-r+v-w)*eta) * zeta*pow(zeta+(y-s-v-w)*eta,y-s-v-w-1)/fact(y-s-v-w)*exp(-zeta-(y-s-v-w)*eta);
          } //end if
        } //for w
      } //for v
    } // for s

    help = (beta1+beta3) * pow(lambda*U*(1-alpha1-alpha3),2) * exp(-(beta1+beta3)-2*(lambda*U*(1-alpha1-alpha3))-y*eta-z*eta);
    int minyz  = std::min(y, z);

    int j;
    dwarf = 0.0;
    for (j=0;j <=minyz; j++) {
      dwarf = dwarf + pow(lambda*U*(1-alpha1-alpha3) + eta*(y-j),y-j-1) / fact(y-j) * pow(lambda*U*(1-alpha1-alpha3) + eta*(z-j),z-j-1) /fact(z-j) * pow(lambda*U*(alpha1+alpha3) + eta*j,j-1) /fact(j) * exp(j*eta);
    }

    pyz = dwarf*help;
    gryz = zwerg*pow(pyz,-1);

    nain = nain + gryz*lambda*pow(lambda+(x-r)*eta,x-r-1)*exp(-lambda-(x-r)*eta) /fact(x-r);
  } //end r

  return (nain);
}


// [[Rcpp::export]]
double dGP1h(int x, int y, double lambda, double alpha, double eta)
{

  double psi = eta*(1-alpha)/lambda;

  int minxy = std::min(x,y);

  int r;
  double nain = 0.0;
  for (r=0;r <=minxy; r++) {
    if (y >= r) {
      nain = nain + fact(y)/fact(y-r)/fact(r) * alpha *(1-alpha) * pow(alpha+psi*r,r-1) *pow(1-alpha+psi*(y-r),y-r-1)/pow(1+psi*y,y-1) *
        lambda*pow(lambda+(x-r)*eta,x-r-1)*exp(-lambda-(x-r)*eta)/fact(x-r);
    } //end if

  } //end r

  return(nain);
} // end of function



// [[Rcpp::export]]
double dR2(int r, int y, int z, double lambda, double alpha1, double alpha2, double alpha3, double eta)
{
  double U = 1/( 1 - alpha1-alpha2-alpha3);

  double beta1 = lambda * U *alpha1;
  double beta2 = lambda * U *alpha2;
  double beta3 = lambda * U *alpha3;
  double zeta = lambda*U*(1-2*alpha1-alpha3);
  double pyz;
  double gryz;
  double help;
  double dwarf;
  double zwerg;

  int s, v, w;
  zwerg=0.0;
  for ( s=0; s <= y; s++)
  {
    for ( v=0; v <= y; v++)
    {
      for ( w=0; w <= y; w++)
      {
        if (r-s-v >= 0 && z-r+v-w >= 0 && y-s-v-w >= 0){
          zwerg = zwerg  + beta3*pow(beta3+s*eta,s-1)/fact(s)*exp(-beta3-s*eta) * beta1*pow(beta1+v*eta,v-1)/fact(v)*exp(-beta1-v*eta) * beta1*pow(beta1+w*eta,w-1)/fact(w)*exp(-beta1-w*eta) * beta2*pow(beta2+(r-s-v)*eta,r-s-v-1)/fact(r-s-v)*exp(-beta2-(r-s-v)*eta) * lambda*pow(lambda+(z-r+v-w)*eta,z-r+v-w-1)/fact(z-r+v-w)*exp(-lambda-(z-r+v-w)*eta) * zeta*pow(zeta+(y-s-v-w)*eta,y-s-v-w-1)/fact(y-s-v-w)*exp(-zeta-(y-s-v-w)*eta);
        } //end if
      } //for w
    } //for v
  } // for s

  help = (beta1+beta3) * pow(lambda*U*(1-alpha1-alpha3),2) * exp(-(beta1+beta3)-2*(lambda*U*(1-alpha1-alpha3))-y*eta-z*eta);
  int minyz  = std::min(y, z);

  int j;
  dwarf = 0.0;
  for (j=0;j <=minyz; j++) {
    dwarf = dwarf + pow(lambda*U*(1-alpha1-alpha3) + eta*(y-j),y-j-1) / fact(y-j) * pow(lambda*U*(1-alpha1-alpha3) + eta*(z-j),z-j-1) /fact(z-j) * pow(lambda*U*(alpha1+alpha3) + eta*j,j-1) /fact(j) * exp(j*eta);
  }
  pyz = dwarf*help;
  gryz = zwerg*pow(pyz,-1);

  return (gryz);
}


// [[Rcpp::export]]
double Pyz(int y, int z, double lambda, double alpha1, double alpha2, double alpha3, double eta)
{
  double U = 1/( 1 - alpha1-alpha2-alpha3);

  double beta1 = lambda * U *alpha1;
  double beta2 = lambda * U *alpha2;
  double beta3 = lambda * U *alpha3;
  double pyz;
  double help;
  double dwarf;


  help = (beta1+beta3) * pow(lambda*U*(1-alpha1-alpha3),2) * exp(-(beta1+beta3)-2*(lambda*U*(1-alpha1-alpha3))-y*eta-z*eta);
  int minyz  = std::min(y, z);

  int j;
  dwarf = 0.0;
  for (j=0;j <=minyz; j++) {
    dwarf = dwarf + pow(lambda*U*(1-alpha1-alpha3) + eta*(y-j),y-j-1) / fact(y-j) * pow(lambda*U*(1-alpha1-alpha3) + eta*(z-j),z-j-1) /fact(z-j) * pow(lambda*U*(alpha1+alpha3) + eta*j,j-1) /fact(j) * exp(j*eta);
  }
  pyz = dwarf*help;


  return (pyz);
}

// [[Rcpp::export]]
std::vector<int> simGP2cov(double sumlimit,  double alpha1, double alpha2, double alpha3,
                           double eta, NumericVector lambdas, int T,int N, int seas1, int seas2, std::vector<int> data, NumericMatrix xreg,
                           std::vector<double> uniform, std::vector<int> innovations)
{
  double U = 1/( 1 - alpha1-alpha2-alpha3);
  double pyz;
  double gryz;
  double nain;
  double help;
  double dwarf;
  double zwerg;
  double lambda;

  int t;

  for (t = N + 1; t <= T;t++ ) {

    int y = data[t-seas1-1];
    int z = data[t-seas2-1];

    NumericVector covar = xreg(t-1,_);
    double mathelp = innerProduct(covar,lambdas);
    lambda = exp(mathelp);

    double beta1 = lambda * U *alpha1;
    double beta2 = lambda * U *alpha2;
    double beta3 = lambda * U *alpha3;
    double zeta = lambda*U*(1-2*alpha1-alpha3);

    double unif = uniform[t-1];

    int r = 0.0;
    double nain = 0.0;
    while (nain <= unif) {
      int s, v, w;
      zwerg=0.0;
      for ( s=0; s <= y; s++)
      {
        for ( v=0; v <= y; v++)
        {
          for ( w=0; w <= y; w++)
          {
            if (r-s-v >= 0 && z-r+v-w >= 0 && y-s-v-w >= 0){
              zwerg = zwerg  + beta3*pow(beta3+s*eta,s-1)/fact(s)*exp(-beta3-s*eta) * beta1*pow(beta1+v*eta,v-1)/fact(v)*exp(-beta1-v*eta) * beta1*pow(beta1+w*eta,w-1)/fact(w)*exp(-beta1-w*eta) * beta2*pow(beta2+(r-s-v)*eta,r-s-v-1)/fact(r-s-v)*exp(-beta2-(r-s-v)*eta) * lambda*pow(lambda+(z-r+v-w)*eta,z-r+v-w-1)/fact(z-r+v-w)*exp(-lambda-(z-r+v-w)*eta) * zeta*pow(zeta+(y-s-v-w)*eta,y-s-v-w-1)/fact(y-s-v-w)*exp(-zeta-(y-s-v-w)*eta);
            } //end if
          } //for w
        } //for v
      } // for s

      help = (beta1+beta3) * pow(lambda*U*(1-alpha1-alpha3),2) * exp(-(beta1+beta3)-2*(lambda*U*(1-alpha1-alpha3))-y*eta-z*eta);
      int minyz  = std::min(y, z);

      int j;
      dwarf = 0.0;
      for (j=0;j <=minyz; j++) {
        dwarf = dwarf + pow(lambda*U*(1-alpha1-alpha3) + eta*(y-j),y-j-1) / fact(y-j) * pow(lambda*U*(1-alpha1-alpha3) + eta*(z-j),z-j-1) /fact(z-j) * pow(lambda*U*(alpha1+alpha3) + eta*j,j-1) /fact(j) * exp(j*eta);
      }

      pyz = dwarf*help;
      gryz = zwerg*pow(pyz,-1);

      nain = nain + gryz;
      r = r+1;
    } //end r

    data[t-1] = r-1 + innovations[t-1];
  } //end t

  return (data);
} //end of function


// [[Rcpp::export]]
std::vector<int> simGP1cov(double sumlimit, double alpha, double eta,
                           NumericVector lambdas, int T,int N, int seas, std::vector<int> data, NumericMatrix xreg,
                           std::vector<double> uniform, std::vector<int> innovations)
{

  double lambda;

  int t;

  for (t = N + 1; t <= T; t++) {


    int y = data[t-seas-1];

    NumericVector covar = xreg(t-1,_);
    double mathelp = innerProduct(covar,lambdas);
    lambda = exp(mathelp);
    double psi = eta*(1-alpha)/lambda;

    double unif = uniform[t-1];
    int r = 0.0;
    double nain = 0.0;
    while (nain <= unif) {
      if (y >= r) {
        nain = nain + fact(y)/fact(y-r)/fact(r) * alpha *(1-alpha) * pow(alpha+psi*r,r-1) *pow(1-alpha+psi*(y-r),y-r-1)/pow(1+psi*y,y-1);
        r = r+1;
      } //end if

    } //end r
    data[t-1] = r-1 + innovations[t-1];

  } // end t

  return(data);
} // end of function


// [[Rcpp::export]]
double likelihoodGP2cov(double sumlimit, double alpha1, double alpha2, double alpha3,
                        double eta, NumericVector lambdas, int T, int seas1, int seas2, std::vector<int> data, NumericMatrix xreg)
{
  double U = 1/( 1 - alpha1-alpha2-alpha3);
  double pyz;
  double gryz;
  double nain;
  double help;
  double dwarf;
  double zwerg;
  double lambda;

  int t;
  double mlef = 0.0;
  for (t = seas2 + 1; t <= T;t++ ) {
    int x = data[t-1];
    int y = data[t-seas1-1];
    int z = data[t-seas2-1];

    NumericVector covar = xreg(t-1,_);
    double mathelp = innerProduct(covar,lambdas);
    lambda = exp(mathelp);

    double beta1 = lambda * U *alpha1;
    double beta2 = lambda * U *alpha2;
    double beta3 = lambda * U *alpha3;
    double zeta = lambda*U*(1-2*alpha1-alpha3);


    int minxyz  = std::min(x, y+z);
    int r;
    double nain = 0.0;
    for (r=0; r <= minxyz;r++) {
      int s, v, w;
      zwerg=0.0;
      for ( s=0; s <= y; s++)
      {
        for ( v=0; v <= y; v++)
        {
          for ( w=0; w <= y; w++)
          {
            if (r-s-v >= 0 && z-r+v-w >= 0 && y-s-v-w >= 0){
              zwerg = zwerg  + beta3*pow(beta3+s*eta,s-1)/fact(s)*exp(-beta3-s*eta) * beta1*pow(beta1+v*eta,v-1)/fact(v)*exp(-beta1-v*eta) * beta1*pow(beta1+w*eta,w-1)/fact(w)*exp(-beta1-w*eta) * beta2*pow(beta2+(r-s-v)*eta,r-s-v-1)/fact(r-s-v)*exp(-beta2-(r-s-v)*eta) * lambda*pow(lambda+(z-r+v-w)*eta,z-r+v-w-1)/fact(z-r+v-w)*exp(-lambda-(z-r+v-w)*eta) * zeta*pow(zeta+(y-s-v-w)*eta,y-s-v-w-1)/fact(y-s-v-w)*exp(-zeta-(y-s-v-w)*eta);
            } //end if
          } //for w
        } //for v
      } // for s

      help = (beta1+beta3) * pow(lambda*U*(1-alpha1-alpha3),2) * exp(-(beta1+beta3)-2*(lambda*U*(1-alpha1-alpha3))-y*eta-z*eta);
      int minyz  = std::min(y, z);

      int j;
      dwarf = 0.0;
      for (j=0;j <=minyz; j++) {
        dwarf = dwarf + pow(lambda*U*(1-alpha1-alpha3) + eta*(y-j),y-j-1) / fact(y-j) * pow(lambda*U*(1-alpha1-alpha3) + eta*(z-j),z-j-1) /fact(z-j) * pow(lambda*U*(alpha1+alpha3) + eta*j,j-1) /fact(j) * exp(j*eta);
      }

      pyz = dwarf*help;
      gryz = zwerg*pow(pyz,-1);

      nain = nain + gryz * lambda*pow(lambda+(x-r)*eta,x-r-1)*exp(-lambda-(x-r)*eta) /fact(x-r);
    } //end r
    
    if (nain <= 0){
      nain = 1 / pow(10, -12);
    }
    mlef = mlef + log(nain);
  } //end t

  return (-mlef);
} //end of function





// [[Rcpp::export]]
double likelihoodGP1cov(double sumlimit, double alpha, double eta, NumericVector lambdas, int T, int seas, std::vector<int> data, NumericMatrix xreg)
{

  double lambda;

  int t;

  double mlef = 0.0;
  for (t = seas + 1; t <= T; t++) {

    int x = data[t-1];
    int y = data[t-seas-1];

    NumericVector covar = xreg(t-1,_);
    double mathelp = innerProduct(covar,lambdas);
    lambda = exp( mathelp);
    double psi = eta*(1-alpha)/lambda;

    int minxy = std::min(x,y);

    int r;
    double nain = 0.0;
    for (r=0;r <=minxy; r++) {
      if (y >= r) {
        nain = nain + fact(y)/fact(y-r)/fact(r) * alpha *(1-alpha) * pow(alpha+psi*r,r-1) *pow(1-alpha+psi*(y-r),y-r-1)/pow(1+psi*y,y-1) *
          lambda*pow(lambda+(x-r)*eta,x-r-1)*exp(-lambda-(x-r)*eta)/fact(x-r);
      } //end if

    } //end r
    mlef = mlef + log(nain);

  } // end t

  return(-mlef);
} // end of function




  // [[Rcpp::export]]
double likelihoodGP2cov_repara(double sumlimit, double z1, double z2, double z3,
                          double eta, NumericVector lambdas, int T, int seas1, int seas2, std::vector<int> data, NumericMatrix xreg)
  {
    double alpha1, alpha2, alpha3;
    double U = 1/( 1 - alpha1-alpha2-alpha3);
    double pyz;
    double gryz;
    double nain;
    double help;
    double dwarf;
    double zwerg;
    double lambda;
    
    alpha3 = z1 / (1 + exp(z2));
    alpha1 = (z1 - alpha3) / 2;
    alpha2 = (1 - alpha1 - alpha3) / (1 + exp(z3));
    
    int t;
    double mlef = 0.0;
    for (t = seas2 + 1; t <= T;t++ ) {
      int x = data[t-1];
      int y = data[t-seas1-1];
      int z = data[t-seas2-1];
      
      NumericVector covar = xreg(t-1,_);
      double mathelp = innerProduct(covar,lambdas);
      lambda = exp(mathelp);
      
      double beta1 = lambda * U *alpha1;
      double beta2 = lambda * U *alpha2;
      double beta3 = lambda * U *alpha3;
      double zeta = lambda*U*(1-2*alpha1-alpha3);
      
      
      int minxyz  = std::min(x, y+z);
      int r;
      double nain = 0.0;
      for (r=0; r <= minxyz;r++) {
        int s, v, w;
        zwerg=0.0;
        for ( s=0; s <= y; s++)
        {
          for ( v=0; v <= y; v++)
          {
            for ( w=0; w <= y; w++)
            {
              if (r-s-v >= 0 && z-r+v-w >= 0 && y-s-v-w >= 0){
                zwerg = zwerg  + beta3*pow(beta3+s*eta,s-1)/fact(s)*exp(-beta3-s*eta) * beta1*pow(beta1+v*eta,v-1)/fact(v)*exp(-beta1-v*eta) *
                  beta1*pow(beta1+w*eta,w-1)/fact(w)*exp(-beta1-w*eta) * beta2*pow(beta2+(r-s-v)*eta,r-s-v-1)/fact(r-s-v)*exp(-beta2-(r-s-v)*eta) *
                  lambda*pow(lambda+(z-r+v-w)*eta,z-r+v-w-1)/fact(z-r+v-w)*exp(-lambda-(z-r+v-w)*eta) *
                  zeta*pow(zeta+(y-s-v-w)*eta,y-s-v-w-1)/fact(y-s-v-w)*exp(-zeta-(y-s-v-w)*eta);
              } //end if
            } //for w
          } //for v
        } // for s
        
        help = (beta1+beta3) * pow(lambda*U*(1-alpha1-alpha3),2) * exp(-(beta1+beta3)-2*(lambda*U*(1-alpha1-alpha3))-y*eta-z*eta);
        int minyz  = std::min(y, z);
        
        int j;
        dwarf = 0.0;
        for (j=0;j <=minyz; j++) {
          dwarf = dwarf + pow(lambda*U*(1-alpha1-alpha3) + eta*(y-j),y-j-1) / fact(y-j) * pow(lambda*U*(1-alpha1-alpha3) + eta*(z-j),z-j-1) /fact(z-j) * pow(lambda*U*(alpha1+alpha3) + eta*j,j-1) /fact(j) * exp(j*eta);
        }
        
        pyz = dwarf*help;
        gryz = zwerg*pow(pyz,-1);
        
        nain = nain + gryz*lambda*pow(lambda+(x-r)*eta,x-r-1)*exp(-lambda-(x-r)*eta) /fact(x-r);
      } //end r
      
      if (nain <= 0){
        nain = 1 / pow(10, -5);
      }
      
      mlef = mlef + log(nain);
    } //end t
    
    return (-mlef);
  } //end of function


  
  
  
    

