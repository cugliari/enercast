#include <Rcpp.h>
using namespace Rcpp;

/*! \file   kerfWAV.c
    \brief  Prediction of time series with kernel smoothing
            and cross-validation for bandwidth choice using a 
            wavelet based dissimilarity.
            (see Antoniadis, Paparoditis and Sapatinas (2006, 2008))
    \author Jairo Cugliari (last modified february 2016)
*/

// --------------------------------------------------------- DistWav
/*! DistWav
   Distance between two series based on the difference (scale by scale)
   of the wavelets coefficients. 
  \param u             vector of difference between wavelet coefficients
  \return              the distance based on the wavelet coefficients
*/

// [[Rcpp::export]]
double DistWav(NumericVector u){
  int k = 0, i0 = 0, i1;
  int nn = u.size() + 1; 
  double dj = 0.0;
  double result = 0.0;

  int J = ceil(log(nn)/log(2));

  for(int j = 0; j < J; j++){
    nn = nn / 2 ;
    i1 = i0 + nn ;
    for(int i = i0; i < i1; i++) { 
      dj += u[k] * u[k]; 
      k++ ;
    }
    result = result + sqrt(dj / nn) ;
    dj =  0.0;
    i0 = i1;
  }
  return result;
}

// ------------------------------------------------------ prevkerfon
/*! Kernel estimation
Kernel estimation of a functional process based on the wavelet
representation using a given window width

\param serie_x   past values of the serie (1 --> n-1)
\param serie_y   past values of the serie (2 --> n)
\param serie_x0  current value of the functional serie (the starting
point)
\param kerneltype type of kernel to pass to evalKern
\param h         window width
\param n         number of observations in the functional serie serie_x
\param p_x       number of time-points for each observation of serie_x
(length of *serie_x0)
\param p_y       number of time-points for each observation of serie_y
\param result    vector containing the estimation K(serie_x0)
\return          the functional kernel estimation of *serie_x1 obtained
by the formula serie_x1=K(serie_x0) where K depends on
the past values of the serie
*/

// [[Rcpp::export]]
NumericVector prevkerfon (NumericVector serie_x, 
                          NumericVector serie_y, 
                          NumericVector serie_x0,
                          double h, double EPS){
  int p = serie_x0.length();
  int n = serie_x.length() / p;
  
  double kval = 0.0;
  NumericVector vide(p) ;
  NumericVector result(p);
  double sim, dist;
    
  for (int j= 0; j< p; j++)  result[j] = 0.0;
  
  for (int i = 0; i < n ; i++) // (*n - 1)
  {
    for (int j= 0; j< p; j++)  vide[j] = serie_x[(i * p) + j] - serie_x0[j];
    
    dist = DistWav(vide);
    sim = exp( - dist * dist / h / h) ;
    
    if(fabs(sim) < EPS) sim = 0.0;
  
    kval += sim;
    
    for (int j= 0; j< p; j++)  result[j] += serie_y[(i * p) + j] * sim;
    
    sim    = 0.0 ;
  }
  
  if( kval > 0.0) 
    for (int j= 0; j< p; j++)   result[j] /= kval;
  else 
    for (int j= 0; j< p; j++)   result[j] = 1 / n;
  
  return(result);
}


// ------------------------------------------------------- CVkerfon
/*! Cross validation of Functional Kernel
Cross validation of the window width used in the functional kernel
estimation of an functional times series process

\param serie_x   past values of the serie (1 --> n-1)
\param serie_y   past values of the serie (2 --> n)
\param h         window width
\param n         number of obs in the  serie_x
\param p_x       number of time-points for each observation of serie_x
(length of *serie_x0)
\param p_y       number of time-points for each observation of serie_y
\param r         length of the cross-validation period
\param result    the criteria of the cross validation
\result          the sum of the square errors obtain with a width of h
*/

// [[Rcpp::export]]
double CVkerfon (NumericVector serie_x, 
                 NumericVector serie_y, 
                 int p,
                 int r, 
                 double h, double EPS){
  
  int n = serie_x.length() / p;
  int tj= n -r; 

  NumericVector vide_court(p); 
  NumericVector res(p);
  NumericVector vide_long(tj * p);
  NumericVector vide_long2(tj * p);
  
  double result = 0.0 ;

  for (int k = 0; k < (tj * p); k++) vide_long[k]  = serie_x[k];
  for (int k = 0; k < (tj * p); k++) vide_long2[k] = serie_y[k];
    
  for (int k = 0; k < r; k++){
    for (int i = 0; i < p; i++)
      vide_court[i] = serie_x[((k + n - r) * p) + i];
      
  res = prevkerfon(vide_long, vide_long2, vide_court, h, EPS);

  for (int i = 0; i < p; i++)
        result += (res[i] - serie_y[((k + n - r) * p) + i]) * 
                  (res[i] - serie_y[((k + n - r) * p) + i]);
    }
    
  result = sqrt( result / (n - r) ) ;

  return(result);
}





/*** R
#DistWav(1:31)
x  <- sin(1:(30*15))
y  <- sin((30*15):1)
x0 <- sin(1:15)
#prevkerfon(x, y, x0, 1,EPS = .Machine$double.eps)
CVkerfon(x, y, p = 15, r = 10, h = 1, EPS = .Machine$double.eps)
*/

