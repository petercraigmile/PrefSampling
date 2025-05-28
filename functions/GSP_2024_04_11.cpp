
#include <Rcpp.h>
using namespace Rcpp;



// [[Rcpp::export]]
NumericMatrix cppEdist (NumericMatrix sites1, NumericMatrix sites2) {

    int n1 = sites1.nrow();
    int n2 = sites2.nrow();
    int p  = sites1.ncol();    
    int i, j, k;
    double ss, diff;
    
    NumericMatrix D(n1, n2);

    for (i=0; i<n1; i++)
      for (j=0; j<n2; j++) {

	ss = 0.0;

	for (k=0; k<p; k++) {

	  diff = sites1(i, k) - sites2(j, k);
	    
	  ss += diff * diff;
	}

	D(i, j) = std::sqrt(ss);
      }

    return D;
}



// [[Rcpp::export]]
NumericMatrix cppGSPexp (NumericMatrix dists,
			 NumericVector theta,
			 double eps=1e-12) {

  int p = theta.size();
  int m = dists.nrow();
  int n = dists.ncol();
  int i, j;
  NumericMatrix Cov(m, n);
  double tau2 = theta(0);
  double phi  = theta(1);

  if (p==2) { // no nugget
    
    for (i=0; i<m; i++) {
      for (j=0; j<n; j++) {
	
	Cov(i,j) = tau2 * std::exp(-dists(i,j)/phi);	
      }
    }
    
  } else { // includes nugget
    
    double sigma2 = theta(2);
    
    for (i=0; i<m; i++) {
      for (j=0; j<n; j++) {	
	
	if (dists(i,j) > eps) {
	  
	  Cov(i,j) = tau2 * std::exp(-dists(i,j)/phi);
	  
	} else {
	  
	  Cov(i,j) = tau2 * std::exp(-dists(i,j)/phi) + sigma2;
	}      
      }
    }    
  }
  
  return (Cov);
}



// [[Rcpp::export]]
NumericMatrix cppGSPDiggle (NumericMatrix dists,
			 NumericVector theta,
			 double eps=1e-12) {

  int p = theta.size();
  int m = dists.nrow();
  int n = dists.ncol();
  int i, j;
  NumericMatrix Cov(m, n);
  double tau2 = theta(0);
  double phi  = theta(1);
  double scaled_h;

  if (p==2) { // no nugget
    
    for (i=0; i<m; i++) {
      for (j=0; j<n; j++) {

	scaled_h = dists(i,j) / phi;
	
	Cov(i,j) = tau2 * (1.0 + scaled_h) * std::exp(-scaled_h);	
      }
    }
    
  } else { // includes nugget
    
    double sigma2 = theta(2);
    
    for (i=0; i<m; i++) {
      for (j=0; j<n; j++) {

	scaled_h = dists(i,j) / phi;

	if (dists(i,j) > eps) {
	  
	  Cov(i,j) = tau2 * (1.0 + scaled_h) * std::exp(-scaled_h);
	  
	} else {
	  
	  Cov(i,j) = tau2 * (1.0 + scaled_h) * std::exp(-scaled_h) + sigma2;
	}      
      }
    }    
  }
  
  return (Cov);
}





// [[Rcpp::export]]
NumericMatrix cppGSPGaussian (NumericMatrix dists,
			      NumericVector theta,
			      double eps=1e-12) {

  int p = theta.size();
  int m = dists.nrow();
  int n = dists.ncol();
  int i, j;
  NumericMatrix Cov(m, n);
  double tau2 = theta(0);
  double phi  = theta(1);
  double scaled_h;

  if (p==2) { // no nugget
    
    for (i=0; i<m; i++) {
      for (j=0; j<n; j++) {

	scaled_h = dists(i,j)/phi;		
	
	Cov(i,j) = tau2 * std::exp(-scaled_h*scaled_h);	
      }
    }
    
  } else { // includes nugget
    
    double sigma2 = theta(2);
    
    for (i=0; i<m; i++) {
      for (j=0; j<n; j++) {

	scaled_h = dists(i,j)/phi;
	
	if (dists(i,j) > eps) {

	  Cov(i,j) = tau2 * std::exp(-scaled_h*scaled_h);
	  
	} else {
	  
	  Cov(i,j) = tau2 * std::exp(-scaled_h*scaled_h) + sigma2;
	}      
      }
    }    
  }
  
  return (Cov);
}





// [[Rcpp::export]]
NumericMatrix cppGSPpowerexp (NumericMatrix dists,
			      NumericVector theta,
			      double eps=1e-12) {

  int p = theta.size();
  int m = dists.nrow();
  int n = dists.ncol();
  int i, j;
  NumericMatrix Cov(m, n);
  double tau2 = theta(0);
  double phi  = theta(1);
  double expo = theta(2);
  double scaled_h;

  if (p==3) { // no nugget
    
    for (i=0; i<m; i++) {
      for (j=0; j<n; j++) {

	scaled_h = dists(i,j)/phi;		
	
	Cov(i,j) = tau2 * std::exp(-std::pow(scaled_h, expo));
      }
    }
    
  } else { // includes nugget
    
    double sigma2 = theta(3);
    
    for (i=0; i<m; i++) {
      for (j=0; j<n; j++) {

	scaled_h = dists(i,j)/phi;
	
	if (dists(i,j) > eps) {

	  Cov(i,j) = tau2 * std::exp(-std::pow(scaled_h, expo));
	  
	} else {
	  
	  Cov(i,j) = tau2 * std::exp(-std::pow(scaled_h, expo)) + sigma2;
	}      
      }
    }    
  }
  
  return (Cov);
}



