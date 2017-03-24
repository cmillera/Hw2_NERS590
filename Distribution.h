#ifndef _DISTRIBUTION_HEADER_
#define _DISTRIBUTION_HEADER_

#include <memory>
#include <iostream>
#include <utility>
#include <vector>
#include <cmath>
#include <cassert>
#include <string>

#include "Random.h"
#include "Point.h"

template <class T>
class distribution {
  private:
    std::string distribution_name;
  protected:

  public:
     distribution( std::string label ) : distribution_name(label) {};          // constructor
    ~distribution() {};   // destructor

    virtual std::string name() final { return distribution_name; };
    virtual T sample() = 0;  // dummy function that must be implemented in each case
};

template <class T>
class arbitraryDelta_distribution : public distribution<T> {
  private:
    T result;
  public:
     arbitraryDelta_distribution( std::string label, T val ) : distribution<T>(label), result(val) {};
    ~arbitraryDelta_distribution() {};

    T sample() { return result; }
};

template <class T>
class arbitraryDiscrete_distribution : public distribution<T> {
  private:
     std::vector< std::pair< T, double > > cdf;
  public:
     arbitraryDiscrete_distribution( std::string label, std::vector< std::pair< T, double > > data );
    ~arbitraryDiscrete_distribution() {};
     T sample();
};

template < class T >
arbitraryDiscrete_distribution<T>::arbitraryDiscrete_distribution( std::string label, std::vector< std::pair< T, double > > data )
 : distribution<T>(label) {
  // convert pdf to cdf
  double c = 0.0;
  for ( auto d : data ) {
    // first is pointer to data type T and second is pdf input
    cdf.push_back( std::make_pair( d.first, d.second + c ) );
    c += d.second;
  }
}

template < class T >
T arbitraryDiscrete_distribution<T>::sample() {
  double   r = Urand() * cdf.back().second;
  for ( auto c : cdf ) {
    // first is pointer to data type T and second is cdf
    if ( r < c.second ) { return c.first; };
  }
  assert( false ); // should not get here
  return cdf.back().first; 
}

class delta_distribution : public distribution<double> {
  private:
    double a;
  public:
     delta_distribution( std::string label, double p1 ) : distribution(label), a(p1) {};
    ~delta_distribution() {};
    double sample() { return a; };
};

class uniform_distribution : public distribution<double> {
  private:
    double a, b;
  public:
     uniform_distribution( std::string label, double p1, double p2 ) : distribution(label), a(p1), b(p2) {};
    ~uniform_distribution() {};
    double sample();
};

class linear_aniso_distribution : public distribution<double> {
  private:
  public:
     linear_aniso_distribution( std::string label ) : distribution(label) {};
    ~linear_aniso_distribution() {};
    double sample();
};

class linear_distribution : public distribution<double> {
  private:
    double a, b, fa, fb;
  public:
    linear_distribution( std::string label, double x1, double x2, double y1, double y2 )  
      : distribution(label), a(x1), b(x2), fa(y1), fb(y2) {};
   ~linear_distribution() {};
   double sample();
};

class exponential_distribution : distribution<double> {
  private:
    double lambda;
  public:
     exponential_distribution( std::string label, double p1 ) : distribution(label), lambda(p1) {};
    ~exponential_distribution() {};
    double sample();
};

class normal_distribution : public distribution<double> {
  private:
    const double twopi = 2.0 * std::acos(-1.0);
    double mu, sigma;
  public:
     normal_distribution( std::string label, double p1, double p2 ) : distribution(label), mu(p1), sigma(p2) {};
    ~normal_distribution() {};
    double sample();
};

class HenyeyGreenstein_distribution : public distribution<double> {
  private:
    double a;
  public:
     HenyeyGreenstein_distribution( std::string label, double p1 ) : distribution(label), a(p1) {};
    ~HenyeyGreenstein_distribution() {};
    double sample();
};

class meanMultiplicity_distribution : public distribution<int> {
  private:
    double nu;
  public:
     meanMultiplicity_distribution( std::string label, double p1 ) : distribution(label), nu(p1) {};
    ~meanMultiplicity_distribution() {};
    int sample();
};

class TerrellFission_distribution : public distribution<int> {
  private:
    double nubar, sigma, b;

    std::vector<double> cdf;
  public:
    TerrellFission_distribution( std::string label, double p1, double p2, double p3 );
    ~TerrellFission_distribution() {};
    int sample();
};

class isotropicDirection_distribution : public distribution<point> {
  private:
    const double twopi = 2.0 * std::acos(-1.0);
  public:
     isotropicDirection_distribution( std::string label ) : distribution(label) {};
    ~isotropicDirection_distribution() {};
    point sample();
};

class anisotropicDirection_distribution : public distribution<point> {
  private:
    const double twopi = 2.0 * std::acos(-1.0);
    double sin_t;

    point axis;
    std::shared_ptr< distribution<double> > dist_mu;
  public:
     anisotropicDirection_distribution( std::string label, point p, std::shared_ptr< distribution<double> > dmu ) 
       : distribution(label), axis(p), dist_mu(dmu) 
       { axis.normalize(); sin_t = std::sqrt( 1.0 - axis.z * axis.z ); };
    ~anisotropicDirection_distribution() {};
    point sample();
};

class independentXYZ_distribution : public distribution<point> {
  private:
    std::shared_ptr< distribution<double> > dist_x, dist_y, dist_z;
  public:
     independentXYZ_distribution( std::string label, std::shared_ptr< distribution<double> > dx, 
       std::shared_ptr< distribution<double> > dy, std::shared_ptr< distribution<double> > dz ) 
       : distribution(label), dist_x(dx), dist_y(dy), dist_z(dz) {};
    ~independentXYZ_distribution() {};

    point sample();
};

class shell_distribution : public distribution<point> {
  private:
    const double twopi = 2.0 * std::acos(-1.0);
    double r_in, r_out;
  public:
     shell_distribution( std::string label, double p1, double p2 ) : distribution(label), r_in(p1), r_out(p2) {};
    ~shell_distribution() {};
    point sample();
};

class zbeam_distribution : public distribution<point>{
  private:
    const double twopi = 2.0 * std::acos(-1.0);
    double x0,y0,z0,Rad;
  public:
     zbeam_distribution( std::string label, double p1, double p2, double p3, double p4 ) 
        : distribution(label), x0(p1), y0(p2), z0(p3), Rad(p4) {};
    ~zbeam_distribution() {};
    point sample();
};
class threePoint_distribution : public distribution<point>{
  private:
    double x1,y1,z1,x2,y2,z2,x3,y3,z3,prob1,prob2,prob3;
  public:
     threePoint_distribution( std::string label, double p1, double p2, double p3, double p4, double p5,
     double p6, double p7, double p8, double p9, double p10, double p11, double p12) 
        : distribution(label), x1(p1), y1(p2), z1(p3), x2(p4), y2(p5), z2(p6), x3(p7), y3(p8), z3(p9),
        prob1(p10), prob2(p11), prob3(p12){};
    ~threePoint_distribution() {};
    point sample();
};
#endif
