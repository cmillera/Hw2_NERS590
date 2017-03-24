#include <cmath>
#include <iostream>

#include "Estimator.h"
#include "Material.h"
#include "Particle.h"

void surface_current_estimator::score( particle* p ) { tally_hist += p->wgt(); }

void counting_estimator::score( particle* p ) { count_hist++; }

void track_length_estimator::score2( particle* p, double track_length, std::shared_ptr< material > medium)
{
    std::shared_ptr< cell > det_cell = p->cellPointer();
    double capxs = medium -> macro_capture_xs();
    double totalxs = medium -> macro_xs();
    
    tally_hist += p->wgt()*track_length*(capxs);
};

void counting_estimator::endHistory() {
  if ( tally.size() < count_hist + 1 ) { tally.resize( count_hist + 1, 0.0 ); }
  tally[ count_hist ] += 1.0;
  nhist++;
  count_hist = 0;
}

void counting_estimator::report() {
  std::cout << name() << std::endl;
  double s1 = 0.0, s2 = 0.0;
  for ( int i = 0 ; i < tally.size() ; i++ ) {
    double p = tally[i] / nhist;
    std::cout << i << " " << p << "   " << std::sqrt( p * ( 1.0 - p ) / nhist ) / p <<  std::endl;
    s1 += p * i;
    s2 += p * i * i;
  }
  std::cout << "   mean = " << s1 << std::endl;
  std::cout << "   var  = " << s2 - s1*s1 << std::endl;
}

void cell_pathLengthFlux_estimator::score2( particle* p, double path_length, std::shared_ptr< material > medium ) {
  
  double capxs = medium -> macro_capture_xs();
    double totalxs = medium -> macro_xs();
  
  tally_hist += p->wgt() * path_length*capxs;
}
