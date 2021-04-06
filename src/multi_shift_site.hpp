//
//  multi_shift_site.hpp
//  MGM
//
//  Created by Karl Gemayel on 7/4/20.
//  Copyright © 2020 Karl Gemayel. All rights reserved.
//

#ifndef multi_shift_site_hpp
#define multi_shift_site_hpp

#include <stdio.h>
#include "logger.h"
#include "common_2.h"
#include "site_2.h"

class MultiShiftSite {
public:
    MultiShiftSite();
    ~MultiShiftSite();
    
    
    // true when this class holds at least one valid site
    bool is_valid;
    
    // gets the site with maximum score. Score is computed as
    // shift_prior * motif_logodds * spacer_logodds
    Site* get_site_with_max_score(std::vector<unsigned char> const & nt,
                                  unsigned int pos,
                                  unsigned int status,
                                  bool include_shift=true) const;
    
    void add_site_with_shift(Site* site, float shift_prob);
    
    std::vector<Site *> get_sites() const;
    
private:
    // sifted sites
    std::vector< Site* > sites;
    
    // shift priors
    std::vector< float > shift_priors;
};

#endif /* multi_shift_site_hpp */
