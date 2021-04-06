//
//  multi_shift_site.cpp
//  MGM
//
//  Created by Karl Gemayel on 7/4/20.
//  Copyright © 2020 Karl Gemayel. All rights reserved.
//

#include "multi_shift_site.hpp"

MultiShiftSite::MultiShiftSite() {
    is_valid = false;
}

MultiShiftSite::~MultiShiftSite() {
    // nothing
//    for (std::vector<Site*>::iterator itr = sites.begin(); itr != sites.end(); itr++) {
//        delete *itr;
//    }
//    
//    sites.clear();
}

void MultiShiftSite::add_site_with_shift(Site *site, float shift_prob) {
    
    if (site == NULL) {
        // should raise warning
        return;
    }
    
    shift_priors.push_back(shift_prob);
    sites.push_back(site);
}

Site* MultiShiftSite::get_site_with_max_score(std::vector<unsigned char> const & nt,
                                              unsigned int pos,
                                              unsigned int status,
                                              bool include_shift) const {
    
    size_t num_shifts = shift_priors.size();
    
    float max_score = -100000;
    Site* site_with_max = NULL;
    
    // for each shift, compute score and keep track of max
    for (size_t shift = 0; shift < num_shifts; shift++) {
        Site* curr_site = sites[shift];
        
        float curr_score = curr_site->GetWithDur(nt, pos, status);
        
        if (include_shift) {
            curr_score *= shift_priors[shift];
        }
        
        if (curr_score > max_score) {
            max_score = curr_score;
            site_with_max = curr_site;
        }
    }
    
    return site_with_max;
}

std::vector<Site *> MultiShiftSite::get_sites() const {
    return sites;
}
