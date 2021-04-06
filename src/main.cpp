// ====================================================
// Author: Alex Lomsadze
// 
// Copyright: Georgia Institute of Technology
// Copyright: GeneProbe Inc.
// 
// Contact information:
//     Mark Borodovsky  borodovsky@gatech.edu
//     Alex Lomsadze    alexl@gatech.edu
// 
// GeneMark.hmm-2 algorithm version 1.0 was published in
// Genome Research 2018
// https://www.ncbi.nlm.nih.gov/pubmed/29773659
// 
// "Modeling leaderless transcription and atypical genes
// results in more accurate gene prediction in prokaryotes."
// 
// by Alexandre Lomsadze, Karl Gemayel, Shiyuyun Tang and Mark Borodovsky
// 
// Project: GeneMark.hmm-2 (no introns)
// File: main.cpp
//
// Projects VERSION is set by first declaration in the main() function
// ====================================================

#include <string>
#include <iostream>

#include "exit.h"
#include "logger.h"
#include "settings_2.h"
#include "sequence_file_2.h"
#include "parameters_2.h"
#include "pset_2.h"
#include "data_2.h"
#include "sequence_map_2.h"
#include "evidence_2.h"
#include "model_2.h"
#include "output_2.h"

#ifdef LICENSE
#include "check.h"
#endif

void set_correct_gms2_start_models(Pset &pset, GMS2_GROUP gms2_group, int bac_arc) {
    
    // archaea
    if (bac_arc == 1) {
        
        // group A
        if (gms2_group == GMS2_A) {
            // set bacteria RBS_A as archaea RBS_A
            for (size_t i = 0; i < pset.first.size(); i++) {
                pset.first[i]->RBS_A = pset.second[i]->RBS_A;
                pset.first[i]->SC_RBS_A = pset.second[i]->SC_RBS_A;
                
                pset.first[i]->pATG_A = pset.second[i]->pATG_A;
                pset.first[i]->pGTG_A = pset.second[i]->pGTG_A;
                pset.first[i]->pTTG_A = pset.second[i]->pTTG_A;
            }
        }
        // group D
        else if (gms2_group == GMS2_D) {
            // set bacteria RBS_A as archaea RBS_A
            for (size_t i = 0; i < pset.first.size(); i++) {
                pset.first[i]->RBS_D = pset.second[i]->RBS_D;
                pset.first[i]->SC_RBS_D = pset.second[i]->SC_RBS_D;
                pset.first[i]->PROMOTER_D = pset.second[i]->PROMOTER_D;
                pset.first[i]->SC_PROMOTER_D = pset.second[i]->SC_PROMOTER_D;
                
                pset.first[i]->pATG_D = pset.second[i]->pATG_D;
                pset.first[i]->pGTG_D = pset.second[i]->pGTG_D;
                pset.first[i]->pTTG_D = pset.second[i]->pTTG_D;
            }
        }
    }
    // bacteria
    if (bac_arc == 0) {
        // group A
        if (gms2_group == GMS2_A) {
            for (size_t i = 0; i < pset.first.size(); i++) {
                pset.second[i]->RBS_A = pset.first[i]->RBS_A;
                pset.second[i]->SC_RBS_A = pset.first[i]->SC_RBS_A;
                
                pset.second[i]->pATG_A = pset.first[i]->pATG_A;
                pset.second[i]->pGTG_A = pset.first[i]->pGTG_A;
                pset.second[i]->pTTG_A = pset.first[i]->pTTG_A;
                
//                pset.second[i]->PROMOTER_C = pset.first[i]->PROMOTER_C;
//                pset.second[i]->SC_PROMOTER_C = pset.first[i]->SC_PROMOTER_C;
            }
        }
        // group B
        if (gms2_group == GMS2_B) {
           for (size_t i = 0; i < pset.first.size(); i++) {
               pset.second[i]->RBS_B = pset.first[i]->RBS_B;
               pset.second[i]->SC_RBS_B = pset.first[i]->SC_RBS_B;
               
               pset.second[i]->pATG_B = pset.first[i]->pATG_B;
               pset.second[i]->pGTG_B = pset.first[i]->pGTG_B;
               pset.second[i]->pTTG_B = pset.first[i]->pTTG_B;
           }
       }
        // group C
        else if (gms2_group == GMS2_C) {
            for (size_t i = 0; i < pset.first.size(); i++) {
                pset.second[i]->RBS_C = pset.first[i]->RBS_C;
                pset.second[i]->SC_RBS_C = pset.first[i]->SC_RBS_C;
                pset.second[i]->PROMOTER_C = pset.first[i]->PROMOTER_C;
                pset.second[i]->SC_PROMOTER_C = pset.first[i]->SC_PROMOTER_C;
                
                pset.second[i]->pATG_C = pset.first[i]->pATG_C;
                pset.second[i]->pGTG_C = pset.first[i]->pGTG_C;
                pset.second[i]->pTTG_C = pset.first[i]->pTTG_C;
            }
        }
        // group X
          if (gms2_group == GMS2_B) {
            for (size_t i = 0; i < pset.first.size(); i++) {
                pset.second[i]->RBS_X = pset.first[i]->RBS_X;
                pset.second[i]->SC_RBS_X = pset.first[i]->SC_RBS_X;
                
                pset.second[i]->pATG_X = pset.first[i]->pATG_X;
                pset.second[i]->pGTG_X = pset.first[i]->pGTG_X;
                pset.second[i]->pTTG_X = pset.first[i]->pTTG_X;
            }
        }
        // group C
        else if (gms2_group == GMS2_AC) {
            for (size_t i = 0; i < pset.first.size(); i++) {
                pset.second[i]->RBS_C = pset.first[i]->RBS_C;
                pset.second[i]->SC_RBS_C = pset.first[i]->SC_RBS_C;
                pset.second[i]->PROMOTER_C = pset.first[i]->PROMOTER_C;
                pset.second[i]->SC_PROMOTER_C = pset.first[i]->SC_PROMOTER_C;
                
                pset.second[i]->pATG_C = pset.first[i]->pATG_C;
                pset.second[i]->pGTG_C = pset.first[i]->pGTG_C;
                pset.second[i]->pTTG_C = pset.first[i]->pTTG_C;
                
                pset.second[i]->RBS_A = pset.first[i]->RBS_A;
                pset.second[i]->SC_RBS_A = pset.first[i]->SC_RBS_A;
                
                pset.second[i]->pATG_A = pset.first[i]->pATG_A;
                pset.second[i]->pGTG_A = pset.first[i]->pGTG_A;
                pset.second[i]->pTTG_A = pset.first[i]->pTTG_A;
            }
        }
    }
}

float compute_logodds_and_fill_in_seqmap(Pset &pset_original, Data &data, SequenceMap& seqmap, Settings &settings, GMS2_GROUP group, int bac_arc) {
    
    // copy parameter set and turn on/off appropriate start models based on gms2 group
    Pset pset = pset_original.deepCopy();
    set_correct_gms2_start_models(pset, group, bac_arc);
    
    seqmap.Init(data.flag, pset.min_gene_length);
    seqmap.CalcGC(data);

//    if (evidence.data.size())
//        seqmap.AddCodingEvidence(data.evi_dir_orf, data.evi_rev_orf);
    
    if (pset.native[0])
        seqmap.CalcStarts(pset.native[0], data.nt, group);
    else
    {
        if (bac_arc == 0) {
            seqmap.CalcStartsGC(pset.first, data.nt, group);
        }
        else
            seqmap.CalcStartsGC(pset.second, data.nt, group);
    }

    if (pset.native[0])
    {
        seqmap.CalcLogP(pset.native, data.nt, NATIVE_TYPE);
        seqmap.CalcLogodd(pset.native, data.nt, NATIVE_TYPE, group);
    }

    if (pset.first[0])
    {
        seqmap.CalcLogP(pset.first, data.nt, ATYPICAL_TYPE_1);
        seqmap.CalcLogodd(pset.first, data.nt, ATYPICAL_TYPE_1, group);
    }

    if (pset.second[0])
    {
        seqmap.CalcLogP(pset.second, data.nt, ATYPICAL_TYPE_2);
        seqmap.CalcLogodd(pset.second, data.nt, ATYPICAL_TYPE_2, group);
    }

    seqmap.Run(settings.hmm.best_start_before_dp, settings.hmm.delta);

    if (pset.native[0])
        seqmap.AddSiteInfo(pset.native[0], data.nt);
    
    if (bac_arc == 0)
        seqmap.AddSiteInfoAtypical(pset.first, data.nt, group);
    else
        seqmap.AddSiteInfoAtypical(pset.second, data.nt, group);
    
    return seqmap.final_logodd;
}


int get_most_common_type(std::vector< BestValue > & predictions) {
    
    size_t count_bac = 0, count_arc = 0;
    for( std::vector< BestValue >::reverse_iterator itr = predictions.rbegin(); itr != predictions.rend(); ++itr )
    {

        switch (itr->gtype & 7) {
            case ATYPICAL_TYPE_1:
                count_bac += 1;
                break;
            case ATYPICAL_TYPE_2:
                count_arc += 1;
                break;
        }
    }
    
    if (count_bac >= count_arc)
        return 0;
    else
        return 1;
}

// ----------------------------------------------------
int main( int argc, char** argv )
{
	try
	{
		std::string VERSION = "1.23";

#ifdef LICENSE
		VERSION += "_lic";
#endif

		Logger        logger;
		Settings      settings( argc, argv, &logger, VERSION );

#ifdef LICENSE
		// to do:
		// move to new version of OS independent key implementation
		// function should return true/false and support verbose mode
		// function should be called from settings after parsing of parameters and before usage

		char path_to_key[] = "";
		check_timekey(path_to_key);
#endif

		SequenceFile  sequence( settings, &logger );
		Evidence      evidence( settings, &logger );

		if ( sequence.AllAtOnce() )
			evidence.SyncWithSequence( sequence.data );

		Parameters    parameters( settings, &logger );
		Pset          pset( settings, parameters, &logger );

		parameters.Initialize( parameters.model_set );

		Data          data( settings, pset, &logger );
		SequenceMap   seqmap( settings, &logger );
		Output        output( settings, &logger, VERSION );

		output.evi_dir_orf = & data.evi_dir_orf;
		output.evi_rev_orf = & data.evi_rev_orf;

		output.Header1(pset);

		bool run_on_all_contigs = true;
		if (!settings.in.run_on_records.empty() || !settings.in.ignore_records.empty())
			run_on_all_contigs = false;

		// progress bar variables
		int records_in = 0;
		int current_record = 0;
		int progress_bar_at = 0;

		if (settings.progress_bar)
			records_in = sequence.RecordsIn();

		for( FastaVectorItr itr = sequence.data.begin() ; itr != sequence.data.end(); itr = sequence.Next( itr ) )
		{
			if (settings.progress_bar)
			{
				if (current_record == progress_bar_at)
				{
					int progress_percent = int(100 * progress_bar_at / records_in);
					std::cout << "Progress:" << progress_percent << "%" << std::endl;
					progress_bar_at += int(records_in / 10);
				}

				++current_record;
			}

			if (!run_on_all_contigs)
			{
				if (!sequence.RunOnThis(itr->name))
					continue;

				if (sequence.IgnoreThis(itr->name))
					continue;
			}

			output.Header2( itr );

			// skip short sequences
			if (itr->data.size() < pset.min_gene_length)
				continue;

			do
			{
				data.Set(itr->data);

				if (evidence.data.size())
					data.ApplyEvidence(evidence.data, itr->name);

                seqmap.Init(data.flag, pset.min_gene_length);
                seqmap.CalcGC(data);

                // test multiple groups
                GMS2_GROUP best_group = GMS2_NONE;
                char best_label = 'N';
                float best_score = -10000000000;
                GMS2_GROUP all_groups [] = {
                    GMS2_A, GMS2_B, GMS2_C, GMS2_D, GMS2_X, GMS2_AC
                };
                char group_labels []  = {'A', 'B', 'C', 'D', 'X', 'Z'};
                int best_type = 0;
                int num_groups = 5;
                
                // determine if archaea or bacteria
                float score = compute_logodds_and_fill_in_seqmap(pset, data, seqmap, settings, GMS2_NONE, 0);
                //std::cout << "Score NONE: " << score << std::endl;
                int genome_type = get_most_common_type(seqmap.predictions);

                for (int bac_arc = 0; bac_arc < 2; bac_arc+=1) {
                    if (bac_arc != genome_type) {
                        continue;
                    }
                    for (int group_idx = 0; group_idx < num_groups; group_idx++) {
                        
                        // bacteria cannot be group D
                        if (bac_arc == 0 && (all_groups[group_idx] == GMS2_D))
                            continue;
                        // archaea only groups A and D
                        else if (bac_arc == 1 && (all_groups[group_idx] != GMS2_A && all_groups[group_idx] != GMS2_D && all_groups[group_idx] != GMS2_NONE))
                            continue;
                        
                        float current_score = compute_logodds_and_fill_in_seqmap(pset, data, seqmap, settings, all_groups[group_idx], bac_arc);
                        if (current_score > best_score) {
                            best_score = current_score;
                            best_group = all_groups[group_idx];
                            best_label = group_labels[group_idx];
                            best_type = bac_arc;
                        }
                        //std::cout << (bac_arc == 0 ? "Bacteria" : "Archaea  ") << "\t" << group_labels[group_idx] << "\t" << current_score << std::endl;
                        
                        // // output to temp file
                        // Settings      settings_tmp( argc, argv, &logger, VERSION );
                        // std::stringstream s;
                        // s << settings_tmp.out.filename << "_" << bac_arc << "_" << group_labels[group_idx];
                        // settings_tmp.out.filename = s.str();
                        // 
                        // Output        output_tmp( settings_tmp, &logger, VERSION );
                        // output_tmp.Header1(pset);
                        // output_tmp.Header2( itr );
                        // output_tmp.PrintGenes(seqmap.predictions, itr, pset.genetic_code);
                        // output_tmp.Stat(itr, seqmap.final_logodd, seqmap.predictions);
                        // output_tmp.Footer();

//                        output_tmp.PrintGenes(seqmap.predictions, itr, pset.genetic_code);

                    }
                }
                
                //std::cout << (best_type == 0 ? "Bacteria" : "Archaea") << "\t" << "Best group: " << best_label << "\t" << best_score << std::endl;
                
                // rerun with best group
                compute_logodds_and_fill_in_seqmap(pset, data, seqmap, settings, best_group, best_type);

			} while (0);

            
			output.PrintGenes(seqmap.predictions, itr, pset.genetic_code);
			output.Stat(itr, seqmap.final_logodd, seqmap.predictions);
            
		}

		output.Footer();

		logger.Print( 0, "# GeneMark.hmm ... done" );
	}
	catch( std::bad_alloc const & )
	{
		Exit("error, out of memory");
	}
	catch( std::exception const & e )
	{
		Exit("error,", e.what());
	}

	return 0;
};
// ----------------------------------------------------

