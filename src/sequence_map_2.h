// ====================================================
// Author: Alex Lomsadze
// Copyright: Georgia Institute of Technology
// Copyright: GeneProbe Inc.
// Release: 2018
// File: sequence_map_2.h
// Project: GeneMark.hmm-2 (no introns)
// ====================================================

#ifndef SEQUENCE_MAP_2_H
#define SEQUENCE_MAP_2_H

#include <string>
#include <vector>
#include <set>
#include <map>

#include "common_2.h"
#include "logger.h"
#include "model_2.h"
#include "data_2.h"
#include "exit.h"

// ----------------------------------------------------
class MapValue
{
public:
    MapValue()
    {
        is_good = false;
        status = 0;
        pos = -1;
        end_pos = -1;
        gc = -1;
        next = 0;
        prev = 0;
        gtype = 0;
        logodd_dur = 0;
        logodd_orf = 0;
        logodd_start = 0;
        prob = LOG_ZERO;
        stype = 0;
        vprob = LOG_ZERO;
        path = 0;

        IsInforced = false;

        logP = LOG_ZERO;
        logP_type = 0;

        logodd_total_native   = LOG_ZERO;
        logodd_total_atypical = LOG_ZERO;
        logodd_orf_native     = LOG_ZERO;
        logodd_orf_atypical   = LOG_ZERO;
        logodd_dur_native     = LOG_ZERO;
        logodd_dur_atypical   = LOG_ZERO;
        logodd_tr_native      = LOG_ZERO;
        logodd_tr_atypical    = LOG_ZERO;

        logodd_start_stop     = LOG_ZERO;
        logodd_RBS            = LOG_ZERO;
        logodd_Promoter       = LOG_ZERO;
        logodd_RBS_SC         = LOG_ZERO;
        logodd_Promoter_SC    = LOG_ZERO;
        
        // break up into motif/spacer
        logodd_RBS_motif      = LOG_ZERO;
        logodd_RBS_spacer      = LOG_ZERO;
        logodd_Promoter_motif      = LOG_ZERO;
        logodd_Promoter_spacer      = LOG_ZERO;
        logodd_RBS_max_spacer = LOG_ZERO;
        logodd_Promoter_max_spacer = LOG_ZERO;

        bayes = 0;
    }

    ~MapValue(){};

    bool is_good;
    unsigned int  status;
    unsigned int  pos;
    unsigned int  end_pos;
    unsigned int  gc;
    MapValue*     next;
    MapValue*     prev;
    char          gtype;
    double        logodd_dur;
    double        logodd_orf;
    double        logodd_start;
    double        prob;
    int           stype;
    double          vprob;
    MapValue*     path;

    bool IsInforced;

    double  logP;
    char    logP_type;

    double  logodd_total_native;
    double  logodd_total_atypical;
    double  logodd_orf_native;
    double  logodd_orf_atypical;
    double  logodd_dur_native;
    double  logodd_dur_atypical;
    double  logodd_tr_native;
    double  logodd_tr_atypical;

    double  logodd_start_stop;
    double  logodd_RBS;
    double  logodd_Promoter;
    double  logodd_RBS_SC;
    double  logodd_Promoter_SC;
    
    
    double logodd_RBS_motif      ;
    double logodd_RBS_spacer     ;
    double logodd_Promoter_motif ;
    double logodd_Promoter_spacer;
    double logodd_RBS_max_spacer;
    double logodd_Promoter_max_spacer;

    double bayes;

private:
    ;
};
// ----------------------------------------------------
class BestValue
{
public:
    BestValue(){ L=0; R=0; complete_L=false; complete_R=false; strand=STRAND_TYPES::NON; gtype=0; origin=0; };
    ~BestValue(){}

    unsigned int  L;
    unsigned int  R;
    bool  complete_L;
    bool  complete_R;
    STRAND_TYPE  strand;
    char  gtype;
    char  stype;
    unsigned int id;

    SiteInfo si;

    MapValue*  origin;

private:
    ;
};
// ----------------------------------------------------
class SequenceMap
{
public:
    SequenceMap( Settings & settings, Logger * const logger );
    ~SequenceMap(){};

    std::vector< MapValue >  data;
    std::vector< BestValue > predictions;
    
    GMS2_GROUP gms2_group;

    void Init( std::vector<unsigned int> const & flags, unsigned int const min_gene_length );
    void CalcGC( Data & d );
    void CalcLogP( std::vector< Model* > & mod, std::vector<unsigned char> & nt, char const gtype );
    void CalcLogodd( std::vector< Model* > & mod, std::vector<unsigned char> & nt, char const gtype, GMS2_GROUP gms2_group = GMS2_NONE);
    void AddIniTermProb(double p_ini_non);

    void Run(bool best_start_before_dp, double delta);

    void CalcStartScoreForPositionNative(Model* m, std::vector<unsigned char> & nt, std::vector< MapValue >::iterator itr, GMS2_GROUP gms2_group=GMS2_NONE);
    void CalcStartScoreForPositionAtypical(Model* m, std::vector<unsigned char> & nt, std::vector< MapValue >::iterator itr, GMS2_GROUP gms2_group=GMS2_NONE);

    void CalcStarts( Model* mod, std::vector<unsigned char> & nt, GMS2_GROUP gms2_group=GMS2_NONE);
    void CalcStartsGC(std::vector< Model* > & mod, std::vector<unsigned char> & nt, GMS2_GROUP gms2_group=GMS2_NONE);
    void AddSiteInfo( Model* mod, std::vector<unsigned char> & nt );
    void AddSiteInfoAtypical( std::vector< Model* > & mod, std::vector<unsigned char> & nt, GMS2_GROUP gms2_group);

    void AddCodingEvidence( std::map<int, INTERVAL_EVI > & d, std::map<int, INTERVAL_EVI > & r );

    double final_logodd;
    

private:

    void FindMaxScoredORFs(void);
    void Dynamic(void);
    void BestPath(void);
    void TestORFsDir(void);
    void TestORFsRev(void);

    void SetBaseData(std::vector<unsigned int> const & flags);
    void SetInitAndTerm(unsigned int const seq_length);

    void LinkDirORF(unsigned int const frame);
    void LinkRevORF(unsigned int const frame);

    void MarkValidDirORFs(unsigned int const min_length);
    void MarkValidRevORFs(unsigned int const min_length);

    std::string PrintMap(void);

    void FindMaxScoredDirORFs(void);
    void FindMaxScoredRevORFs(void);

    bool FindIsland( std::vector< MapValue >::iterator & L, std::vector< MapValue >::iterator & R );
    void ParseIsland( std::vector< MapValue >::iterator L, std::vector< MapValue >::iterator R );

    void SelectBestN( std::vector< MapValue >::iterator p, std::vector< MapValue >::iterator c );
    
    void SelectBestCodDir( MapValue* p, std::vector< MapValue >::iterator c );
    void SelectBestCodRev( MapValue* p, std::vector< MapValue >::iterator c );

    void SelectBestOverlapCodDir( MapValue* p, std::vector< MapValue >::iterator o, std::vector< MapValue >::iterator c );
    void SelectBestOverlapCodRev( MapValue* p, std::vector< MapValue >::iterator o, std::vector< MapValue >::iterator c );

    void DP_initialization_of_noncoding_state(std::vector< MapValue >::iterator itr);
    void DP_initialization_of_coding_state(std::vector< MapValue >::iterator itr);
    void DP_ORF(std::vector< MapValue >::iterator itr);
    void DP_intergenic(std::vector< MapValue >::iterator itr);

    void ExcludeNegativeScoredDirORFs(double delta);
    void ExcludeNegativeScoredRevORFs(double delta);

    void CalcBayesDirORFs(void);
    void CalcBayesRevORFs(void);

    // settings

    double overlap_penalty;
    double noncoding_state_initiation_and_termination_probability;

    Logger * const logger;
};
//================================================
#endif // SEQUENCE_MAP_2_H

