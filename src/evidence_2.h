// ====================================================
// Author: Alex Lomsadze
// Copyright: Georgia Institute of Technology
// Copyright: GeneProbe Inc.
// File: evidence_2.h
// Code was tested with boost 1.48 <boost/algorithm/string.hpp>
// Project: GeneMark.hmm-2 (no introns)
// ====================================================

#ifndef GMHMMP2_EVIDENCE_2_H
#define GMHMMP2_EVIDENCE_2_H

#include <string>
#include <vector>
#include <map>

#include "common_2.h"
#include "settings_2.h"
#include "logger.h"
#include "sequence_file_2.h"

// ----------------------------------------------------

typedef std::map< std::string, std::vector< INTERVAL_EVI > >  EvidenceMap;
typedef std::map< std::string, std::vector< INTERVAL_EVI > >::iterator  EvidenceMapItr;

// ----------------------------------------------------
class Evidence
{
public:

	Evidence( Settings const & settings, Logger * const logger );
	~Evidence(){};

	// key   - seqname from first column of GFF line
	// value - array of evidence intervals; each interval is a parse of GFF line

	EvidenceMap  data;

	void SyncWithSequence( FastaVector const & fasta );
	void HardMaskSequence( FastaVector & fasta );
	
private:

	void SetEvidenceLabels(void);
	void ParseFile( std::string const & name, EvidenceMap & target );
	bool ParseGFFLine( std::string const & source, std::string & seq_name, INTERVAL_EVI & target );
	bool IsBadCodInterval( INTERVAL_EVI const & rec );
	
	void SyncSeqNameWithSizeOne( FastaVector const & fasta, EvidenceMap & evi );
	void VerifySeqNames( FastaVector const & fasta, EvidenceMap const & evi );
	void CheckIntervalsRightBound( FastaVector const & fasta );

	std::map< std::string, EVIDENCE_TYPE > str_to_evidence;

	std::string Status(unsigned int const debug_level);

	// settings

	std::string filename;

	Logger * const logger;

	// drafts
	void SoftMaskToEvidence( FastaVector const & fasta );
	void SoftMaskToIntervals( FastaRecord const & fasta );
};
// ----------------------------------------------------
#endif // GMHMMP2_EVIDENCE_2_H

