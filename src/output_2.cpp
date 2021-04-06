// ====================================================
// Author: Alex Lomsadze
// Copyright: Georgia Institute of Technology
// Copyright: GeneProbe Inc.
// File: output_2.cpp
// Project: GeneMark.hmm-2 (no introns)
// ====================================================

#include <ctime>

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <iomanip>

using std::string;
using std::vector;
using std::endl;

#include "output_2.h"
#include "common_2.h"
#include "exit.h"

// ----------------------------------------------------
Output::Output( Settings const & settings, Logger * const logger, std::string const version ) : logger(logger),version(version)
{
	st = & settings;

	if (settings.out.gid_start)
		gene_count = settings.out.gid_start - 1;
	else
		gene_count = 0;

	fp_format = settings.out.format_t;
	gid_per_contig = settings.out.gid_per_contig;
	gid_label = settings.out.gid_label;

	OpenFile( out, settings.out.filename );

	OpenFile( out_nt, settings.out.nt_filename );
	OpenFile( out_aa, settings.out.aa_filename );

	if (out_nt.is_open())
		out_nt_ref = &out_nt;
	else
		out_nt_ref = 0;

	if (out_aa.is_open())
		out_aa_ref = &out_aa;
	else
		out_aa_ref = 0;

	if (settings.out.add_nt && settings.out.nt_filename.empty())
		out_nt_ref = &out;

	if (settings.out.add_aa && settings.out.aa_filename.empty())
		out_aa_ref = &out;

	gms2_output = false;
	gms2_training_output = false;
	mgm_output = false;

	if ( !settings.hmm.native_filename.empty() )
	{
		gms2_output = true;
	}
	else if ( !settings.hmm.mgm_filename.empty() )
	{
		mgm_output = true;
	}

	if ( fp_format == FILE_FORMATS::TRAIN )
	{
		gms2_training_output = true;
		gms2_output = false;
		mgm_output = false;

		fp_format = FILE_FORMATS::LST;
	}

	complement = SetComplement();
	index = NT2index();

	if ( logger->verbose ) logger->Print( 0, "# Output streams are open ..." );
}
// ----------------------------------------------------
void Output::OpenFile( std::ofstream & ofst, std::string const & name)
{
	if (name.empty())
		return;

	ofst.open( name.c_str() );

	if ( !ofst.is_open() )
		Exit("error on open file:", name);
}
// ----------------------------------------------------
std::string Output::GetTime(void)
{
	time_t now = time(0);
	struct tm * localtm = localtime(&now);
	
	return asctime(localtm);
}
// ----------------------------------------------------
void Output::Header1_LST(void)
{
	out << "# GeneMark.hmm-2 LST format" << endl;
}
// ----------------------------------------------------
void Output::Header1_GFF3(void)
{
	out << "##gff-version 3" << endl;
}
// ----------------------------------------------------
void Output::Header1_GFF(void)
{
	out << "##gff-version 2" << endl;
}
// ----------------------------------------------------
void Output::Header1_shared(Pset const & pset)
{
	out << "# GeneMark.hmm-2 prokaryotic version: " << version << endl;
	out << "# File with sequence: " << st->in.filename << endl;

	if (!st->hmm.native_filename.empty())
	{
		out << "# File with native parameters: " << st->hmm.native_filename << endl;
		out << "# Native species name and build: " << pset.native.at(0)->model_name << " " << pset.native.at(0)->build << endl;
	}

	if (!st->hmm.mgm_filename.empty())
	{
		out << "# File with MetaGeneMark parameters: " << st->hmm.mgm_filename << endl;
	}

	if (!st->in.evidence.empty())
	{
		out << "# File with evidence data: " << st->in.evidence << endl;
	}

	out << "# translation table: " << pset.genetic_code << endl;

	out << "# output date start: " << GetTime(); // endl is part of GetTime
}
// ----------------------------------------------------
void Output::Header2_LST( FastaVectorItr itr )
{
	out << endl;

	if ( ! itr->data.empty() )
		out << "# sequence-region 1 " << itr->data.size() << endl;
	else
		out << "# sequence-region 0 0"<< endl;

	out << "SequenceID: " <<  itr->name << endl;
}
// ----------------------------------------------------
void Output::Header2_GFF3( FastaVectorItr itr )
{
	out << endl;

	if ( ! itr->data.empty() )
		out << "##sequence-region " << itr->name << " 1 " << itr->data.size() << endl;
	else
		out << "##sequence-region " << itr->name << " 0 0"<< endl;
}
// ----------------------------------------------------
void Output::Footer_LST(void)
{
	out << endl;
	out << "# output date end: " << GetTime(); // endl is part of GetTime
	out << "# Done with file: " << st->in.filename << endl;
}
// ----------------------------------------------------
unsigned int Output::AssignGeneIDs(std::vector< BestValue > & predictions, unsigned int last_id)
{
	for (vector< BestValue >::reverse_iterator itr = predictions.rbegin(); itr != predictions.rend(); ++itr)
	{
		last_id += 1;
		itr->id = last_id;
	}

	return last_id;
}
// ----------------------------------------------------
void Output::PrintGenes_LST( std::vector< BestValue > & predictions )
{
	// buffer for sprintf; print into buffer and output buffer to file

	char * buffer = 0;
	unsigned int buffer_size = 2000;
	buffer = ResizeBuffer(buffer, buffer_size);
	int snprintf_size = 0;

	char const * label = 0;

	char const * a = TAGS::NATIVE;
	char const * b = TAGS::ATYPICAL;
	char const * c = TAGS::ATYPICAL;

	if (mgm_output)
	{
		b = TAGS::BACTERIA;
		c = TAGS::ARCHAEA;
	}

	const char* all_labels[5] = {0,0,0,0,0};

	char const native[]   = "native";
	char const atypical[] = "atypical";
	char const bac[] = "bac";
	char const arc[] = "arc";

	if ( gms2_training_output )
	{
		all_labels[0] = native;
		all_labels[1] = bac;
		all_labels[2] = arc;
		all_labels[3] = native;
		all_labels[4] = atypical;
	}

	string attr;
	string gtype_label;
	string training_gtype_label;

	char strand;

	char prefix_L = ' ';
	char prefix_R = ' ';

	for( vector< BestValue >::reverse_iterator itr = predictions.rbegin(); itr != predictions.rend(); ++itr )
	{
		if      (itr->strand == STRAND_TYPES::DIRECT)   strand = '+';
		else if (itr->strand == STRAND_TYPES::REVERSE)  strand = '-';
		else                                            strand = '.';

		// set the attribute values : gtype

		switch (itr->gtype & 7)
		{
		case NATIVE_TYPE:
			label = a;
			break;
		case ATYPICAL_TYPE_1:
			label = b;
			break;
		case ATYPICAL_TYPE_2:
			label = c;
			break;
		}

		gtype_label.assign(label);

		if (!gms2_training_output)
		{
			attr = gtype_label;
		}
		else
		{
			// this string lists all models which predict this gene

			training_gtype_label.clear();

			char gtype = itr->gtype;
			gtype >>= 4;

			if ( gtype & NATIVE_TYPE )      training_gtype_label.append( all_labels[0] );
			if ( gtype & ATYPICAL_TYPE_1 )  training_gtype_label.append( all_labels[1] );
			if ( gtype & ATYPICAL_TYPE_2 )  training_gtype_label.append( all_labels[2] );

			attr = training_gtype_label;
		}

		// why if?
		if ( itr->si.stype )
		{
			attr += " " + itr->si.seq + " " + IntToString(itr->si.spacer) + " " + IntToString(itr->si.stype);
		}

		if ( gms2_training_output )
			attr += " " + gtype_label;

		if ( itr->origin->IsInforced )
		{
			attr += " ";
			if ( itr->strand == STRAND_TYPES::DIRECT )
			{
				attr += evi_dir_orf->at( itr->R - 1 ).attr;
			}
			else if ( itr->strand == STRAND_TYPES::REVERSE )
			{
				attr += evi_rev_orf->at( itr->L - 1 ).attr;
			}
		}

		( itr->complete_L ) ? prefix_L = ' ' : prefix_L = '<';
		( itr->complete_R ) ? prefix_R = ' ' : prefix_R = '>';

		snprintf_size = snprintf( buffer, buffer_size, " %5d %3c   %c%u   %c%u %7u %s", itr->id, strand, prefix_L ,itr->L, prefix_R, itr->R, itr->R - itr->L + 1, attr.c_str() );
		if (snprintf_size > buffer_size)
		{
			buffer_size = snprintf_size;
			buffer = ResizeBuffer(buffer, buffer_size);
			snprintf_size = snprintf(buffer, buffer_size, " %5d %3c   %c%u   %c%u %7u %s", itr->id, strand, prefix_L, itr->L, prefix_R, itr->R, itr->R - itr->L + 1, attr.c_str());
		}
		if (snprintf_size < 0)
			Exit("error on copy string in Output::PrintGenes_LST");

		out << buffer << endl;
	}
}
// ----------------------------------------------------
void Output::PrintGenes_GFF( std::vector< BestValue > & predictions, FastaVectorItr sitr )
{
	// buffer for sprintf; print into buffer and output buffer to file

	char * buffer = 0;
	unsigned int buffer_size = 2000;
	buffer = ResizeBuffer(buffer, buffer_size);
	int snprintf_size = 0;

	// GFF - nine TAB delimited columns
	// seqid - sitr->name.c_str()
	char source[] = "GeneMark.hmm2";
	// type - gene and CDS
	// start - BestValue.L
	// end - BestValue.R
	// score - "." for gene and total logodd for CDS BestValue.origin->prob
	char strand = '.';
	// phase - 0 by design for all genes, including incomplete
	// attributes - combined from 2 fields; gene id and hmm specific information
	string ID;
	string attr;

	char const * label = 0;

	char const * a = TAGS::NATIVE;
	char const * b = TAGS::ATYPICAL;
	char const * c = TAGS::ATYPICAL;

	if (mgm_output)
	{
		b = TAGS::BACTERIA;
		c = TAGS::ARCHAEA;
	}

	for (vector< BestValue >::reverse_iterator itr = predictions.rbegin(); itr != predictions.rend(); ++itr)
	{
		if      (itr->strand == STRAND_TYPES::DIRECT)   strand = '+';
		else if (itr->strand == STRAND_TYPES::REVERSE)  strand = '-';
		else                                            strand = '.';

		if (gid_per_contig)
		{
			// name is composed from sequence ID and gene ID in contig separated by underscore
			ID = "gene_id " + sitr->name + "_" + IntToString(itr->id) + ";";
		}
		else
		{
			// name is composed from gene ID in file separated by underscore
			ID = "gene_id " + IntToString(itr->id) + ";";
		}

		switch (itr->gtype & 7)
		{
		case NATIVE_TYPE:
			label = a;
			break;
		case ATYPICAL_TYPE_1:
			label = b;
			break;
		case ATYPICAL_TYPE_2:
			label = c;
			break;
		}

		attr = ID + " gene_type " + label + ";";

		if (!itr->complete_L || !itr->complete_R)
		{
			if (itr->complete_L)
				attr += " partial 01;";
			else if (itr->complete_R)
				attr += " partial 10;";
			else
				attr += " partial 11;";
		}

		attr += " gc " + IntToString(itr->origin->gc) + ";";
		attr += " length " + IntToString(itr->R - itr->L + 1) + ";";

		snprintf_size = snprintf(buffer, buffer_size, "%s\t%s\tCDS\t%d\t%d\t%.2f\t%c\t0\t%s", sitr->name.c_str(), source, itr->L, itr->R, itr->origin->prob, strand, attr.c_str());
		if (snprintf_size > buffer_size)
		{
			buffer_size = snprintf_size;
			buffer = ResizeBuffer(buffer, buffer_size);
			snprintf_size = snprintf(buffer, buffer_size, "%s\t%s\tCDS\t%d\t%d\t%.2f\t%c\t0\t%s", sitr->name.c_str(), source, itr->L, itr->R, itr->origin->prob, strand, attr.c_str());
		}
		if (snprintf_size < 0)
			Exit("error on copy string in Output::PrintGenes_GFF");
		out << buffer << endl;
	}
}
// ----------------------------------------------------
void Output::PrintGenes_GTF( std::vector< BestValue > & predictions, FastaVectorItr sitr )
{
	// buffer for sprintf; print into buffer and output buffer to file

	char * buffer = 0;
	unsigned int buffer_size = 2000;
	buffer = ResizeBuffer(buffer, buffer_size);
	int snprintf_size = 0;

	// GTF - nine TAB delimited columns
	// seqid - sitr->name.c_str()
	char source[] = "GeneMark.hmm2";
	// type - gene and CDS
	// start - BestValue.L
	// end - BestValue.R
	// score - "." for gene and total logodd for CDS BestValue.origin->prob
	char strand = '.';
	// phase - 0 by design for all genes, including incomplete
	// attributes - combined from 2 fields; gene id and hmm specific information
	string ID;
	string attr;

	char const * label = 0;

	char const * a = TAGS::NATIVE;
	char const * b = TAGS::ATYPICAL;
	char const * c = TAGS::ATYPICAL;

	if (mgm_output)
	{
		b = TAGS::BACTERIA;
		c = TAGS::ARCHAEA;
	}

	for (vector< BestValue >::reverse_iterator itr = predictions.rbegin(); itr != predictions.rend(); ++itr)
	{
		if (itr->strand == STRAND_TYPES::DIRECT)   strand = '+';
		else if (itr->strand == STRAND_TYPES::REVERSE)  strand = '-';
		else                                            strand = '.';

		if (gid_per_contig)
		{
			// name is composed from sequence ID and gene ID in contig separated by underscore
			ID  = "gene_id \""        + sitr->name + "_" + IntToString(itr->id) + "\";";
			ID += " transcript_id \"" + sitr->name + "_" + IntToString(itr->id) + "\";";
		}
		else
		{
			// name is composed from gene ID in file separated by underscore
			ID  = "gene_id \""        + IntToString(itr->id) + "\";";
			ID += " transcript_id \"" + IntToString(itr->id) + "\";";
		}

		switch (itr->gtype & 7)
		{
		case NATIVE_TYPE:
			label = a;
			break;
		case ATYPICAL_TYPE_1:
			label = b;
			break;
		case ATYPICAL_TYPE_2:
			label = c;
			break;
		}

		attr = ID + " gene_type \"" + label + "\";";

		if (!itr->complete_L || !itr->complete_R)
		{
			if (itr->complete_L)
				attr += " partial \"01\";";
			else if (itr->complete_R)
				attr += " partial \"10\";";
			else
				attr += " partial \"11\";";
		}

		attr += " gc \"" + IntToString(itr->origin->gc) + "\";";
		attr += " length \"" + IntToString(itr->R - itr->L + 1) + "\";";

		snprintf_size = snprintf(buffer, buffer_size, "%s\t%s\tCDS\t%d\t%d\t%.2f\t%c\t0\t%s", sitr->name.c_str(), source, itr->L, itr->R, itr->origin->prob, strand, attr.c_str());
		if (snprintf_size > buffer_size)
		{
			buffer_size = snprintf_size;
			buffer = ResizeBuffer(buffer, buffer_size);
			snprintf_size = snprintf(buffer, buffer_size, "%s\t%s\tCDS\t%d\t%d\t%.2f\t%c\t0\t%s", sitr->name.c_str(), source, itr->L, itr->R, itr->origin->prob, strand, attr.c_str());
		}
		if (snprintf_size < 0)
			Exit("error on copy string in Output::PrintGenes_GTF");
		out << buffer << endl;
	}
}
// ----------------------------------------------------
void Output::PrintGenes_GFF3( std::vector< BestValue > & predictions, FastaVectorItr sitr )
{
	// buffer for sprintf; print into buffer and output buffer to file

	char * buffer = 0;
	unsigned int buffer_size = 2000;
	buffer = ResizeBuffer(buffer, buffer_size);
	int snprintf_size = 0;

	// GFF3 - nine TAB delimited columns
	// seqid - sitr->name.c_str()
	char source[] = "GeneMark.hmm2";
	// type - gene and CDS
	// start - BestValue.L
	// end - BestValue.R
	// score - "." for gene and total logodd for CDS BestValue.origin->prob
	char strand = '.';
	// phase - 0 by design for all genes, including incomplete
	// attributes - combined from 3 fields; two canonical ID and Parent and third hmm specific information
	string ID_gene;
	string ID_cds;
	string Parent_cds;
	string cds_attr;

	char const * label = 0;

	char const * a = TAGS::NATIVE;
	char const * b = TAGS::ATYPICAL;
	char const * c = TAGS::ATYPICAL;

	if (mgm_output)
	{
		b = TAGS::BACTERIA;
		c = TAGS::ARCHAEA;
	}

	for( vector< BestValue >::reverse_iterator itr = predictions.rbegin(); itr != predictions.rend(); ++itr )
	{
		if      (itr->strand == STRAND_TYPES::DIRECT)   strand = '+';
		else if (itr->strand == STRAND_TYPES::REVERSE)  strand = '-';
		else                                            strand = '.';

		if (gid_per_contig)
		{
			// name is composed from sequence ID and gene ID in contig separated by underscore
			ID_gene    = "ID=gene_"      + sitr->name + "_" + IntToString(itr->id) + ";";
			ID_cds     = "ID="           + sitr->name + "_" + IntToString(itr->id) + ";";
			Parent_cds = " Parent=gene_" + sitr->name + "_" + IntToString(itr->id) + ";";
		}
		else
		{
			// name is composed from gene ID in file separated by underscore
			ID_gene    = "ID=gene_"      + IntToString(itr->id) + ";";
			ID_cds     = "ID="           + IntToString(itr->id) + ";";
			Parent_cds = " Parent=gene_" + IntToString(itr->id) + ";";
		}

		switch (itr->gtype & 7)
		{
		case NATIVE_TYPE:
			label = a;
			break;
		case ATYPICAL_TYPE_1:
			label = b;
			break;
		case ATYPICAL_TYPE_2:
			label = c;
			break;
		}

		cds_attr = ID_cds + Parent_cds + " gene_type=" + label + ";";

		if (!itr->complete_L || !itr->complete_R)
		{
			if (itr->complete_L)
				cds_attr += " partial=01;";
			else if (itr->complete_R)
				cds_attr += " partial=10;";
			else
				cds_attr += " partial=11;";
		}

		cds_attr += " gc=" + IntToString(itr->origin->gc) + ";";
		cds_attr += " length=" + IntToString(itr->R - itr->L + 1) + ";";

		snprintf_size = snprintf(buffer, buffer_size, "%s\t%s\tgene\t%d\t%d\t.\t%c\t.\t%s", sitr->name.c_str(), source, itr->L, itr->R, strand, ID_gene.c_str());
		if (snprintf_size > buffer_size)
		{
			buffer_size = snprintf_size;
			buffer = ResizeBuffer(buffer, buffer_size);
			snprintf_size = snprintf(buffer, buffer_size, "%s\t%s\tgene\t%d\t%d\t.\t%c\t.\t%s", sitr->name.c_str(), source, itr->L, itr->R, strand, ID_gene.c_str());
		}
		out << buffer << endl;

		snprintf_size = snprintf(buffer, buffer_size, "%s\t%s\tCDS\t%d\t%d\t%.2f\t%c\t0\t%s", sitr->name.c_str(), source, itr->L, itr->R, itr->origin->prob, strand, cds_attr.c_str());
		if (snprintf_size > buffer_size)
		{
			buffer_size = snprintf_size;
			buffer = ResizeBuffer(buffer, buffer_size);
			snprintf_size = snprintf(buffer, buffer_size, "%s\t%s\tCDS\t%d\t%d\t%.2f\t%c\t0\t%s", sitr->name.c_str(), source, itr->L, itr->R, itr->origin->prob, strand, cds_attr.c_str());
		}
		if (snprintf_size < 0)
			Exit("error on copy string in Output::PrintGenes_GFF3");
		out << buffer << endl;
	}
}
// ----------------------------------------------------
void Output::PrintGenes_EXT(std::vector< BestValue > & predictions, FastaVectorItr sitr)
{
	// buffer for sprintf; print into buffer and output buffer to file

	char * buffer = 0;
	unsigned int buffer_size = 2000;
	buffer = ResizeBuffer(buffer, buffer_size);
	int snprintf_size = 0;

	// GFF - nine TAB delimited columns
	// seqid - sitr->name.c_str()
	char source[] = "GeneMark.hmm2";
	// type - gene and CDS
	// start - BestValue.L
	// end - BestValue.R
	// score - "." for gene and total logodd for CDS BestValue.origin->prob
	char strand = '.';
	// phase - 0 by design for all genes, including incomplete
	// attributes - combined from 2 fields; gene id and hmm specific information
	string ID;
	string attr;

	char const * label = 0;

	char const * a = TAGS::NATIVE;
	char const * b = TAGS::ATYPICAL;
	char const * c = TAGS::ATYPICAL;

	if (mgm_output)
	{
		b = TAGS::BACTERIA;
		c = TAGS::ARCHAEA;
	}

	for (vector< BestValue >::reverse_iterator itr = predictions.rbegin(); itr != predictions.rend(); ++itr)
	{
		if (itr->strand == STRAND_TYPES::DIRECT)        strand = '+';
		else if (itr->strand == STRAND_TYPES::REVERSE)  strand = '-';
		else                                            strand = '.';

		if (gid_per_contig)
		{
			// name is composed from sequence ID and gene ID in contig separated by underscore
			ID = "gene_id " + sitr->name + "_" + IntToString(itr->id) + ";";
		}
		else
		{
			// name is composed from gene ID in file separated by underscore
			ID = "gene_id " + IntToString(itr->id) + ";";
		}

		switch (itr->gtype & 7)
		{
		case NATIVE_TYPE:
			label = a;
			break;
		case ATYPICAL_TYPE_1:
			label = b;
			break;
		case ATYPICAL_TYPE_2:
			label = c;
			break;
		}

		attr = ID + " gene_type " + label + ";";

		if (!itr->complete_L || !itr->complete_R)
		{
			if (itr->complete_L)
				attr += " partial 01;";
			else if (itr->complete_R)
				attr += " partial 10;";
			else
				attr += " partial 11;";
		}
		else
		{
			attr += " partial 00;";
		}

		attr += " gc " + IntToString(itr->origin->gc) + ";";
		attr += " length " + IntToString(itr->R - itr->L + 1) + ";";
		attr += " start_lod_score " + DoubleToString(itr->origin->logodd_start, 2) + ";";
		attr += " start_bys_score " + DoubleToString(itr->origin->bayes, 2) + ";";
		attr += " rbs_score " + DoubleToString(itr->origin->logodd_RBS, 2) + ";";
        attr += " promoter_score " + DoubleToString(itr->origin->logodd_Promoter, 2) + ";";
        attr += " rbs_sc_score " + DoubleToString(itr->origin->logodd_RBS_SC, 2) + ";";
        attr += " promoter_sc_score " + DoubleToString(itr->origin->logodd_Promoter_SC, 2) + ";";
        
        attr += " rbs_motif_score " + DoubleToString(itr->origin->logodd_RBS_motif, 2) + ";";
        attr += " promoter_motif_score " + DoubleToString(itr->origin->logodd_Promoter_motif, 2) + ";";
        attr += " rbs_spacer_score " + DoubleToString(itr->origin->logodd_RBS_spacer, 2) + ";";
        attr += " promoter_spacer_score " + DoubleToString(itr->origin->logodd_Promoter_spacer, 2) + ";";
        attr += " start_stop_score " + DoubleToString(itr->origin->logodd_start_stop, 2) + ";";

        attr += " " + itr->si.seq + " " + IntToString(itr->si.spacer) + " " + IntToString(itr->si.stype);
        
        

		if (itr->origin->status & isATG)       attr += " start_codon ATG;";
		else if (itr->origin->status & isGTG)  attr += " start_codon GTG;";
		else if (itr->origin->status & isTTG)  attr += " start_codon TTG;";
		else                                   attr += " start_codon NNN;";

		MapValue* ptr_to_end = itr->origin;

		if (itr->strand == STRAND_TYPES::DIRECT)
		{
			while (ptr_to_end->next)
				ptr_to_end = ptr_to_end->next;
		}
		else if (itr->strand == STRAND_TYPES::REVERSE)
		{
			while (ptr_to_end->prev)
				ptr_to_end = ptr_to_end->prev;
		}

		if (ptr_to_end->status & isTAA)       attr += " stop_codon TAA;";
		else if (ptr_to_end->status & isTAG)  attr += " stop_codon TAG;";
		else if (ptr_to_end->status & isTGA)  attr += " stop_codon TGA;";
		else                                  attr += " stop_codon NNN;";

		snprintf_size = snprintf(buffer, buffer_size, "%s\t%s\tCDS\t%d\t%d\t%.2f\t%c\t0\t%s", sitr->name.c_str(), source, itr->L, itr->R, itr->origin->prob - st->hmm.delta, strand, attr.c_str());
		if (snprintf_size > buffer_size)
		{
			buffer_size = snprintf_size;
			buffer = ResizeBuffer(buffer, buffer_size);
			snprintf_size = snprintf(buffer, buffer_size, "%s\t%s\tCDS\t%d\t%d\t%.2f\t%c\t0\t%s", sitr->name.c_str(), source, itr->L, itr->R, itr->origin->prob - st->hmm.delta, strand, attr.c_str());
		}
		if (snprintf_size < 0)
			Exit("error on copy string in Output::PrintGenes_GFF");
		out << buffer << endl;
	}
}
// ----------------------------------------------------
std::string Output::DoubleToString(double x, unsigned int i)
{
	std::ostringstream ss;
	ss << std::setprecision(i) << std::fixed << x;
	return ss.str();
}
// ----------------------------------------------------
std::string Output::IntToString(int i)
{
	std::ostringstream ss;
	ss << i;
	return ss.str();
}
// ----------------------------------------------------
char * Output::ResizeBuffer(char * buf, unsigned int size)
{
	delete buf;
	buf = 0;
	buf = new char[size + 1];
	if (!buf)
		Exit("out of memory in Output::ResizeBuffer");
	buf[0] = '\0';

	return buf;
}
// ----------------------------------------------------
void Output::Header1(Pset const & pset)
{
	switch( fp_format )
	{
		case FILE_FORMATS::LST:
			Header1_LST();
			break;
		case FILE_FORMATS::GFF:
			Header1_GFF();
			break;
		case FILE_FORMATS::GTF:
			break;
		case FILE_FORMATS::GFF3:
			Header1_GFF3();
			break;
		case FILE_FORMATS::EXT:
			;
			break;
		default:
			Exit("error, unexpected file format found: ", fp_format);
	}

	Header1_shared(pset);
}
// ----------------------------------------------------
void Output::Header2( FastaVectorItr itr )
{	
	switch( fp_format )
	{
		case FILE_FORMATS::LST:
			Header2_LST(itr);
			break;
		case FILE_FORMATS::GFF:
		case FILE_FORMATS::GTF:
		case FILE_FORMATS::GFF3:
		case FILE_FORMATS::EXT:
			Header2_GFF3(itr);
			break;
		default:
			Exit( "error, unexpected file format found: ", fp_format);
	}
}
// ----------------------------------------------------
void Output::Stat(FastaVectorItr sitr, double logodd, std::vector< BestValue > & predictions )
{
	int gene_count_in_contig = 0;
	int total_gene_length_in_contig = 0;

	for (vector< BestValue >::reverse_iterator itr = predictions.rbegin(); itr != predictions.rend(); ++itr)
	{
		gene_count_in_contig += 1;
		total_gene_length_in_contig += itr->R - itr->L + 1;
	}

	double average_gene_count = 0;
	double average_gene_length = 0;
	
	if (gene_count_in_contig)
	{
		average_gene_count = 1000.0 * gene_count_in_contig / sitr->data.size();
		average_gene_length = 1.0 *total_gene_length_in_contig / gene_count_in_contig;
	}

	out << "# " << sitr->name << "\ttotal_logodd\t" << logodd << "\taverage_length\t" << DoubleToString(average_gene_length, 0) << "\taverage_density\t" << DoubleToString(average_gene_count, 2) << endl;
}
// ----------------------------------------------------
void Output::Footer(void)
{
	switch( fp_format )
	{
		case FILE_FORMATS::LST:
			Footer_LST();
			break;
		case FILE_FORMATS::GFF:
		case FILE_FORMATS::GTF:
		case FILE_FORMATS::GFF3:
		case FILE_FORMATS::EXT:
			break;
		default:
			Exit("error, unexpected file format found: ", fp_format);
	}
}
// ----------------------------------------------------
void  Output::PrintGenes( std::vector< BestValue > & predictions, FastaVectorItr sitr, unsigned int gcode)
{
	if (gid_per_contig)
		gene_count = 0;

	gene_count = AssignGeneIDs(predictions, gene_count);

	switch( fp_format )
	{
		case FILE_FORMATS::LST:
			PrintGenes_LST( predictions );
			break;
		case FILE_FORMATS::GFF:
			PrintGenes_GFF( predictions, sitr );
			break;
		case FILE_FORMATS::GTF:
			PrintGenes_GTF( predictions, sitr );
			break;
		case FILE_FORMATS::GFF3:
			PrintGenes_GFF3( predictions, sitr );
			break;
		case FILE_FORMATS::EXT:
			PrintGenes_EXT(predictions, sitr);
			break;
		default:
			Exit("error, unexpected file format found", fp_format);
	}

	if ( out_aa_ref || out_nt_ref )
		PrintGeneSeqToFile( predictions, sitr, gcode );
}
// ----------------------------------------------------
std::string Output::GetNT( vector< BestValue >::reverse_iterator itr, std::string & data )
{
	string seq = data.substr(itr->L - 1, itr->R - itr->L + 1);

	if (itr->strand == STRAND_TYPES::DIRECT)
	{
		;
	}
	else if (itr->strand == STRAND_TYPES::REVERSE)
	{
		seq = RevCompNT(seq);
	}
	else
		Exit("error, unexpected strand value found", itr->strand);

	return seq;
}
// ----------------------------------------------------
std::string Output::RevCompNT(std::string & s)
{
	std::reverse(s.begin(), s.end());

	for (string::iterator i = s.begin(); i != s.end(); ++i)
	{
		*i = complement[ *i];
	}

	return s;
}
// ----------------------------------------------------
std::vector<unsigned char> Output::SetComplement(void)
{
	vector<unsigned char> v(256,'N');

	v['A'] = 'T';
	v['a'] = 't';
	v['T'] = 'A';
	v['t'] = 'a';
	v['C'] = 'G';
	v['c'] = 'g';
	v['G'] = 'C';
	v['g'] = 'c';
	v['N'] = 'N';
	v['n'] = 'n';
	v['R'] = 'Y';
	v['r'] = 'y';
	v['Y'] = 'R';
	v['y'] = 'r';
	v['U'] = 'T';
	v['u'] = 't';

	return v;
}
// ----------------------------------------------------
std::vector<unsigned char> Output::NT2index(void)
{
	vector<unsigned char> v(256, 4);

	v['A'] = 0;
	v['a'] = 0;
	v['T'] = 3;
	v['t'] = 3;
	v['C'] = 1;
	v['c'] = 1;
	v['G'] = 2;
	v['g'] = 2;

	return v;
}
// ----------------------------------------------------
void Output::PrintGeneSeqToFile(std::vector< BestValue > & predictions, FastaVectorItr sitr, unsigned int gcode)
{
	string seq;

	char const * label = 0;

	char const * a = TAGS::NATIVE;
	char const * b = TAGS::ATYPICAL;
	char const * c = TAGS::ATYPICAL;

	if (mgm_output)
	{
		b = TAGS::BACTERIA;
		c = TAGS::ARCHAEA;
	}

	char strand = '.';
	string defline;

	for (vector< BestValue >::reverse_iterator itr = predictions.rbegin(); itr != predictions.rend(); ++itr)
	{
		defline = ">";

		if (!gid_label.empty())
			defline += gid_label;

		if (gid_per_contig)
			defline += sitr->name + "_";

		defline += IntToString(itr->id) + " " + sitr->name + " " + IntToString(itr->L) + " " + IntToString(itr->R);

		if      (itr->strand == STRAND_TYPES::DIRECT)   strand = '+';
		else if (itr->strand == STRAND_TYPES::REVERSE)  strand = '-';
		else                                            strand = '.';

		defline += " ";
		defline += strand;

		switch (itr->gtype & 7)
		{
		case NATIVE_TYPE:
			label = a;
			break;
		case ATYPICAL_TYPE_1:
			label = b;
			break;
		case ATYPICAL_TYPE_2:
			label = c;
			break;
		}

		defline += " gene_type=";
		defline += label;

		if (!itr->complete_L || !itr->complete_R)
		{
			if (itr->complete_L)
				defline += " partial=01";
			else if (itr->complete_R)
				defline += " partial=10";
			else
				defline += " partial=11";
		}

		seq = GetNT(itr, sitr->data);

		if (out_nt_ref)
		{
			*out_nt_ref << defline << endl;

			for (unsigned int pos = 0; pos < seq.length(); pos += 60)
			{
				*out_nt_ref << seq.substr(pos, 60) << endl;
			}
		}

		if (out_aa_ref)
		{
			*out_aa_ref << defline << endl;

			seq = Translate(seq, gcode);

			for (unsigned int pos = 0; pos < seq.length(); pos += 60)
			{
				*out_aa_ref << seq.substr(pos, 60) << endl;
			}
		}
	}
}
// ----------------------------------------------------
std::string Output::Translate(std::string s, unsigned int gcode)
{
	string prot;

	const char * translation_tbl = 0;

	switch (gcode)
	{
	case 11:
	case 1:
		translation_tbl = TRANSLATION_TABLES::tbl_11;
		break;
	case 4:
		translation_tbl = TRANSLATION_TABLES::tbl_4;
		break;
	case 25:
		translation_tbl = TRANSLATION_TABLES::tbl_25;
		break;
	case 15:
		translation_tbl = TRANSLATION_TABLES::tbl_15;
		break;
	default:
		Exit("error, unsupported genetic code specifies: ", gcode);
	}

	for (unsigned int i = 0; i < s.length(); i += 3)
	{
		unsigned char p1 = index[s[i]];
		unsigned char p2 = index[s[i + 1]];
		unsigned char p3 = index[s[i + 2]];

		if ((p1 < 4) && (p2 < 4) && (p3 < 4))
		{
			unsigned int idx = (p1 << 4) + (p2 << 2) + p3;
			char aa = translation_tbl[idx];
			prot.push_back(aa);
		}
		else
			prot.append("X");
	}

	if ((prot.size() > 0)&&(prot[prot.length() - 1] == '*'))
	{
		prot.erase( prot.size() - 1);
	}

	for (unsigned int i = 0; i < prot.length(); i++)
	{
		if (prot[i] == '*')
			Exit("error, in frame stop codon found in protein translation: ", prot);
	}

	return prot;
}
// ----------------------------------------------------

