// ====================================================
// Author: Alex Lomsadze
// Copyright: Georgia Institute of Technology
// Copyright: GeneProbe Inc.
// File: common_2.h
// Project: GeneMark.hmm-2 (no introns)
// ====================================================

#ifndef GMHMMP2_COMMON_2_H
#define GMHMMP2_COMMON_2_H

// ----------------------------------------------------
namespace TAGS
{
	// general

	static char const AUTO[]    = "auto";
	static char const NO[]      = "no";

	// file format

	static char const FASTA[]   = "fasta";
	static char const PLAIN[]   = "plain";
	static char const GENBANK[] = "genbank";
	static char const EMBL[]    = "embl";
	static char const LST[]     = "lst";
	static char const GTF[]     = "gtf";
	static char const GFF[]     = "gff";
	static char const GFF3[]    = "gff3";
	static char const TRAIN[]   = "train";
	static char const EXT[]     = "ext";

	// strand

	static char const BOTH[]    = "both";
	static char const DIRECT[]  = "direct";
	static char const REVERSE[] = "reverse";
	static char const DOT[]     = ".";
	static char const PLUS[]    = "plus";
	static char const MINUS[]   = "minus";

	// shape

	static char const LINEAR[]       = "linear";
	static char const CIRCULAR[]     = "circular";
	static char const INCOMPLETE[]   = "incomplete";

	// mgm type

	static char const MGM_BACTERIA[] = "bac";
	static char const MGM_ARCHAEA[]  = "arc";

	// evidence 
	// no TABS in lables

	static char const CDS_EVIDENCE[]  = "CDS";
	static char const COD_EVIDENCE[]  = "cod";
	static char const NON_EVIDENCE[]  = "non";
	static char const MASK_EVIDENCE[] = "mask";
	static char const RDNA_EVIDENCE[] = "rDNA";
	static char const TRNA_EVIDENCE[] = "tRNA";
	static char const NRNA_EVIDENCE[] = "ncRNA";

	static char const C_EVIDENCE[] = "c";  // the same as "cod"
	static char const N_EVIDENCE[] = "n";  // the same as "non"
	static char const M_EVIDENCE[] = "m";  // the same as "mask"

	static char const FSH_EVIDENCE[]  = "FSH"; // what to do with this ?

	// output labels

	static char const NATIVE[]   = "native";
	static char const ATYPICAL[] = "atypical";
	static char const BACTERIA[] = "bacteria";
	static char const ARCHAEA[]  = "archaea";
}
// ----------------------------------------------------
namespace FILE_FORMATS
{
	enum FILE_FORMAT
	{
		NON,
		AUTO,
		FASTA,
		PLAIN,
		GENBANK,
		EMBL,
		LST,
		GTF,
		GFF,
		GFF3,
		TRAIN,
		EXT
	};
}
typedef FILE_FORMATS::FILE_FORMAT FILE_FORMAT;

namespace STRAND_TYPES
{
	enum STRAND_TYPE
	{
		NON,
		BOTH,
		DIRECT,
		REVERSE
	};
}
typedef STRAND_TYPES::STRAND_TYPE STRAND_TYPE;

namespace SHAPE_TYPES
{
	enum SHAPE_TYPE
	{
		NON,
		INCOMPLETE,
		LINEAR,
		CIRCULAR
	};
}
typedef SHAPE_TYPES::SHAPE_TYPE SHAPE_TYPE;

namespace MGM_TYPES
{
	enum MGM_TYPE
	{
		NON,
		AUTO,
		BAC,
		ARC
	};
}
typedef MGM_TYPES::MGM_TYPE MGM_TYPE;

namespace EVIDENCE_TYPES
{
	enum EVIDENCE_TYPE
	{
		NON,
		CDS,
		COD,
		MASK,
		RDNA,
		TRNA,
		NRNA,
		FSH
	};
}
typedef EVIDENCE_TYPES::EVIDENCE_TYPE EVIDENCE_TYPE;

// ----------------------------------------------------
namespace TRANSLATION_TABLES
{
	//                            AAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTT
	//                            AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT
	//                            ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
	static char const tbl_11[] = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";
	static char const tbl_4[]  = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF";
	static char const tbl_25[] = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSGCWCLFLF";
	static char const tbl_15[] = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YQYSSSS*CWCLFLF";
}
// ----------------------------------------------------
struct INTERVAL
	{
		int left;
		int right;
		STRAND_TYPE strand;
	};
struct INTERVAL_EVI : INTERVAL
	{
		EVIDENCE_TYPE key;
		double score;
		std::string attr;
	};
// ----------------------------------------------------
// do not change the order: A, C, G, T, N, R, Y
// enum NUCLEOTIDE { A, C, G, T, N, R, Y };

namespace NT
{
	static char const A = 0;
	static char const C = 1;
	static char const G = 2;
	static char const T = 3;
	static char const N = 4;
	static char const R = 5;
	static char const Y = 6;
}

const char LETTER[7]    = { 'A', 'C', 'G', 'T', 'N', 'R', 'Y' };
const char RC_LETTER[7] = { 'T', 'G', 'C', 'A', 'N', 'R', 'Y' };

// ----------------------------------------------------

static double const LOG_ZERO = -1000;
static double const VALUE_ZERO = 0;

// ----------------------------------------------------

static unsigned int const dirStart = 0x1;			// 0000 0000 0000 0000 0000 0000 0000 ---1
static unsigned int const dirEnd   = 0x2;			// 0000 0000 0000 0000 0000 0000 0000 --1-
static unsigned int const revStart = 0x4;			// 0000 0000 0000 0000 0000 0000 0000 -1--
static unsigned int const revEnd   = 0x8;			// 0000 0000 0000 0000 0000 0000 0000 1---

static unsigned int const iniCod = 0x10;			// 0000 0000 0000 0000 0000 0000 ---1 0000
static unsigned int const iniNon = 0x20;			// 0000 0000 0000 0000 0000 0000 --1- 0000
static unsigned int const terCod = 0x40;			// 0000 0000 0000 0000 0000 0000 -1-- 0000
static unsigned int const terNon = 0x80;			// 0000 0000 0000 0000 0000 0000 1--- 0000

													// 0000 0000 0000 0000 0000 ---1 0000 0000
													// 0000 0000 0000 0000 0000 --1- 0000 0000
													// 0000 0000 0000 0000 0000 -1-- 0000 0000
													// 0000 0000 0000 0000 0000 1--- 0000 0000

static unsigned int const frame1 = 0x1000;			// 0000 0000 0000 0000 ---1 0000 0000 0000
static unsigned int const frame2 = 0x2000;			// 0000 0000 0000 0000 --1- 0000 0000 0000
static unsigned int const frame3 = 0x4000;			// 0000 0000 0000 0000 -1-- 0000 0000 0000
static unsigned int const isGap  = 0x8000;			// 0000 0000 0000 0000 1--- 0000 0000 0000

static unsigned int const isATG = 0x10000;			// 0000 0000 0000 ---1 0000 0000 0000 0000
static unsigned int const isGTG = 0x20000;			// 0000 0000 0000 --1- 0000 0000 0000 0000
static unsigned int const isTTG = 0x40000;			// 0000 0000 0000 -1-- 0000 0000 0000 0000
													// 0000 0000 0000 1--- 0000 0000 0000 0000

static unsigned int const isTAA = 0x100000;			// 0000 0000 ---1 0000 0000 0000 0000 0000
static unsigned int const isTAG = 0x200000;			// 0000 0000 --1- 0000 0000 0000 0000 0000
static unsigned int const isTGA = 0x400000;			// 0000 0000 -1-- 0000 0000 0000 0000 0000
													// 0000 0000 1--- 0000 0000 0000 0000 0000

static unsigned int const emptyMark = 0xFF8FFF;	    // 0000 0000 1111 1111 1000 1111 1111 1111

static unsigned int const codMark    = dirEnd   | revStart;
static unsigned int const noncodMark = dirStart | revEnd;

static unsigned int const righOrftMark = dirEnd   | revStart;
static unsigned int const leftOrfMark  = dirStart | revEnd;

static unsigned int const incomplete_ORF  = iniCod | terCod | isGap;

static unsigned int const dirORF = dirEnd | dirStart;
static unsigned int const revORF = revEnd | revStart;

static unsigned int const startMark = dirStart | revStart;

static unsigned int const RESERVED_DURATION = 2000;
static unsigned int const RESERVED_INCOMPLETE_DURATION = 2000;

static unsigned char const NATIVE_TYPE     = 0x1;	// 0000 0001
static unsigned char const ATYPICAL_TYPE_1 = 0x2;	// 0000 0010
static unsigned char const ATYPICAL_TYPE_2 = 0x4;	// 0000 0100

enum GMS2_GROUP {GMS2_A, GMS2_B, GMS2_C, GMS2_D, GMS2_X, GMS2_NONE, GMS2_AC};

// ----------------------------------------------------
class SiteInfo
{
public:
	SiteInfo(){spacer=0;score=0;stype=0;}
	~SiteInfo(){;}

	std::string seq;
	int spacer;
	double score;
	int stype;
};
// ----------------------------------------------------
#endif // GMHMMP2_COMMON_2_H

