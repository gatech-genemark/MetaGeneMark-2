/*****************************************************
* (c) Copyright 1997, Gene Probe, Inc. All rights reserved.
*
* from original GeneMark distr modified by AL
*
* check.h
*****************************************************/

#pragma warning(disable : 4996)

#ifndef CHECK_H
#define CHECK_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

static long checksum_timekey (const char *tk);
static time_t decode_timekey(const char *s);
static long check_timekey( char* path );
static char* find_path(const char* file, const char* path);

const char LABEL    [] = "GeneMark.hmm prokaryotic 2";
const char KEYF     [] = ".gmhmmp2_key";
const char ENVLABEL [] = "GMHMMP2_KEY"; 

static long checksum_timekey (const char *tk)
{
	int t, c, a, g, i;

	for (i=t=c=a=g=0; i<strlen(tk); ++i)
	switch (tk[i])
	{
		case 't' :
		case 'T' : ++t; break;
		case 'c' :
		case 'C' : ++c; break;
		case 'a' :
		case 'A' : ++a; break;
		case 'g' :
		case 'G' : ++g; break;
	}

	return ((t<<24) | (c<<16) | (a<<8) | g);
}

static time_t decode_timekey(const char *s)
{
	time_t  result = 0;
	int     i, j;

	for (i=2; i<strlen(s)-(sizeof(time_t)*4); ++i)
	if ((s[i] == 'G') && (s[i-1] == 'T') && (s[i-2] == 'A'))
	{
		for (j=0; j<sizeof(time_t)*4; ++j)
		{
			result <<= 2;
			switch (s[i+j+1])
			{
				case 't' :
				case 'T' : break;
				case 'c' :
				case 'C' : result |= 0x1; break;
				case 'a' :
				case 'A' : result |= 0x2; break;
				case 'g' :
				case 'G' : result |= 0x3; break;
			}
		}
		return result;
	}

	return 0;
}

static long check_timekey( char* path )
{
	char*  keyfilename = 0;
	FILE*  keyfile;
	char   keydata[255];
	char   keycrcstring[80];
	long   keycrc;
	time_t keytime, currenttime;
	double time_difference;

	if( path && path[0] != '\0' )
		keyfilename = find_path(KEYF, path);

	if ( !keyfilename )
		keyfilename = find_path(KEYF, getenv(ENVLABEL));
	if ( !keyfilename )
		keyfilename = find_path(KEYF, getenv("HOME"));
	if ( !keyfilename )
		keyfilename = find_path(KEYF, "./");

	if ( !keyfilename )
	{
		fprintf(stderr, "%s\n", LABEL);
		fprintf(stderr, "License key \"%s\" not found.\n", KEYF);
		fprintf(stderr, "This file is neccessary in order to use %s.\n", LABEL);
		exit(1);
	}
  
	if (!(keyfile = fopen(keyfilename, "r")))
	{
		fprintf(stderr, "%s\n", LABEL);
		fprintf(stderr, "Can't read %s\n", keyfilename);
		fprintf(stderr, "Check that the file is readable and not already in use.\n");
		exit(1);
	}

	fgets(keydata, 255, keyfile);
	fgets(keycrcstring, 80, keyfile);
	fclose(keyfile);

	keycrc = atol(keycrcstring);
	if (keycrc != checksum_timekey(keydata))
	{
		fprintf(stderr, "%s\n", LABEL);
		fprintf(stderr, "License key \"%s\" has been altered.\n", keyfilename);
		fprintf(stderr, "Replace the license key with a valid key.\n");
		exit(1);
	}

	keytime     = decode_timekey(keydata);
	currenttime = time(NULL);

	time_difference = difftime(keytime,currenttime);
	if (time_difference > 0)
	{
		time_difference /= (24.0 * 60.0 * 60.0);
		if ( time_difference < 30 ) 
		{
			fprintf(stderr, "%s\n", LABEL);
			fprintf(stderr, "%i %s remaining in the license period.\n\n",
				(int) time_difference+1, (time_difference < 1.0) ? "day" : "days");
		}

		return keycrc;
	}
	else
	{
		fprintf(stderr, "%s\n", LABEL);
		fprintf(stderr, "Your license period has ended. We hope that you found this\n");
		fprintf(stderr, "software useful. If you would like to renew this license,\n");
		fprintf(stderr, "please contact GeneProbe, Inc. at custserv@genepro.com\n");
		fprintf(stderr, "\n");
		exit(1);
	}

	return 0;
}

char* find_path(const char* file, const char* path)
{
	char* tmp;
	char* pathpart;
	char* complete;
	int   flen, plen;
	FILE* fp;

#define PATHSEP ":"
#define PATHTER "/"
#define PATHTERCHAR '\\'

	if (!file) return (char *) 0;
	flen = strlen(file);

	if (!path || (file[0] == PATHTERCHAR))
	{
		if ((fp = fopen(file, "r")))
		{
			fclose(fp);
			complete = (char*) malloc(flen+1);
			strcpy(complete, file);
			return complete;
		}
		return (char *) 0;
	}

	tmp = (char*) malloc(strlen(path)+1);
	strcpy(tmp, path);

	for (pathpart = strtok(tmp, PATHSEP); pathpart; pathpart=strtok(0, PATHSEP))
	{
		plen = strlen(pathpart);
		complete = (char*) malloc(flen + plen + 5);
		strcpy(complete, pathpart);
		if (complete[plen-1] != PATHTERCHAR) strcat(complete, PATHTER);
		strcat(complete, file);
		if ((fp = fopen(complete, "r")))
		{
			fclose(fp);
			free(tmp);
			return complete;
		}
	}

	free(tmp);

	if ((fp = fopen(file, "r")))
	{
		fclose(fp);
		complete =  (char*) malloc(flen+1);
		strcpy(complete, file);
		return complete;
	}

	return (char *) 0;
}

#endif

