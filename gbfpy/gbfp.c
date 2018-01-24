#define __EXTENSIONS__

#include <stdio.h>
#include <stdlib.h>
#include <regex.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <sys/types.h>
 
#include "gbfp.h"

const char sVer[] = "0.5.0";
const char sNorBase[] = "ACGTRYMKWSBDHVNacgtrymkwsbdhvn";
const char sComBase[] = "TGCAYRKMWSVHDBNtgcayrkmwsvhdbn";
const unsigned int iBaseLen = 30;
char sTempLine[LINELEN] = {'\0',};
regex_t ptRegExLocus;
regex_t ptRegExOneLine;
regex_t ptRegExAccession;
regex_t ptRegExVersion;
regex_t ptRegExRegion;
regex_t ptRegExGI;

#define skipSpace( x ) for (; isspace(*x); x++)
#define putLine( x ) strcpy(sTempLine, x)
#define getLine_w_rtrim( x, y ) \
    getLine(x, y); \
    rtrim(x)

static gb_string getLine(gb_string sLine, FILE *FSeqFile) {
    gb_string sReturn;

    if (*sTempLine != '\0') {
        sReturn = strcpy(sLine, sTempLine);
        *sTempLine = '\0';
    } else {
        sReturn = fgets(sLine, LINELEN, FSeqFile);
    }

    return sReturn;
}

static void rtrim(gb_string sLine) {
    register int i;
    
    for (i = (strlen(sLine) - 1); i >= 0; i--) {
        if (! isspace(sLine[i])) {
            sLine[++i] = '\0';
            break;
        }
    }
}

static void removeRChar(gb_string sLine, char cRemove) {
    register int i;
    
    for (i = (strlen(sLine) - 1); i >= 0; i--) {
        if (sLine[i] == cRemove) {
            sLine[i] = '\0';
            break;
        }
    }
}

static gb_string joinLines(FILE *FSeqFile, unsigned int iSpaceLen) {
    char sLine[LINELEN];
    gb_string sTemp, sJoinedLine;

    sJoinedLine = malloc(sizeof(char) * LINELEN);

    getLine_w_rtrim(sLine, FSeqFile);
    strcpy(sJoinedLine, sLine + iSpaceLen);

    while (fgets(sLine, LINELEN, FSeqFile)) {
        sTemp = sLine;
        skipSpace(sTemp);
        if ((sTemp - sLine) < iSpaceLen) break;
        rtrim(sTemp);
        sJoinedLine = strcat(sJoinedLine, sTemp - 1); /* '- 1' in order to insert a space character at the juncation */
    }

    putLine(sLine); 

    return realloc(sJoinedLine, sizeof(char) * (strlen(sJoinedLine) + 1));
}

void initRegEx(void) {
    const char sLocus[] = "^LOCUS +([a-z|A-Z|0-9|_]+) +([0-9]+) bp +([a-z|A-Z|-]+) +([a-z]+) +([A-Z]{3}) (.+)";
    const char sOneLine[] = "^ *([A-Z]+) +(.+)";
    const char sAccession[] = "^ACCESSION +([a-z|A-Z|0-9|_]+) ?";
    const char sRegion[] = " +REGION: ?([0-9]+)\\.\\.([0-9]+)";
    const char sVersion[] = "^VERSION +([a-z|A-Z|0-9|_.]+) ?";
    const char sGI[] = " +GI: ?([0-9]+)";
 
    regcomp(&ptRegExLocus, sLocus, REG_EXTENDED | REG_ICASE);
    regcomp(&ptRegExOneLine, sOneLine, REG_EXTENDED | REG_ICASE);
    regcomp(&ptRegExAccession, sAccession, REG_EXTENDED | REG_ICASE);
    regcomp(&ptRegExVersion, sVersion, REG_EXTENDED | REG_ICASE);
    regcomp(&ptRegExRegion, sRegion, REG_EXTENDED | REG_ICASE);
    regcomp(&ptRegExGI, sGI, REG_EXTENDED | REG_ICASE);
}

static int Pos2Num(gb_string sPositions, int aiPositions[]) {
    register int i;
    int iNum = 0;

    for (i = strlen(sPositions); i >= 0; i--) {
        if (isdigit(*(sPositions + i))) aiPositions[(aiPositions[iNum] - 1 == i) ? iNum : ++iNum] = i;
        else *(sPositions + i) = '\0';
    }

    return iNum;
}

static int Positions2Numbers(gb_string sPositions, unsigned long *lStart, unsigned long *lEnd) {
    int aiPositions[16] = {-2,};
    int iNum;

    iNum = Pos2Num(sPositions, aiPositions);
   
    if (iNum == 2) {
        *lStart = atol(sPositions + aiPositions[2]);
        *lEnd = atol(sPositions + aiPositions[1]);

        return 1;
    } else if (iNum == 1) {
        *lStart = *lEnd = atol(sPositions + aiPositions[1]);
       
        return 1;
    } else {
        fprintf(stderr, "Warning: cannot parse '%s'\n", sPositions);

        return 0;
    }

}

static void parseLocus(gb_string sLocusStr, gb_data *ptGBData) {
    /*    
    01-05      'LOCUS'
    06-12      spaces
    13-28      Locus name
    29-29      space
    30-40      Length of sequence, right-justified
    41-41      space
    42-43      bp
    44-44      space
    45-47      spaces, ss- (single-stranded), ds- (double-stranded), or
               ms- (mixed-stranded)
    48-53      NA, DNA, RNA, tRNA (transfer RNA), rRNA (ribosomal RNA),
               mRNA (messenger RNA), uRNA (small nuclear RNA), snRNA,
               snoRNA. Left justified.
    54-55      space
    56-63      'linear' followed by two spaces, or 'circular'
    64-64      space
    65-67      The division code (see Section 3.3)
    68-68      space
    69-79      Date, in the form dd-MMM-yyyy (e.g., 15-MAR-1991)
    */
    
    char sTemp[LINELEN];
    unsigned int i, iErr, iLen;
    
    regmatch_t ptRegMatch[7];
    
    struct tData {
        char cType;
        void *Pointer;
    } tDatas[] = {
        {STRING, NULL},
        {LONG, NULL},
        {STRING, NULL},
        {STRING, NULL},
        {STRING, NULL},
        {STRING, NULL}};

    tDatas[0].Pointer = ptGBData->sLocusName;
    tDatas[1].Pointer = &(ptGBData->lLength);
    tDatas[2].Pointer = ptGBData->sType;
    tDatas[3].Pointer = ptGBData->sTopology;
    tDatas[4].Pointer = ptGBData->sDivisionCode;
    tDatas[5].Pointer = ptGBData->sDate;
    
    rtrim(sLocusStr);
        
    if ((iErr = regexec(&ptRegExLocus, sLocusStr, 7, ptRegMatch, 0)) == 0) {
        for (i = 0; i < 6; i++) {
            iLen = ptRegMatch[i + 1].rm_eo - ptRegMatch[i + 1].rm_so;
            switch (tDatas[i].cType) {
            case STRING:
                memcpy(tDatas[i].Pointer, (sLocusStr + ptRegMatch[i + 1].rm_so), iLen);
                *((gb_string ) tDatas[i].Pointer + iLen) = '\0';
                break;
            case LONG:
                memcpy(sTemp, (sLocusStr + ptRegMatch[i + 1].rm_so), iLen);
                sTemp[iLen] = '\0';
                *((unsigned long *) tDatas[i].Pointer) = atol(sTemp);
                break;
            default:
                perror("Unknown Data Type!!!");
            }
        }
    } else {
        regerror(iErr, &ptRegExLocus, sTemp, LINELEN);
        fprintf(stderr, "%s\n", sTemp);
        exit(1);
    }
}

static gb_string checkComplement(gb_string sLocation) {
    gb_string sPosition;

    skipSpace(sLocation);

    for (sPosition = sLocation; *sPosition; sPosition++) {
        /* Check the 1st and the 2nd characters of 'complement' */
        if (*sPosition == 'c' && *(sPosition + 1) == 'o') {
            removeRChar(sLocation, ')');
            return sPosition + 11;
        }
    }

    return sLocation;
}

static gb_string checkJoin(gb_string sLocation) {
    gb_string sPosition;

    skipSpace(sLocation);

    for (sPosition = sLocation; *sPosition; sPosition++) {
        /* Check the 1st and the 2nd characters of 'complement' */
        if (*sPosition == 'j' && *(sPosition + 1) == 'o') {
            removeRChar(sLocation, ')');
            return sPosition + 5;
        }
    }

    return sLocation;
}

/* Parsing a gb_string that contains gb_location information */
static void LocationParser(gb_string sLocation, gb_feature *pFeature) {
    gb_string sTemp;
    gb_string sString = NULL;
    
    unsigned int iLocationNum = 1;
    
    /* Evalue sequence direction
    sString has gb_location and join informations
    */

    sString = checkComplement(sLocation);

    if (sLocation == sString) pFeature->cDirection = NORMAL;
    else pFeature->cDirection = REVCOM;

    /* Remove 'join' gb_string
    sString has gb_location informations
    */

    sString = checkJoin(sString);

    sTemp = sString - 1;
    while((sTemp = strchr((sTemp + 1), ','))) iLocationNum++;
    pFeature->ptLocation = malloc(iLocationNum * sizeof(*(pFeature->ptLocation)));

    iLocationNum = 0;
    sLocation = strtok_r(sString, ",", &sTemp);
    if (Positions2Numbers(sLocation,
        &(((pFeature->ptLocation)+iLocationNum)->lStart),
        &(((pFeature->ptLocation)+iLocationNum)->lEnd)) == 1) iLocationNum++;

    while((sLocation = strtok_r(NULL, ",", &sTemp))) {
        if (Positions2Numbers(sLocation,
            &(((pFeature->ptLocation)+iLocationNum)->lStart),
            &(((pFeature->ptLocation)+iLocationNum)->lEnd))) iLocationNum++;
    }
    
    pFeature->lStart = (pFeature->ptLocation)->lStart;
    pFeature->lEnd = ((pFeature->ptLocation)+(iLocationNum - 1))->lEnd;
    pFeature->iLocationNum = iLocationNum;
}

static gb_string parseQualifier(gb_string sQualifier, gb_string *psValue) {
    gb_string sPosition;

    skipSpace(sQualifier);

    if ((sPosition = strchr(sQualifier, '=')) == NULL) {
        *psValue = sQualifier + strlen(sQualifier);
        return sQualifier;
    }

    *sPosition = '\0';

    sPosition++;

    skipSpace(sQualifier);

    if (*sPosition == '"') {
        sPosition++;
        removeRChar(sPosition, '"');
    }

    *psValue = sPosition;

    return sQualifier;
}

static void QualifierParser(gb_string sQualifier, gb_feature *pFeature) {
    gb_string sValue;
    gb_string sTemp = NULL;
    gb_string sString = NULL;
    gb_qualifier *ptQualifier;
    
    pFeature->ptQualifier = malloc(INITQUALIFIERNUM * sizeof(gb_qualifier));
    ptQualifier = pFeature->ptQualifier;

    /* Parse the 1st gb_qualifier gb_string */
    sString = strtok_r(sQualifier, "\n", &sTemp);

    sQualifier = parseQualifier(sString, &sValue);

    ptQualifier->sQualifier = sQualifier;
    ptQualifier->sValue = sValue;
    ptQualifier++;
    
    /* Parse the rest gb_qualifier gb_string */
    while((sString = strtok_r(NULL, "\n", &sTemp)) != NULL) {
        sQualifier = parseQualifier(sString, &sValue);
        ptQualifier->sQualifier = sQualifier;
        ptQualifier->sValue = sValue;
        ptQualifier++;
    }

    pFeature->iQualifierNum = ptQualifier - pFeature->ptQualifier;
    pFeature->ptQualifier =  realloc(pFeature->ptQualifier, pFeature->iQualifierNum * sizeof(gb_qualifier));
}

static unsigned int SequenceParser(gb_string sSequence, gb_string sSequence2) {
    register unsigned int i = 0;
    register unsigned int j = 0;
    register char c;
    
    while((c = *(sSequence + i++)) != '\0')
        if (isalpha(c) != 0) *(sSequence2 + j++) = c;

    *(sSequence2 + j) = '\0';
        
    return j;
}

static void RevCom(gb_string sSequence) {
    char c;
    unsigned int k;
    unsigned long i, j;
    
    for (i = 0, j = strlen(sSequence) - 1; i < j; i++, j--) {
	c = *(sSequence + i);
	*(sSequence + i) = 'X';
	for (k = 0; k < iBaseLen; k++)
	    if (*(sNorBase + k) == *(sSequence + j)) {
		*(sSequence + i) = *(sComBase + k);
		break;
	    }
	*(sSequence + j) = 'X';
	for (k = 0; k < iBaseLen; k++)
	    if (*(sNorBase + k) == c) {
		*(sSequence + j) = *(sComBase + k);
		break;
	}
    }
}

gb_string getSequence(gb_string sSequence, gb_feature *ptFeature) {
    unsigned long lSeqLen = 1; /* For the '\0' characher */
    unsigned long lStart, lEnd;
    unsigned int i;
    gb_string sSequenceTemp;
    
    for (i = 0; i < ptFeature->iLocationNum; i++)
        lSeqLen += (((ptFeature->ptLocation) + i)->lEnd - ((ptFeature->ptLocation) + i)->lStart + 1);
    
    sSequenceTemp = malloc(lSeqLen * sizeof(char));
    
    lSeqLen = 0;
    
    for (i = 0; i < ptFeature->iLocationNum; i++) {
        lStart = ((ptFeature->ptLocation) + i)->lStart;
        lEnd = ((ptFeature->ptLocation) + i)->lEnd;
        strncpy(sSequenceTemp + lSeqLen, sSequence + lStart - 1, lEnd - lStart + 1);
        lSeqLen += (lEnd - lStart + 1);
    }
    
    *(sSequenceTemp + lSeqLen) = '\0';
    
    if (ptFeature->cDirection == REVCOM) RevCom(sSequenceTemp);

    return sSequenceTemp;
}

static void parseFeature(FILE *FSeqFile, gb_data *ptGBData) {
    char sLine[LINELEN] = {'\0',};
    char sLocation[LINELEN] = {'\0',};
    gb_string sQualifier = NULL;
    gb_string sQualifierTemp = NULL;
    unsigned int iReadPos = INELSE;
    unsigned int iFeatureNum = 0;
    unsigned int iFeatureMem = INITFEATURENUM;
    unsigned int i = 0;
    gb_feature *pFeatures = NULL;
    gb_feature *pFeature = NULL;
    
    pFeatures = (gb_feature *) malloc(iFeatureMem * sizeof(gb_feature));

    /* Parse FEATURES */
    while(fgets(sLine, LINELEN, FSeqFile)) {
        if (! isspace(*sLine)) {
            putLine(sLine);
            break;
        }
        
        rtrim(sLine);
        if (memcmp(sLine + 5, "               ", 15) != 0) {
            if (iFeatureNum == iFeatureMem) {
                iFeatureMem += INITFEATURENUM;
                pFeatures = realloc(pFeatures, sizeof(gb_feature) * iFeatureMem);
            }

            if (strlen(sLocation) != 0) LocationParser(sLocation, (pFeatures + iFeatureNum - 1));
            if (sQualifier < sQualifierTemp) {
                *sQualifierTemp++ = '\n';
                *sQualifierTemp = '\0';
                /* printf("=====\n%s=====", sQualifier); */
                sQualifier = realloc(sQualifier, (sQualifierTemp - sQualifier + 1) * sizeof(*sQualifier));
                QualifierParser(sQualifier, (pFeatures + iFeatureNum - 1));
            }

            *sLocation = '\0';
            sQualifier = malloc(sizeof(*sQualifier) * LINELEN);
            sQualifierTemp = sQualifier;
            
            iReadPos = INFEATURE;
            
            memcpy((pFeatures + iFeatureNum)->sFeature, (sLine + 5), 15);
            *(((pFeatures + iFeatureNum)->sFeature) + 15) = '\0';
            rtrim((pFeatures + iFeatureNum)->sFeature);
            strcpy(sLocation, (sLine + 21));

            /* Feature Initalize */
            pFeature = pFeatures + iFeatureNum;
            pFeature->iNum = iFeatureNum;
            pFeature->cDirection = NORMAL;
            pFeature->iLocationNum = 0;
            pFeature->lStart = 0;
            pFeature->lEnd = 0;
            pFeature->iQualifierNum = 0;
            pFeature->ptLocation = NULL;
            pFeature->ptQualifier = NULL;

            iFeatureNum++;
        } else if (*(sLine + QUALIFIERSTART) == '/') {
            iReadPos = INQUALIFIER;
            if (sQualifier < sQualifierTemp) *sQualifierTemp++ = '\n';
            i = strlen(sLine) - (QUALIFIERSTART + 1);
            memcpy(sQualifierTemp, sLine + (QUALIFIERSTART + 1), i);
            sQualifierTemp += i;
        } else {
            if (iReadPos == INFEATURE) {
                strcpy((sLocation + strlen(sLocation)), (sLine + QUALIFIERSTART));
            } else if (iReadPos == INQUALIFIER) {
                i = strlen(sLine) - QUALIFIERSTART;
                memcpy(sQualifierTemp, sLine + QUALIFIERSTART, i);
                sQualifierTemp += i;
            }
        }
    }

    /* Finishing of the parsing */

    if (iFeatureNum == iFeatureMem) {
        iFeatureMem += INITFEATURENUM;
        pFeatures = realloc(pFeatures, sizeof(gb_feature) * iFeatureMem);
    }

    if (strlen(sLocation) != 0) LocationParser(sLocation, (pFeatures + iFeatureNum - 1));
    if (sQualifier < sQualifierTemp) {
        *sQualifierTemp++ = '\n';
        *sQualifierTemp = '\0';
        sQualifier = realloc(sQualifier, (sQualifierTemp - sQualifier + 1) * sizeof(*sQualifier));
        QualifierParser(sQualifier, (pFeatures + iFeatureNum - 1));
    }

    ptGBData->iFeatureNum = iFeatureNum;
    ptGBData->ptFeatures = pFeatures;
}

/* Parse sequences */
static void parseSequence(FILE *FSeqFile, gb_data *ptGBData) {
    char sLine[LINELEN] = {'\0',};
    unsigned long lSeqLen;

    lSeqLen = 0;
    ptGBData->sSequence = malloc((ptGBData->lLength + 1) * sizeof(char));
    
    while(fgets(sLine, LINELEN, FSeqFile)) {
        if (*sLine == '/' && *(sLine + 1) == '/') {
            putLine(sLine);
            break;
        }
        lSeqLen += SequenceParser(sLine, ptGBData->sSequence + lSeqLen);
    }
}

static void parseDef(FILE *FSeqFile, gb_data *ptGBData) {
    char sLine[LINELEN];
    regmatch_t ptRegMatch[3];

    getLine_w_rtrim(sLine, FSeqFile);

    regexec(&ptRegExOneLine, sLine, 3, ptRegMatch, 0);
    ptGBData->sDef = strdup(sLine + ptRegMatch[2].rm_so);
}

static void parseKeywords(FILE *FSeqFile, gb_data *ptGBData) {
    char sLine[LINELEN];
    regmatch_t ptRegMatch[3];

    getLine_w_rtrim(sLine, FSeqFile);

    regexec(&ptRegExOneLine, sLine, 3, ptRegMatch, 0);
    ptGBData->sKeywords = strdup(sLine + ptRegMatch[2].rm_so);
}

static void parseAccession(FILE *FSeqFile, gb_data *ptGBData) {
    char sLine[LINELEN];
    regmatch_t ptRegMatch[3];

    getLine_w_rtrim(sLine, FSeqFile);

    if (regexec(&ptRegExAccession, sLine, 2, ptRegMatch, 0) == 0) {
        *(sLine + ptRegMatch[1].rm_eo) = '\0';
        ptGBData->sAccession = strdup(sLine + ptRegMatch[1].rm_so);
    }

    if (regexec(&ptRegExRegion, sLine + ptRegMatch[1].rm_eo + 1, 3, ptRegMatch, 0) == 0) {
        *(sLine + ptRegMatch[1].rm_eo) = '\0';
        (ptGBData->lRegion)[0] = atol(sLine + ptRegMatch[1].rm_so);
        *(sLine + ptRegMatch[2].rm_eo) = '\0';
        (ptGBData->lRegion)[1] = atol(sLine + ptRegMatch[2].rm_so);
    }
}

static void parseVersion(FILE *FSeqFile, gb_data *ptGBData) {
    char sLine[LINELEN];
    regmatch_t ptRegMatch[2];

    getLine_w_rtrim(sLine, FSeqFile);

    if (regexec(&ptRegExVersion, sLine, 2, ptRegMatch, 0) == 0) {
        *(sLine + ptRegMatch[1].rm_eo) = '\0';
        ptGBData->sVersion = strdup(sLine + ptRegMatch[1].rm_so);
    }
    if (regexec(&ptRegExGI, sLine + ptRegMatch[1].rm_eo + 1, 2, ptRegMatch, 0) == 0) {
        *(sLine + ptRegMatch[1].rm_eo) = '\0';
        ptGBData->sGI = strdup(sLine + ptRegMatch[1].rm_so);
    }
}

static void parseComment(FILE *FSeqFile, gb_data *ptGBData) {
    ptGBData->sComment = joinLines(FSeqFile, 12);
}

static void parseSource(FILE *FSeqFile, gb_data *ptGBData) {
    char sLine[LINELEN];
    regmatch_t ptRegMatch[3];

    getLine_w_rtrim(sLine, FSeqFile);
    regexec(&ptRegExOneLine, sLine, 3, ptRegMatch, 0);
    ptGBData->sSource = strdup(sLine + ptRegMatch[2].rm_so);

    getLine_w_rtrim(sLine, FSeqFile);
    regexec(&ptRegExOneLine, sLine, 3, ptRegMatch, 0);
    ptGBData->sOrganism = strdup(sLine + ptRegMatch[2].rm_so);

    ptGBData->sLineage = joinLines(FSeqFile, 12);
}

#define processRef( x, y ) \
    y = NULL; \
    getLine_w_rtrim(sLine, FSeqFile); \
    putLine(sLine); \
    if (strstr(sLine, x) != NULL) y = joinLines(FSeqFile, 12)

static void parseReference(FILE *FSeqFile, gb_data *ptGBData) {
    char sLine[LINELEN];
    regmatch_t ptRegMatch[3];
    gb_reference *ptReferences = NULL;
    gb_reference *ptReference = NULL;
    unsigned int iReferenceNum = 0;

    ptReferences = ptGBData->ptReferences;
    iReferenceNum = ptGBData->iReferenceNum;

    ptReferences = realloc(ptReferences, sizeof(gb_reference) * (iReferenceNum + 1));
    ptReference = ptReferences + iReferenceNum;

    getLine_w_rtrim(sLine, FSeqFile);
    regexec(&ptRegExOneLine, sLine, 3, ptRegMatch, 0);
    ptReference->iNum = atoi(sLine + ptRegMatch[2].rm_so);

    processRef("AUTHORS", ptReference->sAuthors);
    processRef("TITLE", ptReference->sTitle);
    processRef("CONSTRM", ptReference->sConsrtm);
    processRef("JOURNAL", ptReference->sJournal);
    processRef("PUBMED", ptReference->sPubMed);

    ptGBData->ptReferences = ptReferences;
    ptGBData->iReferenceNum = iReferenceNum + 1;
}

static void initGBData(gb_data *ptGBData) {
    ptGBData->sAccession = NULL;
    ptGBData->sComment = NULL;
    ptGBData->sDef = NULL;
    ptGBData->sGI = NULL;
    ptGBData->sKeywords = NULL;
    ptGBData->sLineage = NULL;
    ptGBData->sOrganism = NULL;
    ptGBData->sSequence = NULL;
    ptGBData->sSource = NULL;
    ptGBData->sVersion = NULL;
    ptGBData->ptReferences = NULL;
    ptGBData->ptFeatures = NULL;
    ptGBData->iFeatureNum = 0;
    ptGBData->iReferenceNum = 0;
    ptGBData->lLength = 0;
    ptGBData->lRegion[0] = 0;
    ptGBData->lRegion[1] = 0;
    ptGBData->sLocusName[0] = '\0';
    ptGBData->sType[0] = '\0';
    ptGBData->sTopology[0] = '\0';
    ptGBData->sDivisionCode[0] = '\0';
    ptGBData->sDate[0] = '\0';
}

static gb_data *_parseGBFF(FILE *FSeqFile) {
    int i;
    char sLine[LINELEN] = {'\0',};
    gb_data *ptGBData = NULL;

    struct tField {
        char sField[FIELDLEN + 1];
        void (*vFunction)(FILE *FSeqFile, gb_data *ptGBData);
    } atFields[] = {
        {"DEFINITION", parseDef},
        {"ACCESSION", parseAccession},
        {"VERSION", parseVersion},
        {"KEYWORDS", parseKeywords},
        {"SOURCE", parseSource},
        {"REFERENCE", parseReference},
        {"COMMENT", parseComment},
        {"FEATURE", parseFeature},
        {"ORIGIN", parseSequence},
        {"", NULL} /* To terminate seeking */
    };

    /* Confirming GBFF File with LOCUS line */
    while(fgets(sLine, LINELEN, FSeqFile)) {
        if (strstr(sLine, "LOCUS") == sLine) {
            ptGBData = malloc(sizeof(gb_data));
            initGBData(ptGBData);
            break;
        }
    }

    /* If there is a no LOCUS line, next statement return NULL value to end parsing */
    if (ptGBData == NULL) return NULL;
   
    /* Parse LOCUS line */ 
    parseLocus(sLine, ptGBData);

    while(getLine(sLine, FSeqFile)) {
        if (strstr(sLine, "//") == sLine) break;
        for(i = 0; *((atFields + i)->sField); i++) {
            if (strstr(sLine, (atFields + i)->sField) == sLine) {
                putLine(sLine);
                ((atFields + i)->vFunction)(FSeqFile, ptGBData);
                break;
            }
        }
    }
    
    return ptGBData;
}

#define freeString( x ) if (x != NULL) free(x)

void freeGBData(gb_data **pptGBData) {
    int i;
    gb_data *ptGBData = NULL;
    gb_feature *ptFeatures = NULL;
    gb_reference *ptReferences = NULL;
    unsigned int iFeatureNum = 0;
    unsigned int iReferenceNum = 0;
    unsigned int iSeqPos = 0;
  
    for (iSeqPos = 0; *(pptGBData + iSeqPos) != NULL; iSeqPos++) {
        ptGBData = *(pptGBData + iSeqPos);

        ptFeatures = ptGBData->ptFeatures;
        iFeatureNum = ptGBData->iFeatureNum;
    
        for (i = 0; i < iFeatureNum; i++) {
            /* printf("%i, %i\n", iFeatureNum, (ptFeatures+iFeatureNum)->iQualifierNum); */
            free((ptFeatures + i)->ptLocation);
            free(((ptFeatures + i)->ptQualifier)->sQualifier);
        }

        free(ptFeatures);


        ptReferences = ptGBData->ptReferences;
        iReferenceNum = ptGBData->iReferenceNum;
        for (i = 0; i < iReferenceNum; i++) {
            freeString((ptReferences + i)->sAuthors);
            freeString((ptReferences + i)->sConsrtm);
            freeString((ptReferences + i)->sTitle);
            freeString((ptReferences + i)->sJournal);
            freeString((ptReferences + i)->sPubMed);
        }

        freeString(ptGBData->sDef);
        freeString(ptGBData->sAccession);
        freeString(ptGBData->sComment);
        freeString(ptGBData->sGI);
        freeString(ptGBData->sKeywords);
        freeString(ptGBData->sLineage);
        freeString(ptGBData->sOrganism);
        freeString(ptGBData->sSequence);
        freeString(ptGBData->sSource);
        freeString(ptGBData->sVersion);
    }
}

gb_data **parseGBFF(gb_string spFileName) {
    int iGBFSeqPos = 0;
    unsigned int iGBFSeqNum = INITGBFSEQNUM;
    gb_data **pptGBDatas;
    FILE *FSeqFile;
 
    if (spFileName == NULL) {
        FSeqFile = stdin;
    } else {
        if (access(spFileName, F_OK) != 0) {
            /* perror(spFileName); */
            return NULL;
        } else {
            FSeqFile = fopen(spFileName, "r");
        }
    }
   
    initRegEx(); /* Initalize for regular expression */

    pptGBDatas = malloc(iGBFSeqNum * sizeof(gb_data *));

    do {
        if (iGBFSeqNum == iGBFSeqPos) {
            iGBFSeqNum += INITGBFSEQNUM;
            pptGBDatas = realloc(pptGBDatas, iGBFSeqNum * sizeof(gb_data *));
        }
        *(pptGBDatas + iGBFSeqPos) = _parseGBFF(FSeqFile);
    } while (*(pptGBDatas + iGBFSeqPos++) != NULL);
    
    if (spFileName) fclose(FSeqFile);

    return pptGBDatas;
}
