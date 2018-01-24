#include "Python.h"
#include "gbfp.h"

static PyObject *ErrorObject;

#define checkString( x ) x != NULL ? x : ""

PyObject* makeLocationList(gb_location *ptLocation, unsigned int iLocationNum) {
    unsigned int i;
    gb_location *ptLocData;
    PyObject *LocationList;
    PyObject *LocationTuple;

    LocationList = PyList_New(0);

    for (i = 0; i < iLocationNum; i++) {
        ptLocData = ptLocation + i;
        
        LocationTuple = PyTuple_New(2);
        PyTuple_SetItem(LocationTuple, 0, PyInt_FromLong((long) ptLocData->lStart));
        PyTuple_SetItem(LocationTuple, 1, PyInt_FromLong((long) ptLocData->lEnd));
        PyList_Append(LocationList, LocationTuple);
    }

    return LocationList;
}

PyObject* makeQualifierList(gb_qualifier *ptQualifier, unsigned int iQualifierNum) {
    unsigned int i;
    gb_qualifier *ptQualData;
    PyObject *QualifierList;
    PyObject *QualifierTuple;

    QualifierList = PyList_New(0);

    for (i = 0; i < iQualifierNum; i++) {
        ptQualData = ptQualifier + i;

        QualifierTuple = PyTuple_New(2);
        PyTuple_SetItem(QualifierTuple, 0, PyString_FromString(ptQualData->sQualifier));
        PyTuple_SetItem(QualifierTuple, 1, PyString_FromString(ptQualData->sValue));
        PyList_Append(QualifierList, QualifierTuple);
    }

    return QualifierList;
}

PyObject* makeFeatureDict(gb_feature *ptFeature) {
    PyObject *FeatureDict;

    FeatureDict =  PyDict_New();

    PyDict_SetItemString(FeatureDict, "feature", PyString_FromString((char *) &(ptFeature->sFeature)));
    PyDict_SetItemString(FeatureDict, "direction", PyString_FromStringAndSize((char *) &(ptFeature->cDirection), 1));
    PyDict_SetItemString(FeatureDict, "start", PyInt_FromLong((long) ptFeature->lStart));
    PyDict_SetItemString(FeatureDict, "end", PyInt_FromLong((long) ptFeature->lEnd));
    PyDict_SetItemString(FeatureDict, "number", PyInt_FromLong((long) ptFeature->iNum));
    PyDict_SetItemString(FeatureDict, "location_num", PyInt_FromLong((long) ptFeature->iLocationNum));
    PyDict_SetItemString(FeatureDict, "qualifier_num", PyInt_FromLong((long) ptFeature->iQualifierNum));
    PyDict_SetItemString(FeatureDict, "location", makeLocationList(ptFeature->ptLocation, ptFeature->iLocationNum));
    PyDict_SetItemString(FeatureDict, "qualifier", makeQualifierList(ptFeature->ptQualifier, ptFeature->iQualifierNum));

    return FeatureDict; 
}

PyObject* makeReferenceDict(gb_reference *ptReference) {
    PyObject *ReferenceDict;

    ReferenceDict =  PyDict_New();

    PyDict_SetItemString(ReferenceDict, "authors", PyString_FromString(checkString(ptReference->sAuthors)));
    PyDict_SetItemString(ReferenceDict, "consrtm", PyString_FromString(checkString(ptReference->sConsrtm)));
    PyDict_SetItemString(ReferenceDict, "title", PyString_FromString(checkString(ptReference->sTitle)));
    PyDict_SetItemString(ReferenceDict, "journal", PyString_FromString(checkString(ptReference->sJournal)));
    PyDict_SetItemString(ReferenceDict, "pubmed", PyString_FromString(checkString(ptReference->sPubMed)));
    PyDict_SetItemString(ReferenceDict, "number", PyInt_FromLong((long) ptReference->iNum));

    return ReferenceDict; 
}

PyObject* makeGBFFDataDict(gb_data *ptGBFFData) {
    int i;

    PyObject *GBFFDataDict;
    PyObject *FeatureList;
    PyObject *ReferenceList;
    PyObject *RegionTuple;

    GBFFDataDict =  PyDict_New();

    PyDict_SetItemString(GBFFDataDict, "locus_name", PyString_FromString((char *) &(ptGBFFData->sLocusName)));
    PyDict_SetItemString(GBFFDataDict, "length", PyInt_FromLong((long) ptGBFFData->lLength));
    PyDict_SetItemString(GBFFDataDict, "type", PyString_FromString((char *) &(ptGBFFData->sType)));
    PyDict_SetItemString(GBFFDataDict, "topology", PyString_FromString((char *) &(ptGBFFData->sTopology)));
    PyDict_SetItemString(GBFFDataDict, "division_code", PyString_FromString((char *) &(ptGBFFData->sDivisionCode)));
    PyDict_SetItemString(GBFFDataDict, "date", PyString_FromString((char *) &(ptGBFFData->sDate)));
    PyDict_SetItemString(GBFFDataDict, "feature_num", PyInt_FromLong((long) ptGBFFData->iFeatureNum));
    PyDict_SetItemString(GBFFDataDict, "sequence", PyString_FromString(checkString(ptGBFFData->sSequence)));
    PyDict_SetItemString(GBFFDataDict, "accession", PyString_FromString(checkString(ptGBFFData->sAccession)));
    PyDict_SetItemString(GBFFDataDict, "comment", PyString_FromString(checkString(ptGBFFData->sComment)));
    PyDict_SetItemString(GBFFDataDict, "definition", PyString_FromString(checkString(ptGBFFData->sDef)));
    PyDict_SetItemString(GBFFDataDict, "gi", PyString_FromString(checkString(ptGBFFData->sGI)));
    PyDict_SetItemString(GBFFDataDict, "keywords", PyString_FromString(checkString(ptGBFFData->sKeywords)));
    PyDict_SetItemString(GBFFDataDict, "lineage", PyString_FromString(checkString(ptGBFFData->sLineage)));
    PyDict_SetItemString(GBFFDataDict, "organism", PyString_FromString(checkString(ptGBFFData->sOrganism)));
    PyDict_SetItemString(GBFFDataDict, "source", PyString_FromString(checkString(ptGBFFData->sSource)));
    PyDict_SetItemString(GBFFDataDict, "version", PyString_FromString(checkString(ptGBFFData->sVersion)));

    /* Create list of region */
    RegionTuple = PyTuple_New(2);

    PyTuple_SetItem(RegionTuple, 0, PyInt_FromLong((long) ptGBFFData->lRegion[0]));
    PyTuple_SetItem(RegionTuple, 1, PyInt_FromLong((long) ptGBFFData->lRegion[1]));

    PyDict_SetItemString(GBFFDataDict, "region", RegionTuple);


    /* Create list of references */
    ReferenceList = PyList_New(0);

    for(i = 0; i < ptGBFFData->iReferenceNum; i++)
        PyList_Append(ReferenceList, makeReferenceDict((ptGBFFData->ptReferences) + i));

    PyDict_SetItemString(GBFFDataDict, "references", ReferenceList);

    /* Create list of features */
    FeatureList = PyList_New(0);

    for(i = 0; i < ptGBFFData->iFeatureNum; i++)
        PyList_Append(FeatureList, makeFeatureDict((ptGBFFData->ptFeatures) + i));

    PyDict_SetItemString(GBFFDataDict, "features", FeatureList);


    return GBFFDataDict;
}

static PyObject* parse(PyObject *self, PyObject *args) {
    int i;
    char *psFileName;

    gb_data **pptGBFFData;

    PyObject *GBFFDataList;

    /* Parsing arguments */
    if (!PyArg_ParseTuple(args, "s", &psFileName)) return NULL;
    
    /* Parsing with C function */
    pptGBFFData = parseGBFF(psFileName);

    if (pptGBFFData != NULL) {
        /* Convert datas from C to Python */
        GBFFDataList = PyList_New(0);

        for (i = 0; *(pptGBFFData + i) != NULL; i++) {
            PyList_Append(GBFFDataList, makeGBFFDataDict(*(pptGBFFData + i)));
        }

        freeGBData(pptGBFFData);
        
        return GBFFDataList;
    } else {
        PyErr_SetString(PyExc_IOError, "File not found !!!");

        return NULL;
    }
}         
 
static struct PyMethodDef gbfpy_methods[] = { 
    {"parse", parse, METH_VARARGS},
    {NULL, NULL}
};

void initgbfpy(void) {
    PyObject *parser; 
    parser = Py_InitModule("gbfpy", gbfpy_methods);
    ErrorObject = Py_BuildValue("s", "GBFF parser error !!!");
} 
