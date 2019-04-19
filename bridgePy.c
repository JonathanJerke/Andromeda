//
//  bridgePy.c
//  Feynman
//
//  Created by Jonathan Jerke on 9/30/18.
//  Copyright © 2018 Jonathan Jerke. All rights reserved.
//

#include "bridgePy.h"


//PyObject* YinitCal(PyObject * self, PyObject* ptr ){
//    printf("%u %u\n", self, &YinitCal);
//    //self = (PyObject*)&YinitCal;
//    double * ic = malloc(sizeof(double));
//    *ic= 5;
//    PyObject * po = malloc(sizeof(PyObject));
//
//    po = PyObject(5);
//    return po;
//}


static PyObject *
andromeda_free(PyObject *self, PyObject *args)
{
    struct calculation *c = PyLong_AsVoidPtr(args);
    fModel(c);
    free(c);
    return PyLong_FromLong(0);;
}

static PyObject *
andromeda_init(PyObject *self, PyObject *args){
    const char *command;
    struct calculation *c = malloc(sizeof(struct calculation));
    if (!PyArg_ParseTuple(args, "s", &command)){
        *c = initCal();
    } else {
        FILE * in = fopen(command,"r");
        if ( in == NULL ){
            printf("no file\n");
            return NULL;
        }
        initCalculation(c);
        if ( readInput(c,in ) ){
            fclose(in);
            return NULL;
        }
        finalizeInit(c);
        fclose(in);
    }
    iModel(c);
    t1BodyConstruction ( c, eigen);
    tSAboot(c);
    tCollect(&c->i.c,c->i.irrep,eigenVectors+c->i.nStates,c->i.qFloor,c->i.seekPower);
    return PyLong_FromVoidPtr(c);

}

static PyObject *
    andromeda_ham(PyObject *self, PyObject *args)
{
    PyObject *retr;
    const char *command;
    struct calculation *c = malloc(sizeof(struct calculation));
    if (!PyArg_ParseTuple(args, "s", &command)){
        *c = initCal();
    } else {
        FILE * in = fopen(command,"r");
        if ( in == NULL ){
            printf("no file\n");
            return NULL;
        }
        initCalculation(c);
        if ( readInput(c,in ) ){
            fclose(in);
            return NULL;
        }
        finalizeInit(c);
        fclose(in);
    }
    iModel(c);
    retr = PyLong_FromVoidPtr(streams(&c->i.c,eigenVectors,0,0));
    fModel(c);
    free(c);
    
    return retr;
}


//static PyObject *
//andromeda_vec(PyObject *self, PyObject *args)
//{
//    struct calculation *c = PyLong_AsVoidPtr(args);
//    
//    return PyLong_FromVoidPtr(streams(&c->i.c,eigenVectors,0,0));
//}


//static PyObject *
//    andromeda_test(PyObject *self, PyObject *args)
//{
//
//    PyObject *eigenList;
//    PyObject *vector;
//
//    INT_TYPE i;
//
//
//    struct calculation *c = PyLong_AsVoidPtr(args);
//
//
//    eigenList = PyTuple_New(c->i.nStates);
//
//    for (i = 0; i < c->i.nStates ; i++ ){
//        vector = PyTuple_New(CanonicalRank(&c->i.c, eigenVectors + i , 0));
//        for ( r = 0 ; r < CanonicalRank(&c->i.c, eigenVectors + i , 0) ; r++){
//            for ( d = 0; d < 3 ; dim++){
//
//
//
//
//
//            }
//
//
//        }
//
//
//        PyTuple_SetItem(t, i, vector);
//        Py_DECREF(vector);
//
//    }
//
//
//
//
//
//
//    PyTuple_SetItem(t, 0, PyLong_FromLong(1L));
//    PyTuple_SetItem(t, 1, PyLong_FromLong(2L));
//    PyTuple_SetItem(t, 2, PyUnicode_FromString("three"));
//    return t;
//}


//PyMODINIT PyInit_andromeda(void){
PyObject* PyInit_andromeda(void){
    static PyMethodDef pym[] = {
        {"free",  (PyCFunction)andromeda_free, METH_O,
            "frees calculation"},
        {"init",  (PyCFunction)andromeda_init, METH_VARARGS,
            "boots calculation"},
        {"ham",  (PyCFunction)andromeda_ham, METH_VARARGS,
            "gives green eggs and Hamiltonians"},

      //  {"sys",  (PyCFunction)andromeda_vec, METH_O,
        //    "play"},
//        {"sys",  (PyCFunction)andromeda_test, METH_O,
//            "play"},

        {NULL, NULL, 0, NULL}        /* Sentinel */
    };
    
    
    static struct PyModuleDef py = {
        PyModuleDef_HEAD_INIT,
        "andromeda",     /* m_name */
        "This is a module",  /* m_doc */
        -1,                  /* m_size */
        pym ,               /* m_methods */
        NULL,                /* m_reload */
        NULL,                /* m_traverse */
        NULL,                /* m_clear */
        NULL,                /* m_free */
    };

    
    
    return PyModule_Create(&py);
}
