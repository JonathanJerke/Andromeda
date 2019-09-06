//
//  spammify.c
//  Berry
//
//  Created by Jonathan Jerke on 12/12/18.
//  Copyright © 2018 Jonathan Jerke. All rights reserved.
//

#include "spammify.h"

static PyObject *
spam_system(PyObject *self, PyObject *args)
{
    const char *command;
    int sts;
    
    if (!PyArg_ParseTuple(args, "s", &command))
        return NULL;
    sts = system(command);
    return PyLong_FromLong(sts);
}
