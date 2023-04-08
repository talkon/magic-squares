#ifndef SEARCH_H
#define SEARCH_H

#include <stdlib.h>

#include "types.h"

#define ROW 0
#define COL 1

void search_aux(global_t g, row_table rows, row_table cols, search_table table,
                size_t minvec);

#endif // !LABELLING_H
