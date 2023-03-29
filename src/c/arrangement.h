#define ROW_SIZE 6

typedef struct row {
    int els[ROW_SIZE];
} row;

typedef struct rowgroupinfo {
    row* start;
    int length;
} rowgroupinfo;

typedef struct rowgroup {
    row* rows;
    rowgroupinfo* infos;
    int numrows;
    int numsums;
} rowgroup;