#define VEC_SIZE 6

typedef struct vec {
    int elts[VEC_SIZE];
} vec;

typedef struct vecgroupinfo {
    vec* start;
    int length;
} vecgroupinfo;

typedef struct vecgroup {
    vec* vecs;
    vecgroupinfo* infos;
    int numvecs;
    int numsums;
} vecgroup;

typedef struct relabeling {
    vec* vecs;
    int* label_to_elt;
} relabeling;

typedef struct search_table {
    vec rows[6];
    vec cols[6];
    int numrows;
    int numcols;
    int* valid_rows;
    int num_valid_rows;
    int* valid_cols;
    int num_valid_cols;
} search_table;

typedef struct indexed_int {
    int elt;
    int freq;
} indexed_int;

int sum(vec r);

vecgroup read_vecs(char *filename);

void search_sum(vecgroup group, int sum);

void test_print_vecs(char *filename);