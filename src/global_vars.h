// used to map from day-specific indices to individuals (d2s is short for
// day-to-subject).  More specifically, if `d2s[k] == i` is true, then the k-th
// index for day-specific data maps to the i-th index for subject-specific data.
extern int* d2s;
