#include <iostream>
#include <fstream>

using std::ifstream;
using std::streampos;
using std::ios;


// Reads a data file with name `filename` into memory and returns a pointer to
// the beginning of the data.  The value pointed to by `len` is changed to show
// the size of the storage.  Note that the calling function must now take
// ownership of the data and is responsible for freeing it.
//
// The code in this function is adapted from
// http://www.cplusplus.com/doc/tutorial/files/.  See the link (near the end of
// the web page) for details about the code.

template <class T>
T* read_data(std::string& filename, int* len) {

    streampos size;
    char* memblock;

    ifstream file(filename.c_str(), ios::in|ios::binary|ios::ate);
    if (file.is_open()) {
	size = file.tellg();
	memblock = new char[size];
	file.seekg (0, ios::beg);
	file.read (memblock, size);
	file.close();
    }

    *len = size / sizeof(T);
    return (T*) memblock;
}
