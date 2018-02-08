#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

int main(int argc, char* argv[]) {
    if (argc < 3) throw invalid_argument("Please specify filename and kmer length.");

    ifstream file;
    file.open(argv[1]);

    if (!file.is_open()) throw runtime_error("Cannot open file.");

    const uint KMER_LEN = atoi(argv[2]);

    string kmer_buffer[KMER_LEN];
    int max_buffer_index = 1;

    char ch;
    while ((ch = file.get()) != EOF) {
        if (ch == '>') {
            file.ignore(UINT32_MAX, '\n');

            for (uint i = 0; i < KMER_LEN; i++) {
                kmer_buffer[i].clear();
            }

            max_buffer_index = 1;
        } else if (ch == '\n' || ch == '\r') {
            // Skip
        } else {
            for (uint i = 0; i < max_buffer_index; i++) {
                kmer_buffer[i] += ch;

                if (kmer_buffer[i].length() == KMER_LEN) {
                    cout << kmer_buffer[i] << endl;
                    kmer_buffer[i].clear();
                }
            }

            if (max_buffer_index < KMER_LEN) max_buffer_index++;
        }
    }

    file.close();

    return 0;
}