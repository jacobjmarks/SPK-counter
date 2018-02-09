#include <iostream>
#include <fstream>
#include <chrono>
#include <unordered_map>
#include <vector>
#include <iomanip>

using namespace std;

const bool NORMALISE = true;
const uint PRECISION = 6;

uint KMER_LEN;

unordered_map <string, uint> counted_kmers;

void output_kmer_counts() {
    for (const auto &keyval : counted_kmers) {
        cout << keyval.first << "\t" << keyval.second << endl;
    }
}

int main(int argc, char* argv[]) {
    if (argc < 3) throw invalid_argument("Please specify filename and kmer length.");

    ifstream file;
    file.open(argv[1]);

    if (!file.is_open()) throw runtime_error("Cannot open file.");

    KMER_LEN = atoi(argv[2]);

    string kmer_buffer[KMER_LEN];
    uint max_buffer_index = 1;

    cerr << "Counting K-mers...";
    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();

    char ch;
    while ((ch = file.get()) != EOF) {
        if (ch == '>') {
            max_buffer_index = 1;

            for (uint i = 0; i < KMER_LEN; i++) kmer_buffer[i].clear();

            file.ignore(UINT32_MAX, '\n');
        } else if (ch == '\n' || ch == '\r') {
            // Skip
        } else {
            for (uint i = 0; i < max_buffer_index; i++) {
                kmer_buffer[i] += ch;

                if (kmer_buffer[i].length() == KMER_LEN) {
                    counted_kmers[kmer_buffer[i]]++;
                    kmer_buffer[i].clear();
                }
            }

            if (max_buffer_index < KMER_LEN) max_buffer_index++;
        }
    }

    file.close();

    chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
    chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
    cerr << time_span.count() << 's' << endl;

    cerr <<  "Writing results...";
    output_kmer_counts();
    cerr << "DONE" << endl;

    return 0;
}