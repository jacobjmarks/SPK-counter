#include <iostream>
#include <fstream>
#include <chrono>
#include <map>

using namespace std;

map <string, uint> counted_kmers;

void count_kmer(const string * kmer_p) {
    const string kmer = *kmer_p;
    map<string, uint>::const_iterator it;
    it = counted_kmers.find(kmer);
    if (it != counted_kmers.end()) {
        counted_kmers[kmer] += 1;
    } else {
        counted_kmers.emplace(kmer, 1);
    }
    counted_kmers.begin();
}

void output_kmer_counts() {
    map<string, uint>::const_iterator it = counted_kmers.begin();
    for (; it != counted_kmers.end(); it++) {
        cout << it->first << " " << it->second << endl;
    }
}

int main(int argc, char* argv[]) {
    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
    if (argc < 3) throw invalid_argument("Please specify filename and kmer length.");

    ifstream file;
    file.open(argv[1]);

    if (!file.is_open()) throw runtime_error("Cannot open file.");

    const uint KMER_LEN = atoi(argv[2]);

    string kmer_buffer[KMER_LEN];
    uint max_buffer_index = 1;

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
                    count_kmer(&kmer_buffer[i]);
                    kmer_buffer[i].clear();
                }
            }

            if (max_buffer_index < KMER_LEN) max_buffer_index++;
        }
    }

    file.close();

    output_kmer_counts();

    chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
    chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);

    cerr << time_span.count() << 's' << endl;

    return 0;
}