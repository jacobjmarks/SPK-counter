#include <iostream>
#include <fstream>
#include <chrono>
#include <unordered_map>
#include <vector>

using namespace std;

unordered_map <string, vector<uint>> counted_kmers;
uint counter_index = 0;

void count_kmer(const string * kmer_p) {
    const string kmer = *kmer_p;
    for (uint i = counted_kmers[kmer].size(); i <= counter_index; i++) {
        counted_kmers[kmer].push_back(0);
    }
    counted_kmers[kmer][counter_index]++;
}

void output_kmer_counts() {
    for (const auto &keyval : counted_kmers) {
        cout << keyval.first;
        for (uint i = 0; i <= counter_index; i++) {
            cout << "\t" << (i < keyval.second.size() ? keyval.second[i] : 0);
        }
        cout << endl;
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
            if (file.tellg() != 1) counter_index++;
            max_buffer_index = 1;

            for (uint i = 0; i < KMER_LEN; i++) {
                kmer_buffer[i].clear();
            }

            file.ignore(UINT32_MAX, '\n');
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

    chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
    chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
    cerr << time_span.count() << 's' << endl;

    cerr <<  "Writing results...";
    output_kmer_counts();
    cerr << "DONE" << endl;

    return 0;
}