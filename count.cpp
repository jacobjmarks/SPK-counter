#include <iostream>
#include <fstream>
#include <chrono>
#include <unordered_map>
#include <vector>
#include <iomanip>

using namespace std;

const bool COUNT_CANONICAL = true;

const string ALPHABET = "ACTG";
uint KMER_LEN;

unordered_map <string, uint> counted_kmers;

char complement(const char * nucleotide) {
    switch (*nucleotide) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'T': return 'A';
        case 'G': return 'C';
        default: throw runtime_error((string)"Unhandled nucleotide: " + *nucleotide);
    }
}

void count_kmer(const string * kmer_p) {
    const string kmer = *kmer_p;

    if (COUNT_CANONICAL) {
        string reverse_complement;
        for (int i = KMER_LEN-1; i >= 0; i--) reverse_complement.push_back(complement(&kmer[i]));
        if (counted_kmers.find(reverse_complement) != counted_kmers.end()) {
            counted_kmers[reverse_complement]++;
            return;
        }
    }
    counted_kmers[kmer]++;
}

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
        if (ch == '\n' || ch == '\r') {
            // Skip
        } else if (ALPHABET.find(ch) == string::npos) {
            // Char is not part of nucleotide alphabet, break current kmers
            if (ch == '>') file.ignore(UINT32_MAX, '\n');
            max_buffer_index = 1;
            for (uint i = 0; i < KMER_LEN; i++) kmer_buffer[i].clear();
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