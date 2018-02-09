#include <iostream>
#include <fstream>
#include <chrono>
#include <unordered_map>
#include <vector>
#include <iomanip>

using namespace std;

const bool COUNT_CANONICAL = false;

const string ALPHABET = "ACTG";
uint KMER_LEN;

unordered_map <string, uint> counted_kmers;
vector<uint> subseq_indices;

char complement(const char nucleotide) {
    switch (nucleotide) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'T': return 'A';
        case 'G': return 'C';
        default: throw runtime_error((string)"Unhandled nucleotide: " + nucleotide);
    }
}

void get_subseq_indices(const char * filename) {
    ifstream file;
    file.open(filename);

    if (!file.is_open) throw runtime_error((string)"Cannot open file: " + filename);

    string line;
    while (getline(file, line)) {
        if (line[0] == '>') subseq_indices.push_back(file.tellg());
        // if (line[0] == '>') cerr << file.tellg() << endl;
    }
    subseq_indices.push_back(-1);
    file.close();
}

void count_kmer(const string * kmer_p) {
    if (COUNT_CANONICAL) {
        string reverse_complement;
        for (int i = KMER_LEN-1; i >= 0; i--) reverse_complement.push_back(complement((*kmer_p)[i]));
        if (counted_kmers.find(reverse_complement) != counted_kmers.end()) {
            counted_kmers[reverse_complement]++;
            return;
        }
    }
    counted_kmers[*kmer_p]++;
}

void * count_kmers(char * filename, uint read_start, int read_stop) {
    ifstream stream;
    stream.open(filename);

    stream.seekg(read_start);

    string kmer_buffer[KMER_LEN];
    uint max_buffer_index = 1;

    char ch;
    while (stream.tellg() != read_stop) {
        ch = stream.get();

        if (ch == '\n' || ch == '\r') {
            // Skip
        } else if (ALPHABET.find(ch) == string::npos) {
            // Char is not part of nucleotide alphabet, break current kmers
            if (ch == '>') break;
            max_buffer_index = 1;
            for (uint i = 0; i < KMER_LEN; i++) kmer_buffer[i].clear();
        } else {
            for (uint i = 0; i < max_buffer_index; i++) {
                kmer_buffer[i] += ch;

                if (kmer_buffer[i].length() == KMER_LEN) {
                    count_kmer(&kmer_buffer[i]);
                    // cerr << kmer_buffer[i] << endl;
                    kmer_buffer[i].clear();
                }
            }

            if (max_buffer_index < KMER_LEN) max_buffer_index++;
        }
    }

    stream.close();
}

void output_kmer_counts() {
    for (const auto &keyval : counted_kmers) {
        cout << keyval.first << "\t" << keyval.second << endl;
    }
}

int main(int argc, char* argv[]) {
    if (argc < 3) throw invalid_argument("Please specify filename and kmer length.");

    get_subseq_indices(argv[1]);

    KMER_LEN = atoi(argv[2]);

    cerr << "Counting K-mers...";
    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();

    for (uint i = 0; i < subseq_indices.size() - 1; i++) {
        count_kmers(argv[1], subseq_indices[i], subseq_indices[i+1]);
    }

    chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
    chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
    cerr << time_span.count() << 's' << endl;

    cerr <<  "Writing results...";
    output_kmer_counts();
    cerr << "DONE" << endl;

    return 0;
}