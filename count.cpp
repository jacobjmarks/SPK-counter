#include <iostream>
#include <fstream>
#include <chrono>
#include <unordered_map>
#include <vector>
#include <iomanip>
#include <string.h>

using namespace std;

const string ALPHABET = "ACTG";
uint KMER_LEN;
bool COUNT_CANONICAL = false;
unordered_map <string, vector<uint>> COUNTED_KMERS;
uint NUM_FILES = 0;

char complement(const char nucleotide) {
    switch (nucleotide) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'T': return 'A';
        case 'G': return 'C';
        default: throw runtime_error((string)"Unhandled nucleotide: " + nucleotide);
    }
}

void count_kmer(const string * kmer_p) {
    if (COUNT_CANONICAL) {
        string reverse_complement;
        for (int i = KMER_LEN-1; i >= 0; i--) reverse_complement.push_back(complement((*kmer_p)[i]));
        if (COUNTED_KMERS.find(reverse_complement) != COUNTED_KMERS.end()) {
            vector<uint> * counts = &COUNTED_KMERS[reverse_complement];
            while ((*counts).size() < NUM_FILES) (*counts).push_back(0);
            (*counts)[NUM_FILES-1]++;
            return;
        }
    }

    vector<uint> * counts = &COUNTED_KMERS[*kmer_p];
    while ((*counts).size() < NUM_FILES) (*counts).push_back(0);
    (*counts)[NUM_FILES-1]++;
}

void output_kmer_counts() {
    for (const auto &keyval : COUNTED_KMERS) {
        cout << keyval.first;
        for (uint i = 0; i < NUM_FILES; i ++) {
            cout << "\t" << (i < keyval.second.size() ? keyval.second[i] : 0);
        }
        cout << endl;
    }
}

void process_file(string filename) {
    ifstream file;
    file.open(filename);
    if (!file.is_open()) throw runtime_error((string)"Cannot open file: " + filename);

    NUM_FILES++;
    cout << '\t' << filename;
    
    string kmer_buffer[KMER_LEN];
    uint max_buffer_index = 1;

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
}

int main(int argc, char* argv[]) {
    auto usage = [&argv](string message){
        cerr << message << endl;
        fprintf(stderr, "Usage: %s -k kmer_len [-C] fileA [fileB ...]\n", argv[0]);
        fprintf(stderr, "  -k kmer_len : Kmer length as positive nonzero integer.\n");
        fprintf(stderr, "  -C          : Count canonical kmers.\n");
        return 1;
    };

    if (argc < 3) return usage("Not enough arguments!");

    vector<string> files;

    for (int i = 1; i < argc; i++) {
        auto arg = [&argv, i](string arg){
            return !strcmp(argv[i], arg.c_str());
        };

        if (arg("-k")) {
            KMER_LEN = atoi(argv[++i]);
            if (KMER_LEN < 1) {
                return usage("Bad kmer length!");
            }
            continue;
        }

        if (arg("-C")) {
            COUNT_CANONICAL = true;
            continue;
        }
        
        files.push_back(argv[i]);
    }

    if (!KMER_LEN) return usage("No kmer length specified!");

    for (uint i = 0; i < files.size(); i++) {
        cerr << files[i] << endl;
        cerr << "\tCounting K-mers...";
        chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();

        process_file(files[i]);

        chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
        chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
        cerr << time_span.count() << 's' << endl;
    }

    cout << '\n';

    cerr <<  "Writing results...";
    output_kmer_counts();
    cerr << "DONE" << endl;

    return 0;
}