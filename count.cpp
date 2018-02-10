#include <iostream>
#include <fstream>
#include <chrono>
#include <unordered_map>
#include <vector>
#include <iomanip>
#include <deque>

using namespace std;

const bool COUNT_CANONICAL = true;

const string ALPHABET = "ACTG";
uint KMER_LEN;

unordered_map <string, uint> counted_kmers;

char complement(const char nucleotide) {
    switch (nucleotide) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'T': return 'A';
        case 'G': return 'C';
        default: return 0;
    }
}

void count_kmer(const deque<char> kmer_queue) {
    if (COUNT_CANONICAL) {
        string reverse_complement;
        for (int i = KMER_LEN-1; i >= 0; i--) {
            char c = complement(kmer_queue[i]);
            if (c == 0) return;
            reverse_complement += c;
        }
        if (counted_kmers.find(reverse_complement) != counted_kmers.end()) {
            counted_kmers[reverse_complement]++;
            return;
        }
    }
    
    string kmer;
    for (uint i = 0; i < KMER_LEN; i++) kmer += kmer_queue[i];
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

    cerr << "Counting K-mers...";
    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();

    deque<char> kmer_queue;

    char ch;
    while ((ch = file.get()) != EOF) {
        if (ch == '\n' || ch == '\r') {
            // Skip
        } else if (ALPHABET.find(ch) == string::npos) {
            // Char is not part of nucleotide alphabet, break current kmers
            if (ch == '>') file.ignore(UINT32_MAX, '\n');
            kmer_queue.erase(kmer_queue.begin(), kmer_queue.end());

            char * first_kmer = new char [KMER_LEN];
            file.read(first_kmer, KMER_LEN);

            for (uint i = 0; i < KMER_LEN; i++) kmer_queue.push_back(first_kmer[i]);

            delete first_kmer;
            count_kmer(kmer_queue);
        } else {
            kmer_queue.pop_front();
            kmer_queue.push_back(ch);
            count_kmer(kmer_queue);
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