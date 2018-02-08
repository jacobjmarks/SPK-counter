#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char * argv[]) {
    if (argc < 3) {
        cerr << "PLEASE SPECIFY FILENAME AND KMER LENGTH" << endl;
        return 1;
    }

    ifstream file;
    file.open(argv[1]);
    if (!file.is_open()) {
        cerr << "FILE CANNOT BE OPENED" << endl;
        return 1;
    }

    file.seekg (0, file.end);
    const int file_length = file.tellg();
    file.seekg (0, file.beg);

    const int KMER_LEN = atoi(argv[2]);

    while(file.peek() != EOF) {
        char c;
        file.get(c);
        if (c == '>') {
            while (c != '\n') file.get(c);
        } else {
            file.putback(c);
        }
        
        char * kmer = new char[KMER_LEN];
        file.read(kmer, KMER_LEN);

        bool valid = true;

        for (int i = 0; i < sizeof(kmer)/sizeof(char); i++) {
            if (kmer[i] == '\n') valid = false;
        }

        if (valid) {
            cout << kmer << endl;

            file.seekg(file.tellg() - KMER_LEN + 1);

            if (file.tellg() + KMER_LEN > file_length) break;
        }
    }

    return 0;
}