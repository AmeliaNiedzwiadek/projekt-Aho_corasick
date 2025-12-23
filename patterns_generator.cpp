// Generator wzorców DNA z dziurami zgodnych z systemem:
//   • '.' – pojedyncza dziura
//   • losowe 'N'
//   • pobieranie fragmentów z prawdziwej sekwencji FASTA
//
// Służy do tworzenia dużych zestawów testowych.
// Automatycznie generuje 3 zestawy:
//   10 wzorców długości 10
//   50 wzorców długości 12
//   200 wzorców długości 20

#include <bits/stdc++.h>
using namespace std;

// Wczytywanie FASTA do pojedynczego stringa DNA
string load_fasta(const string &path){
    ifstream in(path);
    if(!in){
        cerr << "Cannot open " << path << "\n";
        exit(1);
    }
    string line, acc, result;
    while(getline(in, line)){
        if(!line.empty() && line[0] == '>'){
            if(!acc.empty()){
                result += acc;
                acc.clear();
            }
        } else {
            for(char c: line)
                if(!isspace((unsigned char)c))
                    acc.push_back(toupper(c));
        }
    }
    if(!acc.empty()) result += acc;
    return result;
}

// Dodawanie dziur '.' do wzorca z szansą ~gap_frac
string add_gaps(string s, double gap_frac, mt19937 &rng){
    if(gap_frac <= 0)
        return s;

    int to_gap = max(1, (int)round(s.size()*gap_frac));
    uniform_int_distribution<int> dist(0, s.size()-1);

    unordered_set<int> used;
    while((int)used.size() < to_gap)
        used.insert(dist(rng));

    for(int idx: used)
        s[idx] = '.';

    return s;
}

// Zapis listy wzorców do pliku
void save_patterns(const vector<string> &pats, const string &path){
    ofstream out(path);
    for(const string &p: pats)
        out << p << "\n";
}

// MAIN
int main(int argc, char **argv){
    if(argc < 3){
        cerr << "Usage: " << argv[0] << " <fasta> <prefix> [gap_fraction]\n";
        return 1;
    }

    string fasta = argv[1];
    string prefix = argv[2];
    double gap_frac = (argc >= 4) ? stod(argv[3]) : 0.2;

    string text = load_fasta(fasta);
    mt19937 rng(123456);  // deterministyczne generowanie

    // 3 konfiguracje: liczba wzorców i długość
    vector<pair<int,int>> configs = {
        {10, 10},
        {50, 12},
        {200, 20}
    };

    for(auto &cfg : configs){
        int count = cfg.first;
        int len   = cfg.second;

        vector<string> pats;
        uniform_int_distribution<int> pos(0, text.size()-len-1);

        for(int i=0; i<count; i++){
            int p = pos(rng);
            string s = text.substr(p, len);

            // dodajemy dziury
            s = add_gaps(s, gap_frac, rng);

            pats.push_back(s);
        }

        string fname = prefix + "_" + to_string(count) + ".txt";
        save_patterns(pats, fname);

        cerr << "[OK] Saved " << fname << "\n";
    }

    return 0;
}
