// Wyszukiwanie motywów DNA z tzw. "dziurami" (gapami):
//   • '.' – dowolny jeden nukleotyd
//   • '{k}' – dokładnie k dowolnych nukleotydów
//   • 'N' – wildcard dopasowujący A/C/G/T
//
// Ponieważ Aho–Corasick NIE obsługuje dziur ani powtórzeń,
// robimy tzw. podejście SEED-BASED:
//
//   1) Dzielimy wzorzec na tokeny: SEQ i GAP.
//   2) SEQ (ciągłe fragmenty bez dziur) są seedami.
//   3) Dodajemy seedy do AC.
//   4) Gdy AC znajdzie seed, weryfikujemy cały wzorzec
//      cofając się w tekście.
//
// Efekt: szybkie i poprawne dopasowywanie z dziurami.
// Wyniki: liczba dopasowań, czasy, zużycie pamięci.

#include <bits/stdc++.h>
#include <sys/resource.h>
using namespace std;

// Mapowanie liter na indeks 0..4 dla AC
static inline int char_idx(char c){
    switch(c){
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        case 'N': return 4;  // wildcard znak, ale traktujemy jako osobną literę
        default: return 4;
    }
}

// FASTA: łączenie wszystkich sekwencji w jedną
string load_fasta(const string &path){
    ifstream in(path);
    if(!in){
        cerr << "Cannot open FASTA file: " << path << "\n";
        exit(1);
    }
    string line, acc, result;
    while(getline(in, line)){
        if(!line.empty() && line[0] == '>'){
            // nowy rekord FASTA
            if(!acc.empty()){
                result += acc;
                acc.clear();
            }
        } else{
            for(char c: line)
                if(!isspace((unsigned char)c))
                    acc.push_back(toupper(c));
        }
    }
    if(!acc.empty()) result += acc;
    return result;
}

// Wczytywanie wzorców z pliku
vector<string> load_patterns(const string &path){
    ifstream in(path);
    if(!in){
        cerr << "Cannot open patterns: " << path << "\n";
        exit(1);
    }
    string s;
    vector<string> out;
    while(getline(in, s)){
        string clean;
        for(char c: s)
            if(!isspace((unsigned char)c))
                clean.push_back(toupper(c));
        if(!clean.empty())
            out.push_back(clean);
    }
    return out;
}

// Token wzorca: SEQ albo GAP(k)
struct Token{
    bool is_seq;    // czy to SEQ (true) czy GAP (false)
    string seq;     // treść sekwencji (jeśli SEQ)
    int gap;        // długość GAP
};

// Parsowanie wzorca: np. "A..TG{3}C.NN"
vector<Token> parse_pattern(const string &p){
    vector<Token> toks;
    int i = 0, n = p.size();
    while(i < n){
        char c = p[i];

        // 1) SEQ: litery ACGTN
        if(c=='A'||c=='C'||c=='G'||c=='T'||c=='N'){
            int j = i;
            while(j<n && (p[j]=='A'||p[j]=='C'||p[j]=='G'||p[j]=='T'||p[j]=='N'))
                j++;
            toks.push_back({true, p.substr(i,j-i), 0});
            i = j;
        }

        // 2) GAP-y jako krotki '.' → GAP o długości liczby kropek
        else if(c == '.'){
            int j = i;
            while(j<n && p[j]=='.') j++;
            toks.push_back({false, "", j-i});
            i = j;
        }

        // 3) GAP w formie {k}
        else if(c == '{'){
            int j = i+1;
            while(j<n && p[j] != '}') j++;
            if(j<n){
                int k = stoi(p.substr(i+1, j-(i+1)));
                toks.push_back({false, "", k});
                i = j+1;
            } else {
                // błąd składni – traktujemy jako GAP(1)
                toks.push_back({false, "", 1});
                i++;
            }
        }

        // 4) Ignorujemy inne znaki
        else {
            i++;
        }
    }
    return toks;
}

// Seed building: pobieramy tylko SEQ-y długości ≥ min_seed_len
vector<pair<string,int>> build_seeds(const vector<Token> &toks, int min_seed_len){
    vector<pair<string,int>> seeds;
    int offset = 0;  // pozycja w "przetworzonym" wzorcu

    for(const auto &tk: toks){
        if(tk.is_seq){
            // jeśli SEQ jest wystarczająco długi - dodajemy seed
            if((int)tk.seq.size() >= min_seed_len){
                seeds.emplace_back(tk.seq, offset);
            }
            offset += tk.seq.size();
        }
        else {
            offset += tk.gap;
        }
    }
    return seeds;
}

// Aho–Corasick (AC) — wersja uproszczona z 5-literowym alfabetem
struct OutMeta {
    int pat_id;       // ID wzorca
    int seed_offset;  // offset seeda w całym wzorcu
    int seed_len;     // długość seeda
};

struct Aho {
    vector<array<int,5>> next;      // przejścia
    vector<int> fail;                // fail-link
    vector<vector<OutMeta>> out;     // lista dopasowań w danym stanie

    Aho(){
        // stan 0 – root
        next.push_back(array<int,5>{-1,-1,-1,-1,-1});
        fail.push_back(0);
        out.emplace_back();
    }

    // dodawanie słowa (seeda)
    void add_word(const string &s, const OutMeta &meta){
        int v = 0;
        for(char c: s){
            int id = char_idx(c);
            if(next[v][id] == -1){
                next[v][id] = next.size();
                next.push_back(array<int,5>{-1,-1,-1,-1,-1});
                fail.push_back(0);
                out.emplace_back();
            }
            v = next[v][id];
        }
        out[v].push_back(meta);
    }

    // konstrukcja fail-linków BFS-em
    void build_fail(){
        queue<int> q;

        // inicjalizacja poziomu 1
        for(int c=0; c<5; c++){
            int v = next[0][c];
            if(v != -1){
                fail[v] = 0;
                q.push(v);
            } else {
                // trick: brak przejścia z root → przejście do root
                next[0][c] = 0;
            }
        }

        // BFS
        while(!q.empty()){
            int r = q.front(); q.pop();

            for(int c=0; c<5; c++){
                int u = next[r][c];
                if(u == -1) continue;

                q.push(u);

                int v = fail[r];
                while(next[v][c] == -1)
                    v = fail[v];

                fail[u] = next[v][c];

                // dziedziczymy dopasowania ze stanu fail[u]
                for(const auto &m: out[fail[u]])
                    out[u].push_back(m);
            }
        }
    }

    // wyszukiwanie - callback wywoływany przy każdym dopasowaniu seeda
    template<typename F>
    void search_all(const string &text, F &&callback){
        int v = 0;
        for(int i=0; i<(int)text.size(); i++){
            char c = toupper(text[i]);
            if(c!='A' && c!='C' && c!='G' && c!='T' && c!='N'){
                // nieznany znak - reset
                v = 0;
                continue;
            }
            int id = char_idx(c);

            while(next[v][id] == -1)
                v = fail[v];

            v = next[v][id];

            // mamy wyjście?
            for(const auto &m: out[v])
                callback(i, m);
        }
    }
};

// Funkcja oblicza całkowitą długość wzorca (SEQ + GAP)
int total_pattern_length(const vector<Token> &toks){
    int sum = 0;
    for(const auto &tk: toks){
        sum += tk.is_seq ? tk.seq.size() : tk.gap;
    }
    return sum;
}

// Weryfikacja pełnego wzorca po znalezieniu seeda
bool verify_pattern_at(
    const string &text,
    int seed_end,                  // indeks końca seeda w tekście
    int seed_offset,               // offset seeda we wzorcu
    int seed_len,
    const vector<Token> &toks)
{
    // obliczamy start wzorca
    int start = seed_end - (seed_len - 1) - seed_offset;
    if(start < 0) return false;

    int plen = total_pattern_length(toks);
    if(start + plen > (int)text.size()) return false;

    int tpos = start;

    // przechodzimy po tokenach
    for(const auto &tk: toks){
        if(tk.is_seq){
            // dopasowanie litera po literze
            for(int i=0; i<(int)tk.seq.size(); i++){
                char tc = toupper(text[tpos+i]);
                char pc = tk.seq[i];

                if(pc == 'N') continue;       // wildcard

                if(tc != pc)
                    return false;
            }
            tpos += tk.seq.size();
        }
        else {
            // GAP — przesuwamy się o tk.gap znaków
            tpos += tk.gap;
        }
    }

    return true;
}

// Pomiar pamięci — maksimum RSS
long get_rss_kb(){
    rusage ru{};
    getrusage(RUSAGE_SELF, &ru);
#if defined(__APPLE__)
    return ru.ru_maxrss / 1024;
#else
    return ru.ru_maxrss; 
#endif
}

// MAIN
int main(int argc, char **argv){
    if(argc < 3){
        cerr << "Usage: " << argv[0] << " <fasta> <patterns.txt> [min_seed_len]\n";
        return 1;
    }

    string fasta = argv[1];
    string patfile = argv[2];
    int min_seed = (argc >= 4) ? stoi(argv[3]) : 3;

    auto t0 = chrono::high_resolution_clock::now();

    // ładowanie sekwencji DNA
    string text = load_fasta(fasta);
    auto t1 = chrono::high_resolution_clock::now();
    cerr << "Loaded FASTA length: " << text.size() << "\n";

    // wczytywanie wzorców
    auto patterns = load_patterns(patfile);
    cerr << "Loaded patterns: " << patterns.size() << "\n";

    // parsowanie wzorców
    vector<vector<Token>> ptok(patterns.size());
    vector<int> plen(patterns.size());
    for(int i=0; i<(int)patterns.size(); i++){
        ptok[i] = parse_pattern(patterns[i]);
        plen[i] = total_pattern_length(ptok[i]);
    }

    // budowa automatu AC
    Aho ac;

    for(int pid=0; pid<(int)patterns.size(); pid++){
        auto seeds = build_seeds(ptok[pid], min_seed);

        // fallback jeśli wzorzec nie ma seeda
        if(seeds.empty()){
            for(const auto &tk : ptok[pid]){
                if(tk.is_seq && !tk.seq.empty()){
                    ac.add_word(tk.seq, {pid, 0, (int)tk.seq.size()});
                    break;
                }
            }
        }
        else {
            for(auto &s : seeds){
                ac.add_word(s.first, {pid, s.second, (int)s.first.size()});
            }
        }
    }

    ac.build_fail();
    auto t2 = chrono::high_resolution_clock::now();
    cerr << "Automaton nodes: " << ac.fail.size() << "\n";

    // WYSZUKIWANIE
    unordered_map<int, vector<pair<int,int>>> matches;
    size_t total_hits = 0;

    auto t3 = chrono::high_resolution_clock::now();

    ac.search_all(text, [&](int endpos, const OutMeta &m){
        if(verify_pattern_at(text, endpos, m.seed_offset, m.seed_len, ptok[m.pat_id])){
            int start = endpos - (m.seed_len - 1) - m.seed_offset;
            matches[m.pat_id].push_back({start, start + plen[m.pat_id]});
            total_hits++;
        }
    });

    auto t4 = chrono::high_resolution_clock::now();

    double load_t  = chrono::duration<double>(t1 - t0).count();
    double build_t = chrono::duration<double>(t2 - t1).count();
    double search_t = chrono::duration<double>(t4 - t3).count();

    cerr << "Search time: " << search_t << " s\n";
    cerr << "Total matches: " << total_hits << "\n";
    cerr << "RSS: " << get_rss_kb() << " KB\n";

    // podsumowanie
    cout << "Sequence: " << fasta << "\n";
    cout << "Patterns: " << patfile << "\n";
    cout << "FASTA length: " << text.size() << "\n";
    cout << "Patterns count: " << patterns.size() << "\n";
    cout << "Search time: " << search_t << "\n";
    cout << "Total matches: " << total_hits << "\n";

    return 0;
}
