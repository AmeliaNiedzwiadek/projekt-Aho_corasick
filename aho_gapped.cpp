/**
 * @file aho_gapped.cpp
 * @author Maria Ławniczak (Nr Indeksu: 268544)
 * @brief Wyszukiwanie motywów DNA z dziurami (gapy) i wildcardami (N)
 * @date 2026-01-25
 */

#include <bits/stdc++.h>
#include <sys/resource.h>
using namespace std;

/**
 * @brief Mapowanie znaków DNA na indeksy 0..4 dla automatu
 * Traktujemy A, C, G, T jako standard, a resztę (w tym N) jako indeks 4
 */
static inline int char_idx(char c){
    switch(c){
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        case 'N': return 4;
        default: return 4;
    }
}

/**
 * @brief Wczytywanie pliku FASTA
 * Łączymy wszystkie rekordy w jeden długi ciąg, pomijając nagłówki i białe znaki
 */
string load_fasta(const string &path){
    ifstream in(path);
    if(!in){
        cerr << "Cannot open FASTA file: " << path << "\n";
        exit(1);
    }
    string line, acc, result;
    while(getline(in, line)){
        if(!line.empty() && line[0] == '>'){
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

/**
 * @brief Wczytywanie wzorców tekstowych z pliku (jeden na linię)
 */
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

/**
 * @brief Token wzorca: albo stała sekwencja (SEQ), albo przerwa (GAP)
 */
struct Token{
    bool is_seq;    // true = SEQ, false = GAP
    string seq;     // treść dla SEQ
    int gap;        // długość dla GAP
};

/**
 * @brief Parsowanie wzorca na tokeny (obsługa '.', '{k}' oraz ACGTN)
 */
vector<Token> parse_pattern(const string &p){
    vector<Token> toks;
    int i = 0, n = p.size();
    while(i < n){
        char c = p[i];
        if(c=='A'||c=='C'||c=='G'||c=='T'||c=='N'){
            int j = i;
            while(j<n && (p[j]=='A'||p[j]=='C'||p[j]=='G'||p[j]=='T'||p[j]=='N'))
                j++;
            toks.push_back({true, p.substr(i,j-i), 0});
            i = j;
        }
        else if(c == '.'){
            int j = i;
            while(j<n && p[j]=='.') j++;
            toks.push_back({false, "", j-i});
            i = j;
        }
        else if(c == '{'){
            int j = i+1;
            while(j<n && p[j] != '}') j++;
            if(j<n){
                int k = stoi(p.substr(i+1, j-(i+1)));
                toks.push_back({false, "", k});
                i = j+1;
            } else {
                toks.push_back({false, "", 1});
                i++;
            }
        }
        else { i++; }
    }
    return toks;
}

/**
 * @brief Budowanie seedów (fragmentów SEQ o minimalnej długości)
 */
vector<pair<string,int>> build_seeds(const vector<Token> &toks, int min_seed_len){
    vector<pair<string,int>> seeds;
    int offset = 0;
    for(const auto &tk: toks){
        if(tk.is_seq){
            if((int)tk.seq.size() >= min_seed_len){
                seeds.emplace_back(tk.seq, offset);
            }
            offset += tk.seq.size();
        }
        else { offset += tk.gap; }
    }
    return seeds;
}

/** @brief Dane wyjściowe automatu dla trafionego seeda */
struct OutMeta {
    int pat_id;       // który to wzorzec
    int seed_offset;  // gdzie we wzorcu jest ten seed
    int seed_len;     // jak długi jest ten seed
};

/** @brief Implementacja automatu AC zoptymalizowana pod alfabet DNA */
struct Aho {
    vector<array<int,5>> next;
    vector<int> fail;
    vector<vector<OutMeta>> out;

    Aho(){
        next.push_back(array<int,5>{-1,-1,-1,-1,-1});
        fail.push_back(0);
        out.emplace_back();
    }

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

    void build_fail(){
        queue<int> q;
        for(int c=0; c<5; c++){
            int v = next[0][c];
            if(v != -1){
                fail[v] = 0;
                q.push(v);
            } else { next[0][c] = 0; }
        }
        while(!q.empty()){
            int r = q.front(); q.pop();
            for(int c=0; c<5; c++){
                int u = next[r][c];
                if(u == -1) continue;
                q.push(u);
                int v = fail[r];
                while(next[v][c] == -1) v = fail[v];
                fail[u] = next[v][c];
                for(const auto &m: out[fail[u]]) out[u].push_back(m);
            }
        }
    }

    template<typename F>
    void search_all(const string &text, F &&callback){
        int v = 0;
        for(int i=0; i<(int)text.size(); i++){
            char c = toupper(text[i]);
            if(c!='A' && c!='C' && c!='G' && c!='T' && c!='N'){
                v = 0; continue;
            }
            int id = char_idx(c);
            while(next[v][id] == -1) v = fail[v];
            v = next[v][id];
            for(const auto &m: out[v]) callback(i, m);
        }
    }
};

/** @brief Sumowanie długości wszystkich tokenów we wzorcu */
int total_pattern_length(const vector<Token> &toks){
    int sum = 0;
    for(const auto &tk: toks) sum += tk.is_seq ? tk.seq.size() : tk.gap;
    return sum;
}

/**
 * @brief Weryfikacja naiwna całego wzorca w tekście po trafieniu seeda
 */
bool verify_pattern_at(const string &text, int seed_end, int seed_offset, int seed_len, const vector<Token> &toks)
{
    int start = seed_end - (seed_len - 1) - seed_offset;
    if(start < 0) return false;
    int plen = total_pattern_length(toks);
    if(start + plen > (int)text.size()) return false;

    int tpos = start;
    for(const auto &tk: toks){
        if(tk.is_seq){
            for(int i=0; i<(int)tk.seq.size(); i++){
                char tc = toupper(text[tpos+i]), pc = tk.seq[i];
                if(pc == 'N') continue;
                if(tc != pc) return false;
            }
            tpos += tk.seq.size();
        } else { tpos += tk.gap; }
    }
    return true;
}

/** @brief Funkcja pomocnicza do pomiaru zużycia pamięci */
long get_rss_kb(){
    rusage ru{};
    getrusage(RUSAGE_SELF, &ru);
#if defined(__APPLE__)
    return ru.ru_maxrss / 1024;
#else
    return ru.ru_maxrss; 
#endif
}

int main(int argc, char **argv){
    if(argc < 3){
        cerr << "Usage: " << argv[0] << " <fasta> <patterns.txt> [min_seed_len]\n";
        return 1;
    }

    // Inicjalizacja i ładowanie danych
    string fasta = argv[1], patfile = argv[2];
    int min_seed = (argc >= 4) ? stoi(argv[3]) : 3;
    auto t0 = chrono::high_resolution_clock::now();

    string text = load_fasta(fasta);
    auto t1 = chrono::high_resolution_clock::now();
    auto patterns = load_patterns(patfile);

    // Przygotowanie wzorców i automatu
    vector<vector<Token>> ptok(patterns.size());
    vector<int> plen(patterns.size());
    for(int i=0; i<(int)patterns.size(); i++){
        ptok[i] = parse_pattern(patterns[i]);
        plen[i] = total_pattern_length(ptok[i]);
    }

    Aho ac;
    for(int pid=0; pid<(int)patterns.size(); pid++){
        auto seeds = build_seeds(ptok[pid], min_seed);
        if(seeds.empty()){
            // jeśli wzorzec jest za krótki na min_seed, weź cokolwiek
            for(const auto &tk : ptok[pid])
                if(tk.is_seq && !tk.seq.empty()){
                    ac.add_word(tk.seq, {pid, 0, (int)tk.seq.size()});
                    break;
                }
        } else {
            for(auto &s : seeds) ac.add_word(s.first, {pid, s.second, (int)s.first.size()});
        }
    }

    ac.build_fail();
    auto t2 = chrono::high_resolution_clock::now();

    // Główne wyszukiwanie
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
    double search_t = chrono::duration<double>(t4 - t3).count();

    // Wyświetlanie wyników
    cout << "FASTA length: " << text.size() << "\n"
         << "Patterns count: " << patterns.size() << "\n"
         << "Search time: " << search_t << " s\n"
         << "Total matches: " << total_hits << "\n"
         << "RSS: " << get_rss_kb() << " KB\n";

    return 0;
}
