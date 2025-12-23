// Klasyczna implementacja Aho–Corasick do dopasowywania wielu wzorców jednocześnie.
//
// Ten plik nie obsługuje dziur, wildcardów, {k} itp., jest to plik do:
//  - porównania wydajności z wersją gapped
//  - wizualizacji automatu (eksport do .dot)
//  - demonstrowania działania tries + fail linków

#include <bits/stdc++.h>
using namespace std;

// Alfabet DNA + 'N' (traktowany jak osobna litera)
int char_idx(char c){
    switch(toupper(c)){
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        case 'N': return 4;
        default: return 4; // wszystko inne traktujemy jako N
    }
}

// Wczytywanie wzorców z pliku
vector<string> load_patterns(const string& path){
    ifstream in(path);
    if(!in){
        cerr << "Cannot open " << path << "\n";
        exit(1);
    }
    string line;
    vector<string> out;
    while(getline(in, line)){
        string clean;
        for(char c: line)
            if(!isspace((unsigned char)c))
                clean.push_back(toupper(c));
        if(!clean.empty()) out.push_back(clean);
    }
    return out;
}

// Struktura jednego węzła automatu
struct Node{
    array<int,5> next;  // przejścia
    int fail;           // link fail
    vector<int> out;    // ID wzorców kończących się w tym stanie

    Node(){
        next.fill(-1);
        fail = 0;
    }
};

int main(int argc, char** argv){
    if(argc < 2){
        cerr << "Usage: " << argv[0] << " <patterns.txt> [export.dot]\n";
        return 1;
    }
	
    // Wczytujemy wzorce
    vector<string> pats = load_patterns(argv[1]);

    // TWORZENIE TRIE
    vector<Node> trie(1); // stan 0 = root

    for(size_t pid=0; pid<pats.size(); pid++){
        int v = 0;
        for(char c: pats[pid]){
            int id = char_idx(c);
            if(trie[v].next[id] == -1){
                trie[v].next[id] = trie.size();
                trie.emplace_back();
            }
            v = trie[v].next[id];
        }
        trie[v].out.push_back(pid);
    }

    // BUDOWANIE FAIL-LINKÓW BFS-em
    queue<int> q;

    // inicjalizacja poziomu 1
    for(int c=0; c<5; c++){
        int nxt = trie[0].next[c];
        if(nxt != -1){
            trie[nxt].fail = 0;
            q.push(nxt);
        } else {
            trie[0].next[c] = 0; // brak przejścia - wracamy do root
        }
    }

    // BFS
    while(!q.empty()){
        int r = q.front(); q.pop();
        for(int c=0; c<5; c++){
            int u = trie[r].next[c];
            if(u == -1) continue;

            q.push(u);

            int f = trie[r].fail;
            while(trie[f].next[c] == -1)
                f = trie[f].fail;

            trie[u].fail = trie[f].next[c];

            // odziedzicz wyjścia z fail-linka
            for(int id : trie[ trie[u].fail ].out)
                trie[u].out.push_back(id);
        }
    }

    cerr << "Automaton nodes: " << trie.size() << "\n";

    // Eksport do DOT (graphviz)
    if(argc >= 3){
        string out = argv[2];
        ofstream f(out);

        f << "digraph AC {\n"
          << "rankdir=LR;\n"
          << "node [shape=circle];\n";

        // węzły
        for(size_t i=0; i<trie.size(); i++){
            if(trie[i].out.empty())
                f << "n" << i << " [label=\"" << i << "\"];\n";
            else
                f << "n" << i << " [label=\"" << i
                  << "\\nout=" << trie[i].out.size()
                  << "\", style=filled, fillcolor=lightblue];\n";
        }

        // przejścia
        for(size_t i=0; i<trie.size(); i++){
            for(int c=0; c<5; c++){
                int v = trie[i].next[c];
                if(v==-1) continue;
                char ch = "ACGTN"[c];
                f << "n" << i << " -> n" << v
                  << " [label=\"" << ch << "\"];\n";
            }
        }

        // fail-linki
        for(size_t i=1; i<trie.size(); i++){
            f << "n" << i << " -> n" << trie[i].fail
              << " [style=dashed, color=gray, label=\"fail\"];\n";
        }

        f << "}\n";
        cerr << "Saved DOT graph to: " << out << "\n";
    }

    return 0;
}
