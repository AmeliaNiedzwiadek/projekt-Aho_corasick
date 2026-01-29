/**
 * @file aho_corasick.cpp
 * @author Maria Ławniczak (Nr Indeksu: 268544)
 * @brief Implementacja algorytmu Aho-Corasick dla sekwencji DNA
 * @date 2026-01-25
 */

#include <bits/stdc++.h>

using namespace std;

/**
 * @brief Mapuje znaki alfabetu DNA na indeksy tablicy (0-4)
 * Obsługuje A, C, G, T oraz N (jako błąd lub nieznany nukleotyd)
 * * @param c Znak do zmapowania
 * @return int Indeks z zakresu 0-4
 */
int char_idx(char c) {
    switch (toupper(c)) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        case 'N': return 4;
        default:  return 4; 
    }
}

/**
 * @brief Wczytuje wzorce z pliku tekstowego, czyszcząc je z białych znaków
 * * @param path Ścieżka do pliku
 * @return vector<string> Lista oczyszczonych wzorców
 */
vector<string> load_patterns(const string& path) {
    ifstream in(path);
    if (!in) {
        cerr << "Blad: Nie mozna otworzyc pliku " << path << "\n";
        exit(1);
    }
    string line;
    vector<string> out;
    while (getline(in, line)) {
        string clean;
        for (char c : line) {
            if (!isspace((unsigned char)c))
                clean.push_back(toupper(c));
        }
        if (!clean.empty()) out.push_back(clean);
    }
    return out;
}

/**
 * @brief Struktura reprezentująca pojedynczy węzeł w automacie Aho-Corasick
 */
struct Node {
    array<int, 5> next;  // Przejścia do kolejnych stanów dla alfabetu {A, C, G, T, N}
    int fail;            // Wskaźnik funkcji porażki
    vector<int> out;     ///Lista identyfikatorów wzorców, które kończą się w tym stanie

    Node() {
        next.fill(-1);
        fail = 0;
    }
};

int main(int argc, char** argv) {
    if (argc < 2) {
        cerr << "Uzycie: " << argv[0] << " <wzorce.txt> [eksport.dot]\n";
        return 1;
    }

    // Inicjalizacja danych
    vector<string> pats = load_patterns(argv[1]);
    vector<Node> trie(1); // Stan 0 to root

    // Budowa drzewa Trie
    for (size_t pid = 0; pid < pats.size(); pid++) {
        int v = 0;
        for (char c : pats[pid]) {
            int id = char_idx(c);
            if (trie[v].next[id] == -1) {
                trie[v].next[id] = trie.size();
                trie.emplace_back();
            }
            v = trie[v].next[id];
        }
        trie[v].out.push_back(pid); // Oznaczamy koniec wzorca w danym węźle
    }

    // Budowa funkcji porażki
    queue<int> q;

    // Inicjalizacja poziomu 1 (bezpośredni sąsiedzi roota)
    for (int c = 0; c < 5; c++) {
        int nxt = trie[0].next[c];
        if (nxt != -1 && nxt != 0) {
            trie[nxt].fail = 0;
            q.push(nxt);
        } else {
            trie[0].next[c] = 0; // Jeśli brak przejścia z roota, wracamy do roota
        }
    }

    // Przetwarzanie kolejnych poziomów drzewa
    while (!q.empty()) {
        int r = q.front();
        q.pop();

        for (int c = 0; c < 5; c++) {
            int u = trie[r].next[c];
            if (u == -1 || u == 0) continue;

            q.push(u);
            int f = trie[r].fail;

            // Szukamy najdłuższego właściwego sufiksu, który jest w Trie
            while (trie[f].next[c] == -1) {
                f = trie[f].fail;
            }

            trie[u].fail = trie[f].next[c];
    
            // Jeśli węzeł, do którego prowadzi fail-link, jest końcem wzorca, to obecny węzeł również "zawiera" ten wzorzec
            for (int id : trie[trie[u].fail].out) {
                trie[u].out.push_back(id);
            }
        }
    }

    cerr << "Statystyki automatu:\n";
    cerr << " - Liczba wzorcow: " << pats.size() << "\n";
    cerr << " - Liczba wezlow (stanow): " << trie.size() << "\n";

    // Eksport do formatu DOT 
    if (argc >= 3) {
        ofstream f(argv[2]);
        if (!f) {
            cerr << "Nie mozna utworzyc pliku .dot!\n";
            return 1;
        }

        f << "digraph AC {\n"
          << "  rankdir=LR;\n"
          << "  node [shape=circle];\n";

        for (size_t i = 0; i < trie.size(); i++) {
            if (trie[i].out.empty()) {
                f << "  n" << i << " [label=\"" << i << "\"];\n";
            } else {
                f << "  n" << i << " [label=\"" << i
                  << "\\n(Pats: " << trie[i].out.size()
                  << ")\", style=filled, fillcolor=lightblue];\n";
            }
        }

        // Krawędzie drzewa (przejścia)
        for (size_t i = 0; i < trie.size(); i++) {
            for (int c = 0; c < 5; c++) {
                int v = trie[i].next[c];
                if (v != -1 && v != 0) {
                    f << "  n" << i << " -> n" << v << " [label=\"" << "ACGTN"[c] << "\"];\n";
                }
            }
        }

        // Krawędzie fail-linków (przerywane)
        for (size_t i = 1; i < trie.size(); i++) {
            f << "  n" << i << " -> n" << trie[i].fail
              << " [style=dashed, color=red, label=\"fail\"];\n";
        }

        f << "}\n";
        cerr << "Zapisano graf do: " << argv[2] << "\n";
    }

    return 0;
}
