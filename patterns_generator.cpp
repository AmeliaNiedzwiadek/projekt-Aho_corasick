/**
 * @file pattern_generator.cpp
 * @author Maria Ławniczak (Nr Indeksu: 268544)
 * @brief Generator zestawów testowych wzorców DNA z uwzględnieniem luk (wildcards)
 * Narzędzie ekstrahuje fragmenty z plików FASTA i wprowadza do nich znaki '.' (dowolny nukleotyd)
 * @date 2026-01-25
 */

#include <bits/stdc++.h>
using namespace std;

/**
 * @brief Parsuje plik w formacie FASTA do pojedynczego ciągu znaków
 * Funkcja ignoruje nagłówki (linie zaczynające się od '>') oraz wszelkie białe znaki,
 * normalizując sekwencję do wielkich liter.
 * * @param path Ścieżka do pliku wejściowego .fasta lub .fa
 * @return string Połączona sekwencja nukleotydowa
 */
string load_fasta(const string &path){
    ifstream in(path);
    if(!in){
        cerr << "Cannot open " << path << "\n";
        exit(1);
    }
    string line, acc, result;
    while(getline(in, line)){
        if(!line.empty() && line[0] == '>'){
            // Obsługa wielu rekordów w jednym pliku FASTA
            if(!acc.empty()){
                result += acc;
                acc.clear();
            }
        } else {
            // Czyszczenie linii z białych znaków i dodawanie do akumulatora
            for(char c: line)
                if(!isspace((unsigned char)c))
                    acc.push_back(toupper(c));
        }
    }
    if(!acc.empty()) result += acc;
    return result;
}

/**
 * @brief Losowo wprowadza znaki maskowania (dziury) do wzorca
 * Zamienia wybraną liczbę nukleotydów na kropki ('.') zgodnie z zadaną frakcją
 * * @param s Oryginalna sekwencja DNA
 * @param gap_frac Prawdopodobieństwo/udział luk w sekwencji (0.0 - 1.0)
 * @param rng Referencja do generatora liczb losowych
 * @return string Zmodyfikowany wzorzec z dziurami
 */
string add_gaps(string s, double gap_frac, mt19937 &rng){
    if(gap_frac <= 0)
        return s;

    // Obliczanie docelowej liczby luk (minimum 1, jeśli frakcja > 0)
    int to_gap = max(1, (int)round(s.size()*gap_frac));
    uniform_int_distribution<int> dist(0, s.size()-1);

    // Losowanie unikalnych pozycji dla dziur
    unordered_set<int> used;
    while((int)used.size() < to_gap)
        used.insert(dist(rng));

    for(int idx: used)
        s[idx] = '.';

    return s;
}

/**
 * @brief Zapisuje wygenerowane wzorce do pliku tekstowego (jeden wzorzec w linii)
 * * @param pats Wektor wygenerowanych wzorców
 * @param path Ścieżka docelowa zapisu
 */
void save_patterns(const vector<string> &pats, const string &path){
    ofstream out(path);
    for(const string &p: pats)
        out << p << "\n";
}

/**
 * @brief Punkt wejścia programu - zarządza generowaniem zestawów testowych
 * Obsługuje 3 predefiniowane konfiguracje wielkości i długości wzorców
 * * @param argc Liczba argumentów (wymagane: fasta, prefix)
 * @param argv Tablica argumentów (fasta, prefix, opcjonalnie gap_fraction)
 */
int main(int argc, char **argv){
    if(argc < 3){
        cerr << "Usage: " << argv[0] << " <fasta> <prefix> [gap_fraction]\n";
        return 1;
    }

    string fasta = argv[1];
    string prefix = argv[2];
    // Domyślna frakcja dziur wynosi 0.2 (20%), jeśli nie podano inaczej
    double gap_frac = (argc >= 4) ? stod(argv[3]) : 0.2;

    string text = load_fasta(fasta);
    // Użycie stałego ziarna dla powtarzalności wyników testowych
    mt19937 rng(123456); 

    // Konfiguracje eksperymentów: {liczba_wzorców, długość_pojedynczego_wzorca}
    vector<pair<int,int>> configs = {
        {10, 10},
        {50, 12},
        {200, 20}
    };

    // Generowanie i zapisywanie każdego zestawu do oddzielnego pliku
    for(auto &cfg : configs){
        int count = cfg.first;
        int len   = cfg.second;

        vector<string> pats;
        uniform_int_distribution<int> pos(0, text.size()-len-1);

        for(int i=0; i<count; i++){
            // Pobranie losowego fragmentu z wczytanej sekwencji FASTA
            int p = pos(rng);
            string s = text.substr(p, len);

            // Aplikacja maskowania
            s = add_gaps(s, gap_frac, rng);

            pats.push_back(s);
        }

        // Nazewnictwo plików wynikowych: prefix_liczba.txt
        string fname = prefix + "_" + to_string(count) + ".txt";
        save_patterns(pats, fname);

        cerr << "[OK] Saved " << fname << " (Patterns: " << count << ", Len: " << len << ")\n";
    }

    return 0;
}
