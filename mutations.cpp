/**
 * @file dna_comparator.cpp
 * @author Maria Ławniczak (Nr Indeksu: 268544)
 * @brief Narzędzie do analizy różnic między dwiema sekwencjami DNA
 * Wykrywa SNP (substytucje), insercje, delecje oraz mutacje złożone
 * @date 2026-01-25
 */

#include <bits/stdc++.h>
using namespace std;

/**
 * @brief Funkcja wczytująca dane wejściowe
 * Obsługuje pliki (np. FASTA/TXT) oraz surowe ciągi znaków
 * Automatycznie usuwa białe znaki i normalizuje tekst do wielkich liter
 * * @param path Ścieżka do pliku lub surowa sekwencja
 * @return string Oczyszczona sekwencja nukleotydowa
 */
string load_text(const string &path){
    ifstream in(path);
    if(!in)
        return ""; // Zwraca pusty string, jeśli ścieżka nie jest poprawnym plikiem

    string line, result;
    while(getline(in, line)){
        // Pomijamy nagłówki FASTA (linie zaczynające się od '>')
        if(!line.empty() && line[0] == '>') continue; 
        
        for(char c: line)
            if(!isspace((unsigned char)c))
                result.push_back(toupper(c));
    }
    return result;
}

/**
 * @brief Algorytm porównujący sekwencje przy użyciu dwóch wskaźników
 * Wykorzystuje zachłanną heurystykę (look-ahead) do klasyfikacji zmian
 * * @param a Sekwencja referencyjna (oryginalna)
 * @param b Sekwencja zapytania (zmutowana)
 * @return vector<string> Lista opisowa wykrytych różnic
 */
vector<string> compare_seqs(const string &a, const string &b){
    vector<string> result;

    int i = 0, j = 0; // Wskaźniki pozycji: i dla sekwencji A, j dla sekwencji B
    int n = a.size(), m = b.size();

    while(i < n && j < m){
        if(a[i] == b[j]){
            // Znaki identyczne - przesuwamy oba wskaźniki
            i++;
            j++;
            continue;
        }

        // Analiza typu mutacji

        // SNP / Substytucja - jeśli znaki się różnią, ale następne w obu sekwencjach pasują do siebie
        if(i+1 < n && j+1 < m && a[i+1] == b[j+1]){
            result.push_back(
                "SNP at pos " + to_string(i) +
                ": " + a[i] + " -> " + b[j]
            );
            i++;
            j++;
        }

        // Delecja (usunięcie nukleotydu w sekwencji B) - jeśli następny znak w A pasuje do obecnego w B
        else if(i+1 < n && a[i+1] == b[j]){
            result.push_back(
                "Deletion at pos " + to_string(i) +
                ": removed " + a[i]
            );
            i++;
        }

        // Insercja (wstawienie nukleotydu w sekwencji B)- jeśli obecny znak w A pasuje do następnego w B
        else if(j+1 < m && a[i] == b[j+1]){
            result.push_back(
                "Insertion at pos " + to_string(i) +
                ": inserted " + b[j]
            );
            j++;
        }

        // Zmiana złożona (Complex mutation)- jeśli prosta heurystyka zawodzi - klasyfikujemy jako zmianę grupową
        else {
            result.push_back(
                "Complex mutation near pos A=" +
                to_string(i) + " B=" + to_string(j)
            );
            i++;
            j++;
        }
    }

    // Obsługa końcówek sekwencji - jeśli jedna sekwencja jest dłuższa od drugiej na samym końcu
    while(i < n){
        result.push_back("Deletion at end: removed " + string(1, a[i]));
        i++;
    }
    while(j < m){
        result.push_back("Insertion at end: inserted " + string(1, b[j]));
        j++;
    }

    return result;
}

/**
 * @brief Punkt wejścia programu
 * Obsługuje argumenty wiersza poleceń: ./program plik1 plik2
 */
int main(int argc, char **argv){
    if(argc < 3){
        cerr << "Sposób użycia: " << argv[0] << " <seqA|fileA> <seqB|fileB>\n";
        return 1;
    }

    // Próba wczytania danych jako plików
    string A = load_text(argv[1]);
    string B = load_text(argv[2]);

    // Jeśli load_text zwrócił puste (brak pliku), traktujemy argumenty jako surowe DNA
    if(A.empty()){ A = argv[1]; for(char &c : A) c = toupper(c); }
    if(B.empty()){ B = argv[2]; for(char &c : B) c = toupper(c); }

    // Wykonanie porównania
    auto diffs = compare_seqs(A, B);

    // Prezentacja wyników
    cout << "Detected differences (" << diffs.size() << "):\n";
    if(diffs.empty()){
        cout << " - No differences found. Sequences are identical.\n";
    } else {
        for(const auto &d : diffs)
            cout << " - " << d << "\n";
    }

    return 0;
}

