
// Porównywacz sekwencji DNA. Wykrywa różnice typu:
//   - SNP (pojedyncza zmiana znaku)
//   - delecja
//   - insercja
//   - zmiana złożona (fallback)
//
// Może działać na:
//   - surowych stringach
//   - plikach FASTA / TXT (ignoruje whitespace)

#include <bits/stdc++.h>
using namespace std;

// Wczytywanie sekwencji z pliku tekstowego / FASTA
// Ignoruje whitespace, zamienia na uppercase
string load_text(const string &path){
    ifstream in(path);
    if(!in)
        return ""; // pliku nie ma = zwróci pusty

    string line, result;
    while(getline(in, line)){
        for(char c: line)
            if(!isspace((unsigned char)c))
                result.push_back(toupper(c));
    }
    return result;
}

// Porównanie dwóch sekwencji litera po literze
//   - jeśli różne - próbujemy ustalić typ zmiany:
//       SNP: A!=G, a następne znaki pasują
//       delecja: litera z A pominięta
//       insercja: litera dodatkowa w B
vector<string> compare_seqs(const string &a, const string &b){
    vector<string> result;

    int i = 0, j = 0;
    int n = a.size(), m = b.size();

    while(i < n && j < m){
        if(a[i] == b[j]){
            // zgodne litery
            i++;
            j++;
            continue;
        }

        // 1) SNP (substitution)
        if(i+1 < n && j+1 < m && a[i+1] == b[j+1]){
            result.push_back(
                "SNP at pos " + to_string(i) +
                ": " + a[i] + " -> " + b[j]
            );
            i++;
            j++;
        }

        // 2) Delecja
        else if(i+1 < n && a[i+1] == b[j]){
            result.push_back(
                "Deletion at pos " + to_string(i) +
                ": removed " + a[i]
            );
            i++;
        }

        // 3) Insercja
        else if(j+1 < m && a[i] == b[j+1]){
            result.push_back(
                "Insertion at pos " + to_string(i) +
                ": inserted " + b[j]
            );
            j++;
        }

        // 4) Zmiana złożona
        else {
            result.push_back(
                "Complex mutation near pos A=" +
                to_string(i) + " B=" + to_string(j)
            );
            i++;
            j++;
        }
    }

    // resztki
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

// MAIN
int main(int argc, char **argv){
    if(argc < 3){
        cerr << "Usage: " << argv[0] << " <seqA|fileA> <seqB|fileB>\n";
        return 1;
    }

    // Spróbuj wczytać z pliku — jeśli puste, użyj argumentu jako sekwencji
    string A = load_text(argv[1]);
    string B = load_text(argv[2]);

    if(A.empty()){ A = argv[1]; for(char &c : A) c = toupper(c); }
    if(B.empty()){ B = argv[2]; for(char &c : B) c = toupper(c); }

    auto diffs = compare_seqs(A, B);

    cout << "Differences (" << diffs.size() << "):\n";
    for(const auto &d : diffs)
        cout << " - " << d << "\n";

    return 0;
}
