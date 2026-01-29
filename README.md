Projekt dotyczy implementacji i analizy algorytmów wielowzorcowego wyszukiwania sekwencji w danych genomowych.
Zaimplementowano klasyczny automat Aho–Corasick oraz jego rozszerzenie obsługujące wzorce z lukami (gapped patterns).
Rozwiązanie umożliwia wydajne przeszukiwanie dużych sekwencji DNA przy użyciu wielu wzorców jednocześnie.
Implementacje zostały wykonane w Pythonie oraz C++. Projekt obejmuje również narzędzia pomocnicze do generowania danych testowych, wykrywania mutacji oraz wykonywania benchmarków czasowo-pamięciowych.

Projekt obejmuje:
-budowę automatu Aho–Corasick
-wielowzorcowe wyszukiwanie w sekwencjach DNA
-obsługę wzorców z lukami (wildcards / gaps)
-podział wzorców z lukami na fragmenty (seedy)
-weryfikację pełnego dopasowania po trafieniu seeda
-generowanie wzorców testowych z genomu
-porównywanie sekwencji i wykrywanie mutacji
-pomiar czasu wykonania
-pomiar zużycia pamięci
-eksport struktury automatu do grafu (DOT)

System obsługuje:
-wzorce dokładne (ciągłe)
-wzorce z symbolami wieloznacznymi
-luki o długości stałej
-luki o długości zadanej liczbowo

