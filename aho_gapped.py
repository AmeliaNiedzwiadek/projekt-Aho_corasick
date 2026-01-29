# aho_gapped.py
import time
from collections import deque, defaultdict
import re
import os
import psutil
import matplotlib.pyplot as plt


class Node:
    """
    Reprezentuje węzeł automatu Aho–Corasick.

    Węzeł przechowuje przejścia znakowe, łącze porażki (failure link),
    listę wyjść (dopasowanych wzorców) oraz identyfikator numeryczny.
    """

    __slots__ = ("next", "fail", "out", "id")

    def __init__(self):
        """
        Inicjalizuje pusty węzeł automatu.
        """
        self.next = {}
        self.fail = None
        self.out = []
        self.id = None


class AhoCorasick:
    """
    Implementacja automatu Aho–Corasick używana jako silnik wyszukiwania
    fragmentów (seedów) wzorców z lukami.
    """

    def __init__(self):
        """
        Inicjalizuje pusty automat z węzłem korzenia.
        """
        self.root = Node()
        self._nodes = [self.root]

    def add(self, pattern, pat_id=None, meta=None):
        """
        Dodaje wzorzec (seed) do struktury trie.

        :param pattern: Sekwencja znaków dodawana do automatu.
        :type pattern: str
        :param pat_id: Identyfikator wzorca nadrzędnego.
        :type pat_id: int | None
        :param meta: Metadane powiązane z dopasowaniem (np. offset i długość).
        :type meta: dict | None
        :return: None
        :rtype: None
        """
        node = self.root
        for ch in pattern:
            if ch not in node.next:
                n = Node()
                node.next[ch] = n
                self._nodes.append(n)
            node = node.next[ch]
        node.out.append((pat_id, meta))

    def build(self):
        """
        Buduje automat poprzez wyznaczenie łączy porażki
        oraz propagację list wyjść.

        :return: None
        :rtype: None
        """
        q = deque()
        self.root.fail = self.root

        for ch, node in self.root.next.items():
            node.fail = self.root
            q.append(node)

        while q:
            r = q.popleft()
            for ch, u in r.next.items():
                q.append(u)
                v = r.fail
                while v is not self.root and ch not in v.next:
                    v = v.fail
                u.fail = v.next[ch] if ch in v.next else self.root
                u.out += u.fail.out

        for i, n in enumerate(self._nodes):
            n.id = i

    def search_generator(self, text):
        """
        Generator dopasowań seedów w tekście.

        :param text: Tekst wejściowy do przeszukania.
        :type text: str
        :yield: Krotki (indeks_końca, pat_id, meta).
        :rtype: tuple[int, int, dict]
        """
        node = self.root
        for idx, ch in enumerate(text):
            if ch not in "ACGTN":
                node = self.root
                continue
            while node is not self.root and ch not in node.next:
                node = node.fail
            node = node.next[ch] if ch in node.next else self.root
            if node.out:
                for (pat_id, meta) in node.out:
                    yield idx, pat_id, meta


_pattern_piece_re = re.compile(r'(?:[ACGTN]+)|(?:\.{1,}|\{(\d+)\})', re.IGNORECASE)


def parse_pattern_with_gaps(pat):
    """
    Parsuje wzorzec z lukami do listy tokenów.

    Obsługiwane są:
    - sekwencje ACGTN,
    - luki jako '.' (dowolny znak),
    - luki jako '{n}' (dokładnie n znaków).

    :param pat: Wzorzec wejściowy.
    :type pat: str
    :return: Lista tokenów w postaci (typ, wartość).
    :rtype: list[tuple[str, str|int]]
    :raises ValueError: Gdy składnia nawiasów jest niepoprawna.
    """
    s = pat.strip().upper()
    tokens = []
    i = 0

    while i < len(s):
        if s[i] in 'ACGTN':
            j = i
            while j < len(s) and s[j] in 'ACGTN':
                j += 1
            tokens.append(('seq', s[i:j]))
            i = j

        elif s[i] == '.':
            j = i
            while j < len(s) and s[j] == '.':
                j += 1
            tokens.append(('gap', j - i))
            i = j

        elif s[i] == '{':
            j = i + 1
            while j < len(s) and s[j] != '}':
                j += 1
            if j < len(s) and s[j] == '}':
                num = int(s[i+1:j])
                tokens.append(('gap', num))
                i = j + 1
            else:
                raise ValueError("Bad pattern brace: " + s)

        elif s[i] == 'N':
            tokens.append(('gap', 1))
            i += 1

        else:
            i += 1

    return tokens


def build_seeds_from_pattern(tokens, min_seed_len=3):
    """
    Wyznacza seedy (ciągłe fragmenty) ze wzorca z lukami.

    :param tokens: Tokeny wzorca.
    :type tokens: list[tuple[str, str|int]]
    :param min_seed_len: Minimalna długość seeda.
    :type min_seed_len: int
    :return: Lista seedów i ich offsetów.
    :rtype: list[tuple[str, int]]
    """
    seeds = []
    offset = 0

    for typ, val in tokens:
        if typ == 'seq':
            if len(val) >= min_seed_len:
                seeds.append((val, offset))
            offset += len(val)
        else:
            offset += val

    return seeds


def prepare_automaton_from_patterns(patterns, min_seed_len=3):
    """
    Buduje automat seedów oraz metadane wzorców z lukami.

    :param patterns: Lista wzorców wejściowych.
    :type patterns: list[str]
    :param min_seed_len: Minimalna długość seeda.
    :type min_seed_len: int
    :return: Automat oraz metadane wzorców.
    :rtype: tuple[AhoCorasick, dict]
    """
    ac = AhoCorasick()
    meta = {}

    for pid, pat in enumerate(patterns):
        toks = parse_pattern_with_gaps(pat)
        total_len = sum((len(v) if t == 'seq' else v) for t, v in toks)
        meta[pid] = {'pattern': pat, 'tokens': toks, 'length': total_len}

        seeds = build_seeds_from_pattern(toks, min_seed_len=min_seed_len)

        if not seeds:
            for t, v in toks:
                if t == 'seq' and v:
                    ac.add(v, pat_id=pid, meta={'offset': 0})
                    break
        else:
            for seed, off in seeds:
                ac.add(seed, pat_id=pid, meta={'offset': off})

    ac.build()
    return ac, meta


def verify_pattern_at(text, end_pos, meta_entry):
    """
    Funkcja zarezerwowana do weryfikacji dopasowania wzorca.

    :raises RuntimeError: Zawsze — funkcja nie jest używana.
    """
    raise RuntimeError("Use verify_pattern_with_anchor (needs seed_offset provided)")


def verify_pattern_with_anchor(text, match_end_idx, seed_offset, pat_meta):
    """
    Funkcja zarezerwowana do weryfikacji dopasowania z kotwicą.

    :raises RuntimeError: Zawsze — funkcja nie jest używana.
    """
    raise RuntimeError("Higher-level search uses explicit seed_len from meta mapping.")


def load_fasta(path):
    """
    Wczytuje sekwencję z pliku FASTA.

    :param path: Ścieżka do pliku.
    :type path: str
    :return: Połączona sekwencja.
    :rtype: str
    """
    seqs = []
    with open(path, "r") as f:
        cur = []
        for line in f:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur).upper())
                    cur = []
            else:
                cur.append(line.strip())
        if cur:
            seqs.append("".join(cur).upper())
    return "".join(seqs)


def load_patterns_txt(path):
    """
    Wczytuje wzorce z pliku tekstowego.

    :param path: Ścieżka do pliku.
    :type path: str
    :return: Lista wzorców.
    :rtype: list[str]
    """
    pats = []
    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if s:
                pats.append(s)
    return pats


def search_with_gaps(text, patterns, min_seed_len=3):
    """
    Wyszukuje wzorce z lukami w tekście.

    :param text: Sekwencja wejściowa.
    :type text: str
    :param patterns: Lista wzorców z lukami.
    :type patterns: list[str]
    :param min_seed_len: Minimalna długość seeda.
    :type min_seed_len: int
    :return: Dopasowania pogrupowane według wzorca.
    :rtype: dict[int, list[tuple[int, int]]]
    """
    matches_by_pat = defaultdict(list)

    ac = AhoCorasick()
    meta = {}

    for pid, pat in enumerate(patterns):
        toks = parse_pattern_with_gaps(pat)
        total_len = sum((len(v) if t == 'seq' else v) for t, v in toks)
        meta[pid] = {'pattern': pat, 'tokens': toks, 'length': total_len}

        seeds = build_seeds_from_pattern(toks, min_seed_len=min_seed_len)

        if not seeds:
            for t, v in toks:
                if t == 'seq' and v:
                    ac.add(v, pat_id=pid,
                           meta={'offset': 0, 'seed_len': len(v)})
                    break
        else:
            for seed, off in seeds:
                ac.add(seed, pat_id=pid,
                       meta={'offset': off, 'seed_len': len(seed)})

    ac.build()

    for idx, pid, meta_info in ac.search_generator(text):
        seed_offset = meta_info['offset']
        seed_len = meta_info['seed_len']
        pat_meta = meta[pid]
        pattern_len = pat_meta['length']

        pattern_start = idx - (seed_len - 1) - seed_offset
        pattern_end = pattern_start + pattern_len

        if pattern_start < 0 or pattern_end > len(text):
            continue

        ok = True
        tpos = pattern_start

        for typ, val in pat_meta['tokens']:
            if typ == 'seq':
                if text[tpos:tpos+len(val)] != val:
                    ok = False
                    break
                tpos += len(val)
            else:
                tpos += val

        if ok:
            matches_by_pat[pid].append((pattern_start, pattern_end))

    return matches_by_pat


def compare_sequences(s1, s2):
    """
    Porównuje dwie sekwencje znak po znaku.

    :param s1: Pierwsza sekwencja.
    :type s1: str
    :param s2: Druga sekwencja.
    :type s2: str
    :return: Lista różnic.
    :rtype: list[tuple]
    """
    mismatches = []
    L = min(len(s1), len(s2))

    for i in range(L):
        if s1[i] != s2[i]:
            mismatches.append((i, s1[i], s2[i]))

    if len(s1) != len(s2):
        mismatches.append(('len_diff', len(s1), len(s2)))

    return mismatches


def get_mem_usage_mb():
    """
    Zwraca bieżące zużycie pamięci procesu.

    :return: Zużycie pamięci w MB.
    :rtype: float
    """
    p = psutil.Process(os.getpid())
    return p.memory_info().rss / (1024 * 1024)


def run_benchmark_gapped(fasta_path, patterns_path, min_seed_len=3):
    """
    Uruchamia benchmark wyszukiwania wzorców z lukami.

    :param fasta_path: Plik FASTA.
    :type fasta_path: str
    :param patterns_path: Plik wzorców.
    :type patterns_path: str
    :param min_seed_len: Minimalna długość seeda.
    :type min_seed_len: int
    :return: Statystyki benchmarku.
    :rtype: dict
    """
    text = load_fasta(fasta_path)
    pats = load_patterns_txt(patterns_path)

    t0 = time.perf_counter()
    cmem0 = get_mem_usage_mb()
    matches = search_with_gaps(text, pats, min_seed_len=min_seed_len)
    t1 = time.perf_counter()
    cmem1 = get_mem_usage_mb()

    total_matches = sum(len(v) for v in matches.values())

    return {
        'time': t1 - t0,
        'mem_start_mb': cmem0,
        'mem_end_mb': cmem1,
        'matches': total_matches,
        'per_pattern': {p: len(matches.get(i, [])) for i, p in enumerate(pats)}
    }

if __name__ == "__main__":
    import sys
    if len(sys.argv) == 3:
        fasta = sys.argv[1]
        patterns = sys.argv[2]
    else:
        fasta = "human.fasta"
        patterns = "patterns.txt"
        print("[Spyder mode] Using defaults:", fasta, patterns)
    res = run_benchmark_gapped(fasta, patterns, min_seed_len=3)
    plt.bar([0], [res['time']])
    plt.ylabel("Search time (s)")
    plt.savefig("quick_time.png")
    print("Saved quick_time.png")
