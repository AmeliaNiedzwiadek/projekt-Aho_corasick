import time
from collections import deque

# Node class

class Node:
    """
    Reprezentuje węzeł automatu Aho–Corasick.

    Węzeł przechowuje przejścia znakowe, łącze porażki (failure link),
    listę identyfikatorów wzorców kończących się w tym węźle oraz
    identyfikator numeryczny używany przy wizualizacji grafu.
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


# Aho–Corasick automaton

class AhoCorasick:
    """
    Implementacja automatu Aho–Corasick do wielowzorcowego
    wyszukiwania podciągów w tekście.

    Automat umożliwia dodawanie wzorców, budowę struktury
    przejść oraz efektywne przeszukiwanie tekstu w czasie liniowym.
    """

    def __init__(self):
        """
        Tworzy pusty automat z węzłem korzenia.
        """
        self.root = Node()
        self._nodes = [self.root]

    def add(self, pattern, pat_id=None):
        """
        Dodaje wzorzec do struktury trie automatu.

        :param pattern: Wzorzec znakowy do dodania.
        :type pattern: str
        :param pat_id: Identyfikator wzorca zapisywany w liście wyjść.
                       Jeżeli None, zapisywany jest sam wzorzec.
        :type pat_id: object | None
        :return: None
        :rtype: None
        """
        node = self.root
        for ch in pattern:
            if ch not in node.next:
                new = Node()
                node.next[ch] = new
                self._nodes.append(new)
            node = node.next[ch]
        node.out.append(pat_id if pat_id is not None else pattern)

    def build(self):
        """
        Buduje automat Aho–Corasick na podstawie dodanych wzorców.

        Operacja obejmuje:
        - wyznaczenie łączy porażki (failure links),
        - propagację list wyjść,
        - przypisanie identyfikatorów numerycznych węzłom.

        :return: None
        :rtype: None
        """
        queue = deque()
        self.root.fail = self.root

        for ch, nxt in self.root.next.items():
            nxt.fail = self.root
            queue.append(nxt)

        while queue:
            r = queue.popleft()

            for ch, u in r.next.items():
                queue.append(u)
                v = r.fail

                while v is not self.root and ch not in v.next:
                    v = v.fail

                u.fail = v.next[ch] if ch in v.next else self.root
                u.out += u.fail.out

        for i, n in enumerate(self._nodes):
            n.id = i

    def search(self, text, yield_matches=False):
        """
        Przeszukuje tekst przy użyciu zbudowanego automatu.

        Znaki spoza alfabetu ``ACGTN`` powodują powrót do korzenia
        automatu i są pomijane w dopasowaniach.

        :param text: Tekst wejściowy do przeszukania.
        :type text: str
        :param yield_matches: Jeżeli True, zwracana jest również lista
                              szczegółowych dopasowań.
        :type yield_matches: bool
        :return: Liczba dopasowań lub krotka (liczba, lista_dopasowań).
        :rtype: int | tuple[int, list[tuple[int, object]]]
        """
        node = self.root
        total = 0
        matches = [] if yield_matches else None
        idx = 0

        for ch in text:
            if ch not in "ACGTN":
                node = self.root
                idx += 1
                continue

            while node is not self.root and ch not in node.next:
                node = node.fail

            node = node.next[ch] if ch in node.next else self.root

            if node.out:
                total += len(node.out)
                if yield_matches:
                    for o in node.out:
                        matches.append((idx, o))

            idx += 1

        return (total, matches) if yield_matches else total


def load_fasta(path):
    """
    Wczytuje sekwencję z pliku FASTA i zwraca ją jako pojedynczy łańcuch.

    Linie nagłówkowe (rozpoczynające się od znaku ``>``) są pomijane.

    :param path: Ścieżka do pliku FASTA.
    :type path: str
    :return: Połączona sekwencja w wielkich literach.
    :rtype: str
    """
    seqs = []
    with open(path, "r") as f:
        buf = []
        for line in f:
            if line.startswith(">"):
                if buf:
                    seqs.append("".join(buf).upper())
                    buf = []
            else:
                buf.append(line.strip())
        if buf:
            seqs.append("".join(buf).upper())
    return "".join(seqs)


def load_patterns(path):
    """
    Wczytuje wzorce z pliku tekstowego.

    Pomijane są linie puste oraz zawierające znaki spoza alfabetu ``ACGTN``.

    :param path: Ścieżka do pliku ze wzorcami.
    :type path: str
    :return: Lista poprawnych wzorców.
    :rtype: list[str]
    """
    pats = []
    with open(path, "r") as f:
        for line in f:
            s = line.strip().upper()
            if not s:
                continue
            if any(ch not in "ACGTN" for ch in s):
                continue
            pats.append(s)
    return pats


def export_dot(ac, filename="automaton.dot", limit_nodes=1500):
    """
    Eksportuje strukturę automatu do pliku DOT (Graphviz).

    :param ac: Zbudowany automat Aho–Corasick.
    :type ac: AhoCorasick
    :param filename: Nazwa pliku wyjściowego.
    :type filename: str
    :param limit_nodes: Maksymalna liczba węzłów do zapisania.
    :type limit_nodes: int
    :return: None
    :rtype: None
    """
    nodes = ac._nodes
    N = min(len(nodes), limit_nodes)

    with open(filename, "w") as f:
        f.write("digraph aho {\n")
        f.write("  rankdir=LR;\n")
        f.write("  node [shape=circle,fontname=Helvetica];\n")

        for i in range(N):
            if nodes[i].out:
                f.write(f'  n{i} [label="{i}\\nout={len(nodes[i].out)}",style=filled,fillcolor=lightblue];\n')
            else:
                f.write(f'  n{i} [label="{i}"];\n')

        for i in range(N):
            for ch, nxt in nodes[i].next.items():
                if nxt.id < N:
                    f.write(f'  n{i} -> n{nxt.id} [label="{ch}"];\n')

        for i in range(N):
            fail = nodes[i].fail
            if fail and fail.id < N and nodes[i] is not ac.root:
                f.write(f'  n{i} -> n{fail.id} [style=dashed,color=gray,label="f"];\n')

        f.write("}\n")

    print(f"[OK] DOT automaton saved to {filename}")


def benchmark(fasta_path, patterns_path):
    """
    Wykonuje pełny benchmark działania automatu Aho–Corasick.

    :param fasta_path: Ścieżka do pliku FASTA.
    :type fasta_path: str
    :param patterns_path: Ścieżka do pliku ze wzorcami.
    :type patterns_path: str
    :return: None
    :rtype: None
    """
    print("=== Benchmark: Aho–Corasick ===")
    print("Loading FASTA...")
    text = load_fasta(fasta_path)
    print(f"Text length: {len(text)}")

    print("Loading patterns...")
    pats = load_patterns(patterns_path)
    print(f"Patterns: {len(pats)}")

    ac = AhoCorasick()

    print("Adding patterns...")
    for i, p in enumerate(pats):
        ac.add(p, i)

    print("Building automaton...")
    t0 = time.perf_counter()
    ac.build()
    t1 = time.perf_counter()
    print(f"Build time: {t1 - t0:.3f} s; nodes = {len(ac._nodes)}")

    print("Searching...")
    t0 = time.perf_counter()
    matches = ac.search(text)
    t1 = time.perf_counter()
    print(f"Search time: {t1 - t0:.3f} s")
    print(f"Total matches: {matches}")

    export_dot(ac, "automaton.dot")


