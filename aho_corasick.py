import time
from collections import deque

# ============================================================
# Node class
# ============================================================

class Node:
    __slots__ = ("next", "fail", "out", "id")

    def __init__(self):
        self.next = {}     # transitions: char -> Node
        self.fail = None   # failure link
        self.out = []      # patterns ending here
        self.id = None     # node ID for visualization


# ============================================================
# Aho–Corasick automaton
# ============================================================

class AhoCorasick:
    def __init__(self):
        self.root = Node()
        self._nodes = [self.root]

    def add(self, pattern, pat_id=None):
        node = self.root
        for ch in pattern:
            if ch not in node.next:
                new = Node()
                node.next[ch] = new
                self._nodes.append(new)
            node = node.next[ch]
        node.out.append(pat_id if pat_id is not None else pattern)

    def build(self):
        queue = deque()
        self.root.fail = self.root

        for ch, nxt in self.root.next.items():
            nxt.fail = self.root
            queue.append(nxt)

        # BFS build
        while queue:
            r = queue.popleft()

            for ch, u in r.next.items():
                queue.append(u)
                v = r.fail

                while v is not self.root and ch not in v.next:
                    v = v.fail

                u.fail = v.next[ch] if ch in v.next else self.root
                u.out += u.fail.out

        # Assign numeric IDs
        for i, n in enumerate(self._nodes):
            n.id = i

    def search(self, text, yield_matches=False):
        node = self.root
        total = 0
        matches = [] if yield_matches else None
        idx = 0

        for ch in text:
            if ch not in "ACGTN":   # ignore invalid DNA chars
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
    nodes = ac._nodes
    N = min(len(nodes), limit_nodes)

    with open(filename, "w") as f:
        f.write("digraph aho {\n")
        f.write("  rankdir=LR;\n")
        f.write("  node [shape=circle,fontname=Helvetica];\n")

        # nodes
        for i in range(N):
            if nodes[i].out:
                f.write(f'  n{i} [label="{i}\\nout={len(nodes[i].out)}",style=filled,fillcolor=lightblue];\n')
            else:
                f.write(f'  n{i} [label="{i}"];\n')

        # transitions
        for i in range(N):
            for ch, nxt in nodes[i].next.items():
                if nxt.id < N:
                    f.write(f'  n{i} -> n{nxt.id} [label="{ch}"];\n')

        # failure links
        for i in range(N):
            fail = nodes[i].fail
            if fail and fail.id < N and nodes[i] is not ac.root:
                f.write(f'  n{i} -> n{fail.id} [style=dashed,color=gray,label="f"];\n')

        f.write("}\n")

    print(f"[OK] DOT automaton saved to {filename}")


# ============================================================
# Benchmark
# ============================================================

def benchmark(fasta_path, patterns_path):
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


if __name__ == "__main__":
    import sys

    if len(sys.argv) != 3:
        print("Usage:")
        print("   python3 aho_corasick.py <sequence.fasta> <patterns.txt>")
        sys.exit(1)

    fasta = sys.argv[1]
    pats = sys.argv[2]

    print(f"Sequence file : {fasta}")
    print(f"Patterns file : {pats}\n")

    benchmark(fasta, pats)


