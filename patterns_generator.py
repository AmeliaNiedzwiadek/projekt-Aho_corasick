# patterns_generator.py
import random
import argparse

def load_fasta(path):
    seqs = []
    with open(path) as f:
        cur = []
        for line in f:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur).upper())
                    cur=[]
            else:
                cur.append(line.strip())
        if cur:
            seqs.append("".join(cur).upper())
    return "".join(seqs)

def add_gaps(seq, gap_fraction=0.0):
    """Zamienia część znaków na '.' (gap/wildcard)"""
    if gap_fraction <= 0:
        return seq
    seq = list(seq)
    n = max(1, int(len(seq) * gap_fraction))
    idx = random.sample(range(len(seq)), n)
    for i in idx:
        seq[i] = '.'
    return "".join(seq)

def generate_patterns(text, count, length, gap_fraction):
    """Generuje losowe motywy z genomu z dziurami"""
    patterns = []
    L = len(text)
    for _ in range(count):
        i = random.randint(0, L - length - 1)
        motif = text[i:i+length]
        motif = add_gaps(motif, gap_fraction)
        patterns.append(motif)
    return patterns

def save_patterns(patterns, filename):
    with open(filename, "w") as f:
        for p in patterns:
            f.write(p + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta", help="Plik FASTA z genomem")
    parser.add_argument("--gaps", type=float, default=0.2,
                        help="Ułamek dziur (0.0 – 0.9)")
    parser.add_argument("--prefix", default="patterns",
                        help="Prefix wynikowych plików")
    args = parser.parse_args()

    text = load_fasta(args.fasta)

    configs = [
        (10, 10),
        (50, 12),
        (200, 20)
    ]

    for count, length in configs:
        pats = generate_patterns(text, count, length, args.gaps)
        filename = f"{args.prefix}_{count}.txt"
        save_patterns(pats, filename)
        print(f"[OK] zapisano: {filename}")
