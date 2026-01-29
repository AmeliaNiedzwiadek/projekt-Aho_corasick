# patterns_generator.py
import random
import argparse


def load_fasta(path):
    """
    Wczytuje sekwencję z pliku FASTA i zwraca ją jako pojedynczy łańcuch.

    Linie nagłówkowe (rozpoczynające się od znaku '>') są pomijane.
    Wszystkie fragmenty sekwencji są łączone i konwertowane do wielkich liter.

    :param path: Ścieżka do pliku FASTA.
    :type path: str
    :return: Połączona sekwencja.
    :rtype: str
    """
    seqs = []
    with open(path) as f:
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


def add_gaps(seq, gap_fraction=0.0):
    """
    Wprowadza luki (wildcard '.') do sekwencji z określonym udziałem.

    Wybrane losowo pozycje są zastępowane znakiem '.'.

    :param seq: Sekwencja wejściowa.
    :type seq: str
    :param gap_fraction: Ułamek pozycji do zamiany na luki.
    :type gap_fraction: float
    :return: Sekwencja z wprowadzonymi lukami.
    :rtype: str
    """
    if gap_fraction <= 0:
        return seq

    seq = list(seq)
    n = max(1, int(len(seq) * gap_fraction))
    idx = random.sample(range(len(seq)), n)

    for i in idx:
        seq[i] = '.'

    return "".join(seq)


def generate_patterns(text, count, length, gap_fraction):
    """
    Generuje losowe motywy (patterny) z sekwencji referencyjnej.

    Każdy motyw jest wycinany z losowej pozycji i opcjonalnie
    modyfikowany przez wprowadzenie luk.

    :param text: Sekwencja referencyjna.
    :type text: str
    :param count: Liczba generowanych wzorców.
    :type count: int
    :param length: Długość pojedynczego wzorca.
    :type length: int
    :param gap_fraction: Ułamek pozycji zamienianych na luki.
    :type gap_fraction: float
    :return: Lista wygenerowanych wzorców.
    :rtype: list[str]
    """
    patterns = []
    L = len(text)

    for _ in range(count):
        i = random.randint(0, L - length - 1)
        motif = text[i:i+length]
        motif = add_gaps(motif, gap_fraction)
        patterns.append(motif)

    return patterns


def save_patterns(patterns, filename):
    """
    Zapisuje listę wzorców do pliku tekstowego.

    Każdy wzorzec zapisywany jest w osobnej linii.

    :param patterns: Lista wzorców.
    :type patterns: list[str]
    :param filename: Nazwa pliku wyjściowego.
    :type filename: str
    :return: None
    :rtype: None
    """
    with open(filename, "w") as f:
        for p in patterns:
            f.write(p + "\n")


if __name__ == "__main__":
    """
    Punkt wejścia programu w trybie wiersza poleceń.

    Program generuje zestawy plików ze wzorcami o różnych
    licznościach i długościach na podstawie sekwencji FASTA.
    """
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
