# mutations.py
def compare_sequences(seq1, seq2):
    """
    Porównuje dwie sekwencje DNA i wykrywa:
    - mutacje punktowe (SNP)
    - insercje
    - delecje

    Zwraca listę opisów zmian.
    """
    changes = []
    i = j = 0
    L1, L2 = len(seq1), len(seq2)

    while i < L1 and j < L2:
        if seq1[i] == seq2[j]:
            i += 1
            j += 1
        else:
            # SNP
            if i+1 < L1 and j+1 < L2 and seq1[i+1] == seq2[j+1]:
                changes.append(f"SNP at pos {i}: {seq1[i]} → {seq2[j]}")
                i += 1
                j += 1
            # deletion
            elif i+1 < L1 and seq1[i+1] == seq2[j]:
                changes.append(f"Deletion at pos {i}: deleted {seq1[i]}")
                i += 1
            # insertion
            elif j+1 < L2 and seq1[i] == seq2[j+1]:
                changes.append(f"Insertion at pos {i}: inserted {seq2[j]}")
                j += 1
            else:
                changes.append(f"Complex mutation around pos {i}/{j}")
                i += 1
                j += 1

    # jeśli jedna sekwencja jest dłuższa
    while i < L1:
        changes.append(f"Deletion at end: {seq1[i]}")
        i += 1

    while j < L2:
        changes.append(f"Insertion at end: {seq2[j]}")
        j += 1

    return changes


# prosta funkcja testowa
if __name__ == "__main__":
    seqA = "ATGCCGTA"
    seqB = "ATGACGGA"

    diff = compare_sequences(seqA, seqB)
    print("Differences:")
    for d in diff:
        print(" -", d)

