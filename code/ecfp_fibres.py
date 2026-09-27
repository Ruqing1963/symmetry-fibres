# -*- coding: utf-8 -*-
"""Fibres of circular (Morgan / ECFP) fingerprints, checked with RDKit.

(A) Explicit collision pairs.  Distinct molecules (different canonical SMILES, non-isomorphic graphs) whose
    Morgan fingerprints of radius r coincide: the hashed 2048-bit ECFP, the hashed count fingerprint and the
    unhashed multiset of environment identifiers.  Built by the transposition principle: a linear chain
    A-u-X-u-Y-u-Z-u-B with identical linkers u long enough that no radius-r environment reaches from one
    substituent to the next; exchanging X and Y preserves every environment.  Radius r+1 separates the pair
    when the linker has exactly 2r atoms.
(B) Substituted linear chains as strings.  Molecules CH3-(CH(R_1))-...-(CH(R_N))-CH3 with R_i in
    {H, F, Cl} are strings over a 3-letter alphabet.  The radius-r Morgan identifier of an interior atom
    depends on its (2r+1)-window read UP TO REVERSAL (the neighbours of an atom are unordered), and the whole
    fingerprint on the multiset of these windows plus the end windows.  We enumerate all 3^N strings, group
    them by the unhashed count fingerprint, and compare the partition with the one predicted from the
    reversal-symmetric window spectrum.  Two molecules are the same iff their strings are equal or mirror
    images, so classes are also counted modulo reversal.
(C) Stereochemistry.  Enantiomers and E/Z isomers under the default (achiral) fingerprint and with chirality.

Writes results/ecfp_fibres_output.txt and results/ecfp_pairs.json.
"""
import os, sys, json, itertools, collections, time
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator, AllChem

RDLogger.DisableLog("rdApp.*")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")


class Tee:
    def __init__(self, path):
        self.f, self.o = open(path, "w", encoding="utf-8"), sys.stdout

    def write(self, s):
        self.f.write(s); self.o.write(s)

    def flush(self):
        self.f.flush(); self.o.flush()


def gens(r, chiral=False):
    g = rdFingerprintGenerator.GetMorganGenerator(radius=r, fpSize=2048, includeChirality=chiral)
    return g


def fingerprints(mol, r, chiral=False):
    g = gens(r, chiral)
    bits = g.GetFingerprint(mol).ToBitString()
    counts = tuple(sorted(g.GetCountFingerprint(mol).GetNonzeroElements().items()))
    sparse = tuple(sorted(g.GetSparseCountFingerprint(mol).GetNonzeroElements().items()))
    return bits, counts, sparse


def canon(smiles):
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))


def chain(units, linker):
    """A-u-X-u-Y-u-Z-u-B as SMILES; units are SMILES fragments for one substituted carbon"""
    u = "C" * linker
    return "C" + u + u.join(units) + u + "C"


def compare(name, s1, s2, radii=(1, 2, 3, 4), chiral=False):
    m1, m2 = Chem.MolFromSmiles(s1), Chem.MolFromSmiles(s2)
    c1, c2 = Chem.MolToSmiles(m1), Chem.MolToSmiles(m2)
    same_mol = c1 == c2
    out = dict(name=name, smiles=[c1, c2], same_molecule=same_mol, atoms=m1.GetNumAtoms(), radius={})
    print(f"    {name}: {c1}  vs  {c2}   same molecule: {same_mol}")
    for r in radii:
        f1, f2 = fingerprints(m1, r, chiral), fingerprints(m2, r, chiral)
        eq = [a == b for a, b in zip(f1, f2)]
        out["radius"][r] = dict(bits2048=eq[0], counts2048=eq[1], unhashed=eq[2])
        print(f"        radius {r} (ECFP{2 * r}): 2048-bit {eq[0]}, counts {eq[1]}, unhashed identifier multiset {eq[2]}")
    return out


# ------------------------------------------------------------------ (B) chains as strings

COARSE = {"H": "H", "F": "X", "L": "X", "M": "M", "^": "^"}


def window_class(s, r):
    """the reversal-symmetric radius-r window multiset of the chain string s (letters M = terminal methyl,
    H = CH2, F = CHF, L = CHCl).  A Morgan identifier of radius r sees the full letter of every atom within
    distance r-1 and, at distance exactly r, only the atom's own invariants (degree, H count): CHF and CHCl
    are then both 'X'.  Absent atoms beyond the chain ends are '^'.  The neighbours of an atom are unordered,
    so a window is read up to reversal."""
    t = "^" * r + s + "^" * r
    ws = []
    for i in range(len(s)):
        w = t[i:i + 2 * r + 1]
        w = COARSE[w[0]] + w[1:-1] + COARSE[w[-1]] if r > 0 else w
        ws.append(min(w, w[::-1]))
    return tuple(sorted(ws))


def chain_smiles(s):
    sub = {"H": "C", "F": "C(F)", "L": "C(Cl)"}
    return "C" + "".join(sub[a] for a in s) + "C"


def part_B(N, radii):
    alphabet = "HFL"
    strings = ["".join(p) for p in itertools.product(alphabet, repeat=N)]
    mols = {s: Chem.MolFromSmiles(chain_smiles(s)) for s in strings}
    canon_of = {s: Chem.MolToSmiles(mols[s]) for s in strings}
    n_mol = len(set(canon_of.values()))
    print(f"\n=== (B) substituted chains CH3-(CHR)_{N}-CH3, R in {{H, F, Cl}}: {len(strings)} strings, "
          f"{n_mol} distinct molecules (strings up to reversal) ===")
    print("    predicted partition: multiset of radius-r windows read up to reversal, full letters within distance\n"
          "    r-1 and CHF = CHCl at distance r (an identifier sees a neighbour's substituent one round later)")
    res = {}
    for r in radii:
        g = gens(r)
        t = time.time()
        fp_class = collections.defaultdict(set)
        for s in strings:
            key = tuple(sorted(g.GetSparseCountFingerprint(mols[s]).GetNonzeroElements().items()))
            fp_class[key].add(s)
        # predicted partition: strings with the same reversal-symmetric window spectrum.  The terminal
        # methyl carbons contribute windows too: model them as letters 'H' with the chain extended by one
        # unit on each side, i.e. the molecule string is 'H' + s + 'H' inside end markers.
        pred_class = collections.defaultdict(set)
        for s in strings:
            pred_class[window_class("M" + s + "M", r)].add(s)
        P_fp = sorted(sorted(v) for v in fp_class.values())
        P_pred = sorted(sorted(v) for v in pred_class.values())
        equal = P_fp == P_pred
        # collisions between distinct molecules
        coll = [v for v in fp_class.values() if len({canon_of[s] for s in v}) > 1]
        n_coll_mol = sum(len({canon_of[s] for s in v}) for v in coll)
        largest = max(coll, key=lambda v: len({canon_of[s] for s in v})) if coll else set()
        print(f"    radius {r}: {len(fp_class)} fingerprint classes; partition equals the reversal-symmetric "
              f"{2 * r + 1}-window spectrum: {equal}; classes with >= 2 distinct molecules: {len(coll)}, covering "
              f"{n_coll_mol} of {n_mol} molecules; largest class {len({canon_of[s] for s in largest})} molecules "
              f"[{time.time() - t:.0f}s]")
        if largest:
            ex = sorted({canon_of[s] for s in largest})[:4]
            print("        e.g. " + " | ".join(ex))
        res[r] = dict(classes=len(fp_class), partition_equals_window_spectrum=equal, collision_classes=len(coll),
                      molecules_in_collisions=n_coll_mol, molecules=n_mol, largest_class=len({canon_of[s] for s in largest}))
    return res


def main():
    os.makedirs(RESULTS, exist_ok=True)
    sys.stdout = Tee(os.path.join(RESULTS, "ecfp_fibres_output.txt"))
    from rdkit import __version__ as rdv
    print(f"RDKit {rdv}")
    out = {"A": [], "B": {}, "C": []}

    print("\n=== (A) explicit collision pairs: A-u-X-u-Y-u-Z-u-B with X and Y exchanged ===")
    X, Y, Z = "C(C)", "C(Cl)", "C(O)"        # methyl, chloro, hydroxy substituted carbons
    for linker in (4, 5, 6):
        s1, s2 = chain([X, Y, Z], linker), chain([Y, X, Z], linker)
        out["A"].append(compare(f"linker (CH2)_{linker}, X=CH(CH3), Y=CHCl, Z=CH(OH)", s1, s2))
    # a ring example: two macrocycles with the same units in a different cyclic order (X Y Z W vs X Z Y W)
    W = "C(Br)"
    u = "C" * 4
    ring1 = "C1" + u + X + u + Y + u + Z + u + W + u[:-1] + "C1"
    ring2 = "C1" + u + X + u + Z + u + Y + u + W + u[:-1] + "C1"
    out["A"].append(compare("macrocycle, cyclic order X Y Z W vs X Z Y W, linkers (CH2)_4", ring1, ring2))
    # pharmacophore geometry of the linker-4 pair: topological distances between the heteroatoms
    print("    topological (shortest-path) distances between the substituent heteroatoms, linker (CH2)_4:")
    for label, smi in (("X Y Z", chain([X, Y, Z], 4)), ("Y X Z", chain([Y, X, Z], 4))):
        m = Chem.MolFromSmiles(smi)
        D = Chem.GetDistanceMatrix(m)
        idx = {a.GetSymbol(): a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() in ("Cl", "O")}
        me = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "C" and a.GetDegree() == 1
              and any(n.GetDegree() == 3 for n in a.GetNeighbors())][0]      # the methyl branch
        print(f"        order {label}: d(CH3, Cl) = {int(D[me][idx['Cl']])}, d(Cl, OH) = {int(D[idx['Cl']][idx['O']])}, "
              f"d(CH3, OH) = {int(D[me][idx['O']])} bonds")

    out["B"] = part_B(9, (1, 2, 3))

    print("\n=== (C) stereoisomers ===")
    out["C"].append(compare("L-alanine vs D-alanine, achiral fingerprint", "C[C@H](N)C(=O)O", "C[C@@H](N)C(=O)O", (2,)))
    out["C"].append(compare("L-alanine vs D-alanine, chiral fingerprint", "C[C@H](N)C(=O)O", "C[C@@H](N)C(=O)O", (2,), chiral=True))
    out["C"].append(compare("(E)- vs (Z)-2-butene, achiral fingerprint", "C/C=C/C", "C/C=C\\C", (2,)))
    out["C"].append(compare("(E)- vs (Z)-2-butene, chiral fingerprint", "C/C=C/C", "C/C=C\\C", (2,), chiral=True))
    json.dump(out, open(os.path.join(RESULTS, "ecfp_pairs.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
