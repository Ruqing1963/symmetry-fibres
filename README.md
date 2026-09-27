# Symmetry-constrained spectral fibres

**Reverse-complement group actions on double-stranded DNA and the fibres of circular molecular fingerprints**

Zhengyi Chen<sup>1</sup> and Ruqing Chen<sup>2</sup>

<sup>1</sup> Guangxi Key Laboratory of Drug Discovery and Optimization, Guangxi Engineering Research Center for Pharmaceutical Molecular Screening and Druggability Evaluation, School of Pharmacy, Guilin Medical University, Guilin 541199, P. R. China (chenzhengyi@glmc.edu.cn)
<sup>2</sup> GUT Geoservice Inc., Montreal, Québec, Canada (ruqing@hotmail.com)

DOI: to be assigned (Zenodo)

## About

A representation that depends on an object only through its local content is constant on the object's *spectral fibre*. This paper treats two settings in which a symmetry enters.

- **Double-stranded DNA.** Reverse complementation generates a group acting on sequence space. A strand-agnostic local encoder is blind to the double-stranded fibre DS(S, k), which contains the one-strand fibres of both strands and is a lattice sum of BEST determinants, #P-hard to evaluate (results of the monograph *Algebraic Genomics*). For φX174 at k = 10, |DS| = 31,925,753,246,212: between the single-strand count 8,610,708,632 and the 1.4·10³⁶ Eulerian circuits of the double cover, which a naive count would give.
- **Circular fingerprints.** On substituted linear chains the radius-r Morgan (ECFP<sub>2r</sub>) partition is exactly the spectrum of (2r+1)-windows read up to reversal, with the substituent at distance r coarsened to its atom type (Theorem 5). Verified exhaustively on all 19,683 chains CH₃(CHR)₉CH₃, R ∈ {H, F, Cl}, at radii 1–3. Among the 9,963 distinct molecules, 3.2% share their ECFP4 with another molecule; none share ECFP6.
- **Collision pairs.** Distinct molecules with identical ECFP4 and ECFP6 (2048-bit, count and unhashed), built by exchanging substituents across repeated (CH₂)<sub>ℓ</sub> linkers; the fingerprints coincide iff ℓ ≥ 2r−1, and the bound is attained. In the linker-4 pair the Cl-to-OH distance is 7 bonds in one molecule and 12 in the other.
- **A unified quotient** covers composition encoders, convolutional encoders, strand-agnostic DNA encoders and circular fingerprints.

Companion papers: [10.5281/zenodo.22950095](https://doi.org/10.5281/zenodo.22950095), [10.5281/zenodo.22959211](https://doi.org/10.5281/zenodo.22959211), [10.5281/zenodo.22960439](https://doi.org/10.5281/zenodo.22960439), [10.5281/zenodo.22963882](https://doi.org/10.5281/zenodo.22963882), [10.5281/zenodo.22966020](https://doi.org/10.5281/zenodo.22966020). Monograph: [10.5281/zenodo.22943487](https://doi.org/10.5281/zenodo.22943487).

## Repository layout

```
symmetry-fibres/
├── paper/      symmetry_fibres.tex, symmetry_fibres.pdf     the paper
├── code/       ecfp_fibres.py        collision pairs, exhaustive chain check, census, stereo examples (RDKit)
├── data/       phiX174.fasta         bacteriophage phiX174, RefSeq NC_001422.1
├── results/    ecfp_fibres_output.txt, ecfp_pairs.json
└── README.md, CITATION.cff, requirements.txt, LICENSE (code, MIT), LICENSE-CC-BY-4.0.md (text)
```

The φX174 numbers are computed by the verification scripts of the monograph's repository, <https://github.com/Ruqing1963/algebraic-genomics> (`papers/Ch07_*/ds_lattice.py`, `ds_flips.py`; `papers/Ch09_*/ds_fermion.py`, `ds_count.py`), whose transcripts the paper quotes. They are not recomputed here.

## Reproducing

```
cd code
python -u ecfp_fibres.py     # a few seconds; needs rdkit
cd ../paper
pdflatex symmetry_fibres
pdflatex symmetry_fibres
```

Tested with Python 3.12 and RDKit 2026.03.6.

## License

Text and figures: [CC BY 4.0](LICENSE-CC-BY-4.0.md). Code: [MIT](LICENSE). The φX174 sequence is from RefSeq (public domain).
