# seqlib

`seqlib` is a small, type-safe Rust library for working with DNA and RNA sequences.

It provides a **robust core representation** for biological sequences with:
- explicit DNA / RNA alphabets
- full IUPAC ambiguity support
- compiler-enforced correctness for operations like complementation

`seqlib` is designed as a **library dependency**, not a command-line tool.  
It is developed to support the `scarscape` project, but is usable on its own for
any Rust code that needs reliable nucleotide sequence handling.

---

## Design goals

- **Correctness first**: invalid sequences are rejected at construction time
- **Type safety**: DNA and RNA are distinct types, not runtime flags
- **Strict alphabets**: DNA and RNA reject invalid bases at construction time (e.g. `U` in DNA, `T` in RNA)
- **Explicit ambiguity handling**: ambiguity is modeled, not ignored
- **Ergonomic by default**: core sequence operations use copy-on-modify semantics, making them easy to compose and safe for downstream use
- **Explicit performance opt-ins**: in-place mutation methods are provided for performance- or memory-critical workflows where mutation is acceptable
- **Small surface area**: no I/O, no parsing frameworks, no CLI
- **Composable**: intended to be embedded in larger bioinformatics pipelines

---

## Adding to your project

`seqlib` is currently consumed directly from GitHub:

```bash
cargo add seqlib --git https://github.com/selkamand/seqlib
```

---

## Core types

Most users should work with the concrete sequence types:

- `DnaSeq` — validated DNA sequences (`A, C, G, T` + IUPAC ambiguity codes)
- `RnaSeq` — validated RNA sequences (`A, C, G, U` + IUPAC ambiguity codes)

These are type aliases over a generic `Seq<B>` implementation and are the
**recommended entry points** for sequence creation and manipulation.

The generic `Seq<B>` type exists to support reusable algorithms and future
extensions, but most downstream code should not need to interact with it
directly.

---

## Basic usage

```rust
use seqlib::sequences::DnaSeq;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let seq = DnaSeq::new("ACGTN")?;

    println!("Sequence: {}", seq);
    println!("Length: {}", seq.len());
    println!("Reverse complement: {}", seq.reverse_complement());

    Ok(())
}
```

Invalid input is rejected immediately:

```rust
let bad = DnaSeq::new("ACGTX"); // returns Err(...)
```

---

## Copy-on-modify vs in-place mutation

`seqlib` is designed to be **ergonomic and safe by default**.

All standard sequence operations—such as reverse complementation, subsequence
extraction, and reversal—use **copy-on-modify semantics**. These methods return
new sequences rather than mutating existing ones, making them easy to compose,
store, and pass through downstream code without surprising side effects.

For performance-critical or memory-sensitive workflows, `seqlib` also exposes
explicit **in-place mutation** methods (e.g. `reverse_complement_in_place`,
`subseq_in_place`). These methods are clearly named and opt-in, allowing callers
to trade ergonomics for efficiency when appropriate.

---

## Ambiguity handling

`seqlib` supports IUPAC ambiguity codes (e.g. `N`, `R`, `Y`, `S`), but **does not
treat ambiguous symbols as wildcards** in higher-level analyses.

Ambiguous bases represent *unknown concrete bases*, not flexible matches. As a
result, operations that require certainty (such as palindromic sequence
detection) will conservatively return `false` if ambiguity is present.

---

## Palindromic sequences

`seqlib` defines palindromic sequences using a **strict, certainty-based
definition**.

A sequence is considered palindromic only if:
- it contains **no ambiguous bases**,
- its length is **non-zero and even**, and
- it is identical to its **reverse complement** at the level of concrete bases.

This guarantees that palindromic sequences can be **counted and filtered without
overcounting** symbolically palindromic but ambiguous inputs such as `NNNNNN`.

```rust
use seqlib::sequences::DnaSeq;

assert!(DnaSeq::new("GAATTC").unwrap().is_palindromic());
assert!(!DnaSeq::new("NNNNNN").unwrap().is_palindromic());
assert!(!DnaSeq::new("AAA").unwrap().is_palindromic());
```

This conservative definition is intentional and designed for statistical and
motif-based analyses where false positives must be avoided.

---

## Features

- DNA and RNA alphabets with full IUPAC ambiguity support
- Infallible complement and reverse-complement operations
- Explicit ambiguity detection
- Allocation-aware APIs with both copy-on-modify and in-place variants
- No external dependencies beyond the Rust standard library

---

## Examples

See the [`examples/`](examples) directory for typical usage patterns and
small self-contained demonstrations.

---

## Out of Scope

To keep seqlib focused the following features will not be implemented in this
library.

- No FASTA / FASTQ parsing
- No file I/O
- No CLI
- No alignment, variant calling, or translation logic
- No soft-masking or case preservation

`seqlib` is intended to be a **foundation**, not a full bioinformatics toolkit.
