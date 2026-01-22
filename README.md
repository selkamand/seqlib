
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
- **Explicit ambiguity handling**: ambiguity is modeled, not ignored
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

## Features

- DNA and RNA alphabets with full IUPAC ambiguity support
- Infallible complement and reverse-complement operations
- Explicit ambiguity detection
- Lightweight, allocation-aware design
- No external dependencies beyond the Rust standard library

---

## Examples

See the [`examples/`](examples) directory for typical usage patterns and
small self-contained demonstrations.

---

## Out of Scope

To keep seqlib focused the following features will not be implemented in this library.

- No FASTA / FASTQ parsing
- No file I/O
- No CLI
- No alignment, variant calling, or translation logic

`seqlib` is intended to be a **foundation**, not a full bioinformatics toolkit.
