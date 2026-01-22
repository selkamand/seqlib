use seqlib::sequences::DnaSeq;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let original = DnaSeq::new("ACGT")?;

    // Copy-on-modify (default)
    let rc = original.reverse_complement();
    println!("Original: {}", original);
    println!("Reverse complement: {}", rc);

    // In-place mutation (opt-in)
    let mut seq = original.clone();
    seq.reverse_complement_in_place();
    println!("Mutated in place: {}", seq);

    Ok(())
}
