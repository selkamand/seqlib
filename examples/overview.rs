use seqlib::sequences::{DnaSeq, RnaSeq};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // --- DNA ---
    // Non-palindromic, odd length → shows middle base handling
    let dna = DnaSeq::new("AGACT")?;
    println!("DNA\n{}", dna.describe());

    // Strict DNA: U is rejected
    assert!(DnaSeq::new("ACGU").is_err());

    // --- RNA ---
    let rna = RnaSeq::new("ACGU")?;
    println!("RNA\n{}", rna.describe());

    // Strict RNA: T is rejected
    assert!(RnaSeq::new("ACGT").is_err());

    Ok(())
}
