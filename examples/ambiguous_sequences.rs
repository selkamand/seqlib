use seqlib::sequence::DnaSeq;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let clean = DnaSeq::new("ACGT")?;
    let ambiguous = DnaSeq::new("ACNT")?;

    assert!(clean.all_unambiguous());
    assert!(ambiguous.any_ambiguous());

    println!("Clean: {} (unambiguous)", clean);
    println!("Ambiguous: {} (contains N)", ambiguous);

    Ok(())
}
