use seqlib::sequences::DnaSeq;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let eco_ri = DnaSeq::new("GAATTC")?;
    let ambiguous = DnaSeq::new("NNNNNN")?;

    assert!(eco_ri.is_palindromic());
    assert!(!ambiguous.is_palindromic());

    println!("{} is palindromic", eco_ri);
    println!("{} is NOT palindromic (ambiguous)", ambiguous);

    Ok(())
}
