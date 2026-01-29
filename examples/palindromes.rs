use seqlib::sequence::DnaSeq;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let palindrome = DnaSeq::new("GAATTC")?;
    let ambiguous = DnaSeq::new("NNNNNN")?;

    assert!(palindrome.is_palindromic());
    assert!(!ambiguous.is_palindromic());

    println!("{palindrome} is palindromic");
    println!("{ambiguous} is NOT palindromic (ambiguous)");

    Ok(())
}
