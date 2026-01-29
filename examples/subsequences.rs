use seqlib::sequence::DnaSeq;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let seq = DnaSeq::new("ACGTAC")?;

    let owned = seq.subseq(1, 4)?;
    let view = seq.subseq_slice(1, 4)?;

    println!("Original: {}", seq);
    println!("Owned subseq: {}", owned);
    println!("Slice length: {}", view.len());

    Ok(())
}
