fn main() {
    let seq = seqlib::sequences::DnaSeq::new("AGACT").unwrap();
    eprintln!("{} Sequence: {}", seq.alphabet(), seq);
    eprintln!("Length: {} nucleotides", seq.len());
    eprintln!(
        "Middle Base: {}",
        seq.middlebase().unwrap_or(&seqlib::base::DnaBase::N)
    );

    eprintln!("{}", seq.describe());

    // Complement Base
    eprintln!("{}", seq.complement().describe());

    // Reverse Complement
    eprintln!("{}", seq.reverse_complement().describe());
}
