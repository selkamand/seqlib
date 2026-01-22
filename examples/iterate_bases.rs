use seqlib::sequences::DnaSeq;
use seqlib::base::Base;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let seq = DnaSeq::new("ACGT")?;

    for base in seq.as_slice() {
        println!("{}", base.to_char());
    }

    Ok(())
}
