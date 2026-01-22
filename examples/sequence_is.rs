use seqlib::{errors::SeqError, sequences::DnaSeq};

fn main() -> Result<(), SeqError> {
    let seq = DnaSeq::new("ACCTAGGT")?;
    eprintln!("{}", seq.describe());

    Ok(())
}
