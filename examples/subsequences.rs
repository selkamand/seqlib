use seqlib::{
    basepos,
    coords::{BaseInterval, BasePos},
    sequences::{BaseSliceExt, IupacDnaSeq},
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Define a sequence
    let seq = IupacDnaSeq::new("ACGTAC")?;
    println!("{seq} <- Sequence (original)");

    // Define a interval (1 based start & end, both-end inclusive)
    let start = basepos!(2);
    let end = basepos!(4);
    let interval = BaseInterval::new(start, end)?;

    // Highlight where this interval is on our sequence
    println!(
        "{} <- Sequence (annotated by interval {} [{}bp])",
        seq.format_with_highlight_interval(Some(&interval)),
        interval,
        interval.len()
    );

    // Grab the subsequence (owned copy)
    let subseq = seq.subseq(&interval)?;

    // Print out the slice with
    println!("{subseq} <- sub-sequence");

    // If you just want to borrow a slice, use the subseq_slice method
    let subseq_slice = seq.subseq_slice(&interval)?;
    println!("{} <- sub-sequence", subseq_slice.to_string_upper());

    Ok(())
}
