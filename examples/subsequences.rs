use seqlib::{
    coords::{BaseInterval, BasePos, InterbaseInterval, InterbasePos},
    sequences::{BaseSliceExt, DnaSeq},
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Define a sequence
    println!("\n-------------------------");
    println!("Define Sequence:");
    println!("-------------------------");
    let seq = DnaSeq::new("ACGTAC")?;
    println!("{seq} <- Sequence (original)");

    println!("\n-------------------------");
    println!("Construct an Interval [2-4]:");
    println!("-------------------------");

    // Define either a BaseInterval (1 based start & end, both-end inclusive)
    // or a InterbaseInterval (0 based start half open)
    let start = BasePos::new(2)?;
    let end = BasePos::new(4)?;
    let interval = BaseInterval::new(start, end)?;

    // Highlight where this interval is on our sequence
    println!(
        "{} <- Sequence (annotated by interval {} [{}bp])",
        seq.format_with_highlighted_base_interval(Some(&interval)),
        interval,
        interval.len()
    );

    println!("\n-------------------------");
    println!("Subsequence (owned copy):");
    println!("-------------------------");
    // Grab the subsequence (owned copy)
    let subseq = seq.subseq_by_base_interval(&interval)?;

    // Print out the slice
    println!("{subseq} <- sub-sequence");

    // If you just want to borrow a slice, use the slice_by.. methods
    println!("\n-------------------------");
    println!("Subsequence (borrow a slice):");
    println!("-------------------------");
    let subseq_slice = seq.slice_by_base_interval(&interval)?;
    println!("{} <- sub-sequence", subseq_slice.to_string_upper());

    println!("\n-------------------------");
    println!("All slicing/subsequences also work \nwith interbase intervalse:");
    println!("-------------------------");
    let start2 = InterbasePos::new(2);
    let end2 = InterbasePos::new(4);
    let interval2 = InterbaseInterval::new(start2, end2)?;

    // Highlight where this interval is on our sequence
    println!(
        "{} <- Sequence (annotated by interval {} [{}bp])",
        seq.format_with_highlighted_interbase_interval(Some(&interval2)),
        interval2,
        interval2.len()
    );

    Ok(())
}
