use seqlib::{context, coord::Pos, mutation::SmallMutation, sequence::DnaSeq};

fn main() {
    let seq = DnaSeq::new("ACTGACGTA").unwrap();
    let mutation = SmallMutation::new(
        "chr1".to_owned(),
        Pos::new(50usize).unwrap(),
        DnaSeq::new("A").unwrap(),
        DnaSeq::new("C").unwrap(),
        false,
        false,
    );

    let context = context::ContextWindow::new();
    DnaSeq::new("ACTGATCGATCGAGCATGCTACGGGGCCGATCGATTATCGATCAGTCA")
}
