use seqlib::{
    coord::Pos,
    sequence::{DnaSeq, Seq},
};

fn main() {
    let seq = DnaSeq::new("ACTGATTTT").unwrap();

    // Highlight the 4th element in sequence vector (5th base in sequence since rust vectors are
    // 0-based)
    println!("{}", seq.format_with_highlight_index(Some(4)));

    // Highlight the 4th base by specifying its 1-based position
    let position = Pos::new(4usize).unwrap();
    println!("{}", seq.format_with_highlight_pos(Some(position)));
}
