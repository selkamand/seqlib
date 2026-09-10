use seqlib::coords::{BaseInterval, BasePos, Region, Strand};
use seqlib::mutations::{MutationWithContext, SmallMutation};
use seqlib::sequences::{BaseSliceExt, SourcedSeq};
use seqlib::{basepos, dna};
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    let mutation = SmallMutation::new(
        "Chr1".to_owned(),
        basepos!(2000),
        dna!("A"),
        dna!("C"),
        Some(Strand::Positive),
    );

    let interval = BaseInterval::try_new(basepos!(2000), basepos!(2000))?;

    let context = SourcedSeq::new(
        dna!("ACTGATCGAACGAGCATGCTACGGGGCCGATCGATTATCGATCAGTCA"),
        Region::new("Chr1", interval),
        Some(Strand::Positive),
    );

    let mutation_with_context = MutationWithContext::new(mutation, context)?;

    eprintln!("{mutation_with_context}");

    eprintln!(
        "-----Full Sequence Comparison----\n{}",
        mutation_with_context.to_difference_string()
    );

    let tnc = mutation_with_context.kmer_centered_on_anchor(3);
    eprintln!(
        "\n\nTNC: {}",
        tnc.map(|x| x.to_string_upper())
            .unwrap_or("Could not be found".to_string())
    );

    Ok(())
}
