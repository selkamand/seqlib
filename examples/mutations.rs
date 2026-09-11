use seqlib::{
    basepos,
    coords::{BaseInterval, BasePos, Region, Strand},
    dna,
    mutations::{DnaSmallMutation, MutationWithContext},
    sequences::SourcedSeq,
};

fn main() {
    // Define a small dna mutation (e.g. a snv / insertion / deletion / etc)
    let mutation = DnaSmallMutation::new(
        "chr1".to_owned(),
        basepos!(2004),
        dna!("A"),
        dna!("G"),
        Some(Strand::Positive),
    );

    // Print the DNA mutation
    println!("-------------------------");
    println!("Small DNA mutation:");
    println!("-------------------------");
    println!("{}", mutation.format_with_colour());

    // Define a sequence from a reference genome which contains the mutated site
    let context = SourcedSeq::new(
        dna!("ACGTACGTGCA"),
        Region::new(
            "chr1",
            BaseInterval::new(basepos!(2000), basepos!(2010))
                .expect("example context coordinates are valid"),
        ),
        Some(Strand::Positive),
    );

    // Print the context sequence
    println!("\n-------------------------");
    println!("Sequence from reference genome contextualising the mutation:");
    println!("-------------------------");
    println!("{}:", context.format_with_colour());

    // Create a dedicated MutationWithContext type
    let mutation_with_context = MutationWithContext::new(mutation, context.clone()).unwrap();
    println!("\n-------------------------");
    println!("Mutation with Context");
    println!("-------------------------");
    println!("{}", mutation_with_context.format_with_colour());

    // Representing Indels
    let indel = DnaSmallMutation::new(
        "chr1".to_owned(),
        basepos!(2004),
        dna!("A"),
        dna!(""),
        Some(Strand::Positive),
    );

    let indelcontext = context.clone();
    let indel_with_context = MutationWithContext::new(indel, indelcontext).unwrap();
    println!("\n-------------------------");
    println!("Indel with Context");
    println!("-------------------------");
    println!("{}", indel_with_context.format_with_colour());
}
