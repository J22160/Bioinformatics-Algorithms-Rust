//! Small Parsimony — Fitch's Algorithm for Ancestral State Reconstruction

use std::collections::HashSet;

/// A node in the parsimony tree.
#[derive(Debug)]
pub enum ParsimonyNode {
    Leaf {
        label: String,
        character: char,
    },
    Internal {
        label: String,
        left: Box<ParsimonyNode>,
        right: Box<ParsimonyNode>,
        assigned: Option<char>,
    },
}

/// Fitch's algorithm: bottom-up pass to compute parsimony sets.
/// Returns (parsimony_score, set of possible characters at this node).
pub fn fitch_bottom_up(node: &ParsimonyNode) -> (usize, HashSet<char>) {
    match node {
        ParsimonyNode::Leaf { character, .. } => {
            let mut set = HashSet::new();
            set.insert(*character);
            (0, set)
        }
        ParsimonyNode::Internal { left, right, .. } => {
            let (l_score, l_set) = fitch_bottom_up(left);
            let (r_score, r_set) = fitch_bottom_up(right);

            let intersection: HashSet<char> = l_set.intersection(&r_set).copied().collect();

            if intersection.is_empty() {
                let union: HashSet<char> = l_set.union(&r_set).copied().collect();
                (l_score + r_score + 1, union)
            } else {
                (l_score + r_score, intersection)
            }
        }
    }
}

/// Fitch's algorithm: top-down pass to assign characters.
pub fn fitch_top_down(node: &mut ParsimonyNode, parent_char: Option<char>) {
    match node {
        ParsimonyNode::Leaf { .. } => {}
        ParsimonyNode::Internal { left, right, assigned, .. } => {
            let (_, set) = fitch_bottom_up_immutable(node);

            let chosen = if let Some(pc) = parent_char {
                if set.contains(&pc) { pc } else { *set.iter().next().unwrap() }
            } else {
                *set.iter().next().unwrap()
            };
            *assigned = Some(chosen);

            fitch_top_down(left, Some(chosen));
            fitch_top_down(right, Some(chosen));
        }
    }
}

fn fitch_bottom_up_immutable(node: &ParsimonyNode) -> (usize, HashSet<char>) {
    match node {
        ParsimonyNode::Leaf { character, .. } => {
            let mut s = HashSet::new(); s.insert(*character); (0, s)
        }
        ParsimonyNode::Internal { left, right, .. } => {
            let (ls, lset) = fitch_bottom_up_immutable(left);
            let (rs, rset) = fitch_bottom_up_immutable(right);
            let inter: HashSet<char> = lset.intersection(&rset).copied().collect();
            if inter.is_empty() {
                (ls + rs + 1, lset.union(&rset).copied().collect())
            } else {
                (ls + rs, inter)
            }
        }
    }
}

fn print_tree(node: &ParsimonyNode, prefix: &str, is_last: bool) {
    let connector = if is_last { "└── " } else { "├── " };
    match node {
        ParsimonyNode::Leaf { label, character, .. } => {
            println!("{}{}{} [{}]", prefix, connector, label, character);
        }
        ParsimonyNode::Internal { label, left, right, assigned, .. } => {
            let ch = assigned.map_or('?', |c| c);
            println!("{}{}{} [{}]", prefix, connector, label, ch);
            let new_prefix = format!("{}{}", prefix, if is_last { "    " } else { "│   " });
            print_tree(left, &new_prefix, false);
            print_tree(right, &new_prefix, true);
        }
    }
}

fn main() {
    println!("=== Fitch's Algorithm: Ancestral Nucleotide Reconstruction ===\n");

    // Build a small primate tree with nucleotides at one position
    let mut tree = ParsimonyNode::Internal {
        label: "Root".into(),
        left: Box::new(ParsimonyNode::Internal {
            label: "Apes".into(),
            left: Box::new(ParsimonyNode::Internal {
                label: "GreatApes".into(),
                left: Box::new(ParsimonyNode::Leaf { label: "Human".into(), character: 'A' }),
                right: Box::new(ParsimonyNode::Leaf { label: "Chimp".into(), character: 'A' }),
                assigned: None,
            }),
            right: Box::new(ParsimonyNode::Leaf { label: "Gorilla".into(), character: 'A' }),
            assigned: None,
        }),
        right: Box::new(ParsimonyNode::Internal {
            label: "OldWorld".into(),
            left: Box::new(ParsimonyNode::Leaf { label: "Macaque".into(), character: 'G' }),
            right: Box::new(ParsimonyNode::Leaf { label: "Baboon".into(), character: 'G' }),
            assigned: None,
        }),
        assigned: None,
    };

    let (score, root_set) = fitch_bottom_up(&tree);
    println!("Parsimony score: {}", score);
    println!("Possible root characters: {:?}\n", root_set);

    fitch_top_down(&mut tree, None);

    println!("Tree with ancestral assignments:");
    print_tree(&tree, "", true);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_no_changes() {
        let tree = ParsimonyNode::Internal {
            label: "R".into(),
            left: Box::new(ParsimonyNode::Leaf { label: "A".into(), character: 'T' }),
            right: Box::new(ParsimonyNode::Leaf { label: "B".into(), character: 'T' }),
            assigned: None,
        };
        let (score, _) = fitch_bottom_up(&tree);
        assert_eq!(score, 0);
    }

    #[test]
    fn test_one_change() {
        let tree = ParsimonyNode::Internal {
            label: "R".into(),
            left: Box::new(ParsimonyNode::Leaf { label: "A".into(), character: 'A' }),
            right: Box::new(ParsimonyNode::Leaf { label: "B".into(), character: 'G' }),
            assigned: None,
        };
        let (score, _) = fitch_bottom_up(&tree);
        assert_eq!(score, 1);
    }
}
